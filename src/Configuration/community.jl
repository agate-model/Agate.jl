"""Build a plankton community from a base community spec.

This helper updates the size structure for one or more plankton groups,
leaving all other group fields intact. Size-class counts are inferred from
explicit diameter vectors or read from structured size definitions.

The ordering of groups in the returned `NamedTuple` is the ordering of
`base` (i.e. `keys(base)`).

Keywords
--------
- `diameters=(; ...)`: optional per-group size-structure inputs, keyed by group symbol.

Accepted size-structure inputs
------------------------------
Each group size structure must define class diameters either as a generated
range or as explicit values. Accepted public forms are:

- `(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)` for a generated range.
- `[1.0, 3.2, 10.0]` for explicit class diameters.

Examples
--------
```julia
community = build_plankton_community(base;
    diameters=(grazers=grazer_sizes, microbes=microbe_sizes),
)
```
"""
function build_plankton_community(base::NamedTuple; diameters::NamedTuple=NamedTuple())
    base_keys = keys(base)
    for k in keys(diameters)
        k in base_keys || throw(ArgumentError("diameters: unknown group symbol $k"))
    end

    values_ = ntuple(length(base_keys)) do i
        g = base_keys[i]
        spec = getfield(base, g)
        diameter_definition = if hasproperty(diameters, g)
            getfield(diameters, g)
        else
            getproperty(spec, :diameters)
        end
        normalized_diameters = normalize_diameters(diameter_definition)
        size_class_count = if !isnothing(normalized_diameters.n)
            normalized_diameters.n
        elseif hasproperty(spec, :n)
            getproperty(spec, :n)
        else
            throw(ArgumentError("group $g: size structure must provide `n`"))
        end

        return (; spec..., n=size_class_count, diameters=diameter_definition)
    end

    return NamedTuple{base_keys}(values_)
end

"""Abstract supertype for diameter specifications."""
abstract type AbstractDiameterSpecification end

"""A diameter specification defined by an explicit list of class diameters.

The size-class count is inferred from `length(diameters)`.
"""
struct DiameterListSpecification{T,VT<:AbstractVector{T}} <: AbstractDiameterSpecification
    diameters::VT
end

"""A diameter specification defined by a class count, range, and splitting method.

`n` is the number of size classes. `min_diameter` and `max_diameter` define the
range of equivalent spherical diameters. `splitting` selects the spacing method,
for example `:log_splitting` or `:linear_splitting`.
"""
struct DiameterRangeSpecification{I<:Integer,T1,T2} <: AbstractDiameterSpecification
    n::I
    min_diameter::T1
    max_diameter::T2
    splitting::Symbol
end

"""Construction-time representation of a parsed plankton community.

`CommunityContext` stores the flattened plankton layout, class metadata, and role
axes needed while resolving parameters and interaction matrices. It is setup-time
semantic state and is not passed into tracer-tendency kernels.

Fields
------
- `T`: scalar type used for construction.
- `n_total`: total number of plankton classes.
- `diameters`: flattened diameter vector in global plankton order.
- `pfts`: per-class PFT specifications.
- `class_symbols`: flattened ecological class symbols such as `:P_1`, `:P_2`, or `:diat_1`.
- `group_symbols`: group symbol for each flattened class.
- `group_indices`: mapping from group symbol to flattened class indices.
- `consumer_indices`: flattened indices used for the consumer axis.
- `prey_indices`: flattened indices used for the prey axis.
"""
struct CommunityContext{T<:Real,VT<:AbstractVector{T}}
    scalar_type::Type{T}
    n_total::Int
    diameters::VT
    pfts::Vector{PFTSpecification}
    class_symbols::Vector{Symbol}
    group_symbols::Vector{Symbol}
    group_indices::Dict{Symbol,Vector{Int}}
    consumer_indices::Vector{Int}
    prey_indices::Vector{Int}
end

"""Return normalized diameter input and any size-class count it defines."""
function normalize_diameters(diameters::AbstractVector)
    (; n=length(diameters), specification=DiameterListSpecification(diameters))
end

function normalize_diameters(spec::NamedTuple)
    required = (:n, :min_esd, :max_esd, :splitting)
    all(hasproperty(spec, field) for field in required) || throw(
        ArgumentError(
            "diameter range NamedTuple must define `n`, `min_esd`, `max_esd`, and `splitting`",
        ),
    )
    return (;
        n=spec.n,
        specification=DiameterRangeSpecification(
            spec.n, spec.min_esd, spec.max_esd, spec.splitting
        ),
    )
end

function normalize_diameters(spec::DiameterListSpecification)
    (; n=length(spec.diameters), specification=spec)
end

normalize_diameters(spec::DiameterRangeSpecification) = (; n=spec.n, specification=spec)

normalize_diameters(spec) = throw(ArgumentError("invalid `diameters` specification"))

"""Validate a population-community realization.

Throws a single `ArgumentError` listing all structural issues.
"""
function validate_community(community)
    community isa NamedTuple || throw(ArgumentError("community must be a NamedTuple"))
    issues = String[]

    for k in keys(community)
        spec = getfield(community, k)

        if !hasproperty(spec, :diameters)
            push!(issues, "group $(k): missing required field `diameters`")
        else
            d = getproperty(spec, :diameters)
            normalized_diameters = try
                normalize_diameters(d)
            catch err
                err isa ArgumentError || rethrow()
                push!(issues, "group $(k): invalid `diameters` specification")
                nothing
            end

            if !isnothing(normalized_diameters)
                spec_n = hasproperty(spec, :n) ? getproperty(spec, :n) : nothing
                diameter_n = normalized_diameters.n

                if !isnothing(spec_n) && !isnothing(diameter_n) && spec_n != diameter_n
                    push!(
                        issues,
                        "group $(k): `n` ($(spec_n)) does not match diameters.n ($(diameter_n))",
                    )
                end

                n_source = !isnothing(diameter_n) ? diameter_n : spec_n
                if isnothing(n_source)
                    push!(issues, "group $(k): missing required field `n` for non-explicit diameters")
                elseif !(n_source isa Integer) || n_source < 1
                    push!(issues, "group $(k): `n` must be a positive integer")
                elseif d isa AbstractVector && n_source != length(d)
                    push!(
                        issues,
                        "group $(k): `n` ($(n_source)) does not match length(diameters) ($(length(d)))",
                    )
                end
            end
        end

        if hasproperty(spec, :class_symbols)
            symbols = getproperty(spec, :class_symbols)
            symbols isa Tuple || symbols isa AbstractVector ||
                push!(issues, "group $(k): `class_symbols` must be a tuple or vector")
            if symbols isa Tuple || symbols isa AbstractVector
                all(symbol -> symbol isa Symbol, symbols) ||
                    push!(issues, "group $(k): `class_symbols` must contain only Symbols")
            end
        end

        if !hasproperty(spec, :pft)
            push!(issues, "group $(k): missing required field `pft`")
        else
            pft = getproperty(spec, :pft)
            pft isa Union{PFTSpecification,NamedTuple} ||
                push!(issues, "group $(k): `pft` must be PFTSpecification or NamedTuple")
        end
    end

    isempty(issues) || throw(ArgumentError(join(issues, "\n")))
    return nothing
end

"""Parse `community` into a flattened `CommunityContext`.

Keyword arguments
-----------------
- `biogeochem_tracers`: non-plankton tracer names used to reject class-name collisions.
- `interaction_roles`: optional `NamedTuple` with fields `consumers` and
  `prey`. Each field may be `nothing`, a collection of group symbols, an index
  vector, or a boolean mask.
When `interaction_roles` is omitted, both interaction axes include all classes.
Parameter applicability is resolved separately from these interaction axes and derives directly from process participation.
"""
function parse_community(
    ::Type{T},
    community::NamedTuple;
    biogeochem_tracers::Tuple=(),
    interaction_roles=nothing,
) where {T<:Real}
    if !isnothing(interaction_roles)
        (
            hasproperty(interaction_roles, :consumers) &&
            hasproperty(interaction_roles, :prey)
        ) || throw(
            ArgumentError(
                "interaction_roles must have fields :consumers and :prey (each may be `nothing`, group Symbols, indices, or boolean masks).",
            ),
        )
    end

    group_order = keys(community)
    class_symbols = Symbol[]
    group_of = Symbol[]
    pfts = PFTSpecification[]
    diameters = T[]

    for g in group_order
        spec = getfield(community, g)
        normalized_diameters = normalize_diameters(getproperty(spec, :diameters))
        dspec = normalized_diameters.specification
        n = if dspec isa DiameterListSpecification
            length(dspec.diameters)
        elseif !isnothing(normalized_diameters.n)
            normalized_diameters.n
        elseif hasproperty(spec, :n)
            getproperty(spec, :n)
        else
            throw(ArgumentError("group $g: missing required field `n`"))
        end
        ds = param_compute_diameters(T, n, dspec)
        pft_raw = getproperty(spec, :pft)
        pft = pft_raw isa PFTSpecification ? pft_raw : PFTSpecification(pft_raw)

        realized_class_symbols = if hasproperty(spec, :class_symbols)
            supplied = Tuple(getproperty(spec, :class_symbols))
            length(supplied) == n || throw(
                ArgumentError("group $g: class_symbols must have length $n")
            )
            all(symbol -> symbol isa Symbol, supplied) || throw(
                ArgumentError("group $g: class_symbols must contain only Symbols")
            )
            supplied
        else
            ntuple(i -> Symbol(string(g), "_", i), n)
        end

        for i in 1:n
            push!(class_symbols, realized_class_symbols[i])
            push!(group_of, g)
            push!(pfts, pft)
            push!(diameters, ds[i])
        end
    end

    biogeochem_symbols = Set(biogeochem_tracers)
    length(unique(class_symbols)) == length(class_symbols) || throw(
        ArgumentError("community realizes duplicate ecological class symbols")
    )
    conflicting_symbols = [symbol for symbol in class_symbols if symbol in biogeochem_symbols]
    isempty(conflicting_symbols) || throw(
        ArgumentError(
            "population class names conflict with non-plankton tracers: $(unique(conflicting_symbols))"
        ),
    )

    n_total = length(class_symbols)

    group_indices = Dict{Symbol,Vector{Int}}()
    for (i, g) in enumerate(group_of)
        push!(get!(group_indices, g, Int[]), i)
    end

    interaction_roles_resolved =
        isnothing(interaction_roles) ? (consumers=nothing, prey=nothing) : interaction_roles

    function indices_for_role(role, role_name::Symbol)
        if role === nothing
            return collect(1:n_total)
        elseif role isa AbstractVector{Bool}
            length(role) == n_total ||
                throw(ArgumentError("$role_name mask must have length $n_total"))
            return findall(role)
        elseif role isa AbstractVector{Int}
            idx = collect(role)
            all(1 .<= idx .<= n_total) ||
                throw(ArgumentError("$role_name indices must be in 1:$n_total"))
            return idx
        elseif role isa Tuple || role isa AbstractVector{Symbol}
            requested = Set{Symbol}(role)
            for g in requested
                haskey(group_indices, g) ||
                    throw(ArgumentError("Unknown group symbol $g in $role_name roles"))
            end

            idx = Int[]
            for g in group_order
                g ∈ requested || continue
                append!(idx, group_indices[g])
            end
            return idx
        else
            throw(
                ArgumentError(
                    "$role_name roles must be nothing, a Bool mask, an Int index vector, or a collection of group Symbols",
                ),
            )
        end
    end

    consumer_indices = indices_for_role(
        getproperty(interaction_roles_resolved, :consumers), :consumers
    )
    prey_indices = indices_for_role(getproperty(interaction_roles_resolved, :prey), :prey)

    community_context = CommunityContext{T,typeof(diameters)}(
        T,
        n_total,
        diameters,
        pfts,
        class_symbols,
        group_of,
        group_indices,
        consumer_indices,
        prey_indices,
    )

    return community_context
end

@inline function param_check_length(name::Symbol, expected::Int, got::Int)
    if expected != got
        throw(ArgumentError("$(name) must have length $(expected) but has length $(got)"))
    end
    return nothing
end

function param_compute_diameters(
    ::Type{T}, n::Int, spec::DiameterRangeSpecification
) where {T<:Real}
    min_d = T(spec.min_diameter)
    max_d = T(spec.max_diameter)

    if n == 1
        return T[min_d]
    end

    diameters = Vector{T}(undef, n)

    if spec.splitting === :log_splitting
        log_min = log(min_d)
        log_max = log(max_d)
        step = (log_max - log_min) / T(n - 1)
        @inbounds for i in 1:n
            diameters[i] = exp(log_min + T(i - 1) * step)
        end
    elseif spec.splitting === :linear_splitting
        step = (max_d - min_d) / T(n - 1)
        @inbounds for i in 1:n
            diameters[i] = min_d + T(i - 1) * step
        end
    else
        throw(ArgumentError("Unsupported splitting method: $(spec.splitting)"))
    end

    return diameters
end

function param_compute_diameters(
    ::Type{T}, n::Int, spec::DiameterListSpecification
) where {T<:Real}
    param_check_length(:diameters, n, length(spec.diameters))
    diameters = Vector{T}(undef, n)
    @inbounds for i in 1:n
        diameters[i] = T(spec.diameters[i])
    end
    return diameters
end
