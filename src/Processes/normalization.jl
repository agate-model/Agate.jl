"""Author-facing scientific model definition."""
struct ModelDefinition{C,P,A}
    components::C
    processes::P
    parameters::A
end

function ModelDefinition(; components::NamedTuple, processes::NamedTuple, parameters=nothing)
    return ModelDefinition(components, processes, parameters)
end

function ModelDefinition(family::AbstractModelFamily)
    return ModelDefinition(;
        components=default_components(family),
        processes=default_processes(family),
        parameters=parameter_definitions(family),
    )
end

"""Named, validated process instance in canonical normalized model state.

`facts` contains setup-only decisions that compilation may trust, such as resolved
population-state references and Growth routing ownership.
"""
struct NamedProcess{P<:AbstractProcess,F}
    id::Symbol
    process::P
    facts::F
end

NamedProcess(id::Symbol, process::P) where {P<:AbstractProcess} =
    NamedProcess(id, process, NamedTuple())

process_id(process::NamedProcess) = process.id
formulation(process::NamedProcess) = formulation(process.process)
factors(process::NamedProcess) = factors(process.process)
participants(process::NamedProcess) = participants(process.process)

"""Zero-or-one setup-time qualifier attached to a local parameter slot."""
struct Qualifier
    axis::Symbol
    value::Symbol
end

"""Formulation-local declaration of one semantic parameter slot.

Dimensionality is structural: zero, one, or two declared axes imply scalar, vector, or
matrix values respectively. `qualify` identifies repeated semantic instances without
changing storage dimensionality. For a scalar slot, a qualifier that is also a process
participant role can still provide ecological applicability without becoming a storage axis.
"""
struct ParameterSlot{A<:Tuple}
    name::Symbol
    axes::A
    qualify::Union{Nothing,Symbol}
end

function ParameterSlot(
    name::Symbol,
    axes::Tuple=();
    qualify::Union{Nothing,Symbol}=nothing,
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter slot axes currently support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    return ParameterSlot(name, axes, qualify)
end

parameter_slots(::AbstractFormulation) = ()
parameter_slots(::Smith) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
)
parameter_slots(::Geider) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
    ParameterSlot(:chlorophyll_to_carbon_ratio, (:population,)),
)
parameter_slots(::Monod) = (ParameterSlot(:K, (:population,); qualify=:resource),)
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness),)
parameter_slots(::Q10) = (
    ParameterSlot(:q10),
    ParameterSlot(:reference_temperature),
)
parameter_slots(::IdealizedGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::PreferentialGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::HeterotrophicConsumption) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:resource,)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::LinearMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::QuadraticMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::LinearRemineralization) = (
    ParameterSlot(:rate; qualify=:source),
)
parameter_slots(products::Products) = length(products.targets) == 1 ? () :
    (ParameterSlot(:fraction; qualify=:product),)
parameter_slots(::FixedStoichiometry) = (ParameterSlot(:ratio; qualify=:currency),)

"""Resolved setup-time mapping from one local parameter slot to model storage.

The scientific-tree location and qualifier identify the local slot during compilation.
`storage_axes=nothing` means slot-local storage; explicit axes identify the global runtime
storage coordinate system used by the bound model parameter.
"""
struct ParameterBinding{P,Q,A,S}
    process::Symbol
    path::P
    slot::Symbol
    qualifier::Q
    axes::A
    parameter::Symbol
    storage_axes::S
end

"""Concrete participant applicability of one resolved parameter binding."""
struct ParameterApplicability{B,C,T}
    binding::B
    axis_components::C
    axis_classes::T
end

_parameter_slot_source(node::Union{AbstractFormulation,AbstractStoichiometry,Products}) = node
_parameter_slot_source(node) = formulation(node)

function _slot_qualifier(slot::ParameterSlot, context::NamedTuple)
    isnothing(slot.qualify) && return nothing
    hasproperty(context, slot.qualify) || return nothing
    value = getproperty(context, slot.qualify)
    value isa Symbol || throw(
        ArgumentError("parameter slot qualifier :$(slot.qualify) must identify one Symbol"),
    )
    return Qualifier(slot.qualify, value)
end

function _slot_metadata(
    named::NamedProcess,
    path::Tuple,
    slot::ParameterSlot,
    context::NamedTuple,
)
    all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter binding path must contain only Symbols"),
    )
    return (;
        process=process_id(named),
        path,
        slot=slot.name,
        qualifier=_slot_qualifier(slot, context),
        axes=slot.axes,
    )
end

_qualifier_key(::Nothing) = nothing
_qualifier_key(qualifier::Qualifier) = (qualifier.axis, qualifier.value)
_binding_key(process::Symbol, path::Tuple, slot::Symbol, qualifier) =
    (process, path, slot, _qualifier_key(qualifier))
_binding_key(binding::ParameterBinding) =
    _binding_key(binding.process, binding.path, binding.slot, binding.qualifier)
_binding_key(metadata::NamedTuple) =
    _binding_key(metadata.process, metadata.path, metadata.slot, metadata.qualifier)

function _emit_parameter_slots!(
    uses::Vector{Any},
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
)
    bindings = _validate_binding_slot_names(node)
    for slot in parameter_slots(_parameter_slot_source(node))
        qualifier = _slot_qualifier(slot, context)
        metadata = _slot_metadata(named, path, slot, context)
        parameter, explicit = _binding_value(bindings, slot, qualifier)
        push!(uses, (; metadata..., parameter, explicit))
    end
    return nothing
end

function _visit_factor_slots!(
    uses::Vector{Any}, named::NamedProcess, path::Tuple, factor::AbstractFactor
)
    _emit_parameter_slots!(
        uses, named, path, factor; context=factor_parameter_context(factor)
    )
    for (name, child) in pairs(_validated_factor_children(factor))
        _visit_factor_slots!(
            uses, named, factor_child_path(path, factor, name), child
        )
    end
    return nothing
end

function _visit_product_slots!(
    uses::Vector{Any}, named::NamedProcess, path::Tuple, products::Products
)
    for product in keys(products.fractions)
        _emit_parameter_slots!(uses, named, path, products; context=(product=product,))
    end
    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return nothing
    currencies = keys(first(values(products.targets)))
    for currency in currencies
        currency === stoichiometry.reference && continue
        _emit_parameter_slots!(
            uses,
            named,
            (path..., :stoichiometry),
            stoichiometry;
            context=(currency=currency,),
        )
    end
    return nothing
end

function _visit_process_slots!(uses::Vector{Any}, named::NamedProcess)
    process = named.process
    if process isa Mortality
        for population in process.populations
            _emit_parameter_slots!(uses, named, (), process; context=(population=population,))
        end
    elseif process isa Remineralization
        for source in process.sources
            _emit_parameter_slots!(uses, named, (), process; context=(source=source,))
        end
    else
        _emit_parameter_slots!(uses, named, (), process)
    end

    for (name, factor) in pairs(factors(process))
        _visit_factor_slots!(uses, named, (:factors, name), factor)
    end

    products = process_products(process)
    isnothing(products) || _visit_product_slots!(uses, named, product_path(process), products)

    if process isa Growth && !isnothing(process.stoichiometry)
        nutrients = getproperty(process.factors, named.facts.routing.factor)
        for currency in keys(nutrients.responses)
            _emit_parameter_slots!(
                uses,
                named,
                (:stoichiometry,),
                process.stoichiometry;
                context=(currency=currency,),
            )
        end
    end
    return nothing
end

"""Setup-time normalized scientific model definition.

`parameter_bindings` is the canonical ordered contract; `parameter_lookup` is a transient
setup cache used while lowering processes.
"""
struct NormalizedModelDefinition{C,P,A,O,D,B,L}
    components::C
    processes::P
    parameters::A
    derived_parameter_order::O
    driver_identities::D
    parameter_bindings::B
    parameter_lookup::L
end

"""Return the canonical external-driver identities required by a normalized model."""
driver_identities(definition::NormalizedModelDefinition) = definition.driver_identities

"""Return resolved local-slot-to-model-parameter bindings."""
parameter_bindings(definition::NormalizedModelDefinition) = definition.parameter_bindings

function _parameter_binding(
    definition::NormalizedModelDefinition,
    process::Symbol,
    path::Tuple,
    slot::Symbol,
    qualifier,
)
    key = _binding_key(process, path, slot, qualifier)
    return get(definition.parameter_lookup, key) do
        throw(ArgumentError(
            "no model parameter is bound to slot :$slot at process :$process path $path qualifier $qualifier",
        ))
    end
end

"""Resolve all parameter slots for one scientific node from its formulation schema."""
function parameter_slot_bindings(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
)
    slot_source = _parameter_slot_source(node)
    slots = parameter_slots(slot_source)
    names = Tuple(slot.name for slot in slots)
    bindings = Tuple(
        begin
            qualifier = _slot_qualifier(slot, context)
            _parameter_binding(
                definition, process_id(named), path, slot.name, qualifier
            )
        end for slot in slots
    )
    return NamedTuple{names}(bindings)
end

function _validated_factor_children(factor::AbstractFactor)
    children = factor_children(factor)
    children isa NamedTuple || throw(
        ArgumentError("factor children for $(typeof(factor)) must be a NamedTuple"),
    )
    all(child -> child isa AbstractFactor, values(children)) || throw(
        ArgumentError("factor children for $(typeof(factor)) must be process factors"),
    )
    return children
end

function _factor_component_references(factor::AbstractFactor)
    references = Symbol[
        input.component for input in factor_inputs(factor) if input isa FactorComponent
    ]
    for child in values(_validated_factor_children(factor))
        append!(references, _factor_component_references(child))
    end
    return Tuple(references)
end

function _product_component_references(products::Products)
    return Tuple(Iterators.flatten(
        target isa Symbol ? (target,) : Tuple(values(target)) for target in values(products.targets)
    ))
end

function _process_component_references(process::AbstractProcess)
    references = Symbol[]
    for values_for_role in values(participants(process))
        append!(references, values_for_role)
    end
    for factor in values(factors(process))
        append!(references, _factor_component_references(factor))
    end
    products = process_products(process)
    isnothing(products) || append!(references, _product_component_references(products))
    return Tuple(references)
end

function _scalar_pool(components, name::Symbol, id::Symbol, label::AbstractString)
    pool = getproperty(components, name)
    pool isa Pool || throw(ArgumentError("process :$id $label :$name must be a Pool"))
    isnothing(size_structure(pool)) || throw(
        ArgumentError("process :$id $label :$name must be a scalar Pool"),
    )
    return pool
end

function _population_state(components, name::Symbol, id::Symbol, label::AbstractString)
    population = getproperty(components, name)
    population isa Population || throw(
        ArgumentError("process :$id $label :$name must be a Population"),
    )
    state_names = keys(states(population))
    length(state_names) == 1 || throw(
        ArgumentError("process :$id $label :$name requires explicit state selection"),
    )
    return PopulationStateRef(name, only(state_names))
end

_state_currency(components, ref::PopulationStateRef) =
    state_currency(getproperty(components, ref.population), ref.state)

function _require_currency(actual, expected, id::Symbol, label::AbstractString)
    actual === expected || throw(
        ArgumentError("process :$id $label has currency :$actual, expected :$expected"),
    )
end

function _validate_state_currencies!(components, refs, expected, id, label)
    for ref in refs
        _require_currency(
            _state_currency(components, ref), expected, id,
            "$label :$(ref.population).$(ref.state)",
        )
    end
    return nothing
end

function _validate_products!(id, products, components, reference, label)
    stoichiometry = products.stoichiometry
    if isnothing(stoichiometry)
        for (name, target) in pairs(products.targets)
            target isa Symbol || throw(
                ArgumentError("process :$id $label product :$name requires a scalar target"),
            )
            pool = _scalar_pool(components, target, id, "$label product :$name target")
            _require_currency(currency(pool), reference, id, "$label product :$name target :$target")
        end
        return nothing
    end

    _require_currency(stoichiometry.reference, reference, id, "$label stoichiometric reference")
    for (name, targets) in pairs(products.targets), (target_currency, target) in pairs(targets)
        pool = _scalar_pool(components, target, id, "$label product :$name target")
        _require_currency(
            currency(pool), target_currency, id, "$label product :$name target :$target",
        )
    end
    return nothing
end

function _growth_facts(id::Symbol, process::Growth, components::NamedTuple)
    routing = Tuple(
        (name, factor) for (name, factor) in pairs(process.factors)
        if factor isa Union{NutrientResponse,Nutrients}
    )
    length(routing) == 1 || throw(ArgumentError(
        "process :$id growth must declare exactly one NutrientResponse or Nutrients routing factor",
    ))
    factor_name, nutrient_factor = only(routing)

    rate_owners = Tuple(
        name for (name, factor) in pairs(process.factors)
        if any(slot -> slot.name === :maximum_rate, parameter_slots(formulation(factor)))
    )
    length(rate_owners) == 1 || throw(ArgumentError(
        "process :$id growth must declare exactly one factor that owns the maximum_rate slot",
    ))
    population_states = Tuple(
        _population_state(components, population, id, "population")
        for population in process.populations
    )

    if nutrient_factor isa NutrientResponse
        isnothing(process.source) || throw(ArgumentError(
            "process :$id single-resource growth derives its source from the nutrient response; omit `source`",
        ))
        isnothing(process.stoichiometry) || throw(ArgumentError(
            "process :$id single-resource growth does not take fixed stoichiometry",
        ))
        resource = _scalar_pool(components, nutrient_factor.resource, id, "nutrient factor resource")
        reference = currency(resource)
        _validate_state_currencies!(components, population_states, reference, id, "population state")
        route = (mode=:single_resource, factor=factor_name, source=nutrient_factor.resource)
    else
        process.source isa Symbol || throw(
            ArgumentError("process :$id multi-resource growth requires a source component"),
        )
        process.stoichiometry isa FixedStoichiometry || throw(
            ArgumentError("process :$id multi-resource growth requires FixedStoichiometry"),
        )
        reference = process.stoichiometry.reference
        source = _scalar_pool(components, process.source, id, "growth source")
        _require_currency(currency(source), reference, id, "growth source :$(process.source)")
        _validate_state_currencies!(components, population_states, reference, id, "population state")
        for (target_currency, response) in pairs(nutrient_factor.responses)
            resource = _scalar_pool(
                components, response.resource, id, "nutrient response :$target_currency resource",
            )
            _require_currency(
                currency(resource), target_currency, id,
                "nutrient response :$target_currency resource :$(response.resource)",
            )
        end
        route = (
            mode=:multi_resource, factor=factor_name, source=process.source,
            stoichiometry=process.stoichiometry,
        )
    end
    return (;
        population_states, routing=route, maximum_rate_factor=only(rate_owners),
    )
end

function _consumption_facts(id::Symbol, process::Consumption, components::NamedTuple)
    consumer_states = Tuple(
        _population_state(components, consumer, id, "consumer") for consumer in process.consumers
    )
    reference = _state_currency(components, first(consumer_states))
    _validate_state_currencies!(components, consumer_states, reference, id, "consumer state")

    if uses_living_interactions(process.formulation)
        resources = Tuple(
            _population_state(components, resource, id, "living resource")
            for resource in process.resources
        )
        _validate_state_currencies!(components, resources, reference, id, "resource state")
    else
        resources = process.resources
        for resource in resources
            pool = getproperty(components, resource)
            pool isa Pool || throw(
                ArgumentError("process :$id heterotrophic resource :$resource must be a Pool"),
            )
            _require_currency(currency(pool), reference, id, "heterotrophic resource :$resource")
        end
    end
    isnothing(process.products) || _validate_products!(
        id, process.products, components, reference, "unassimilated products"
    )
    return (; consumer_states, resources)
end

function _mortality_facts(id::Symbol, process::Mortality, components::NamedTuple)
    population_states = Tuple(
        _population_state(components, population, id, "mortality population")
        for population in process.populations
    )
    if !isnothing(process.products)
        reference = _state_currency(components, first(population_states))
        _validate_state_currencies!(
            components, population_states, reference, id, "mortality population state"
        )
        _validate_products!(id, process.products, components, reference, "mortality products")
    end
    return (; population_states)
end

function _remineralization_facts(id::Symbol, process::Remineralization, components::NamedTuple)
    destination = _scalar_pool(components, process.destination, id, "remineralization destination")
    reference = currency(destination)
    for source in process.sources
        pool = _scalar_pool(components, source, id, "remineralization source")
        _require_currency(currency(pool), reference, id, "remineralization source :$source")
    end
    return NamedTuple()
end

_normalized_process_facts(::Symbol, ::AbstractProcess, ::NamedTuple) = NamedTuple()
_normalized_process_facts(id::Symbol, process::Growth, components::NamedTuple) =
    _growth_facts(id, process, components)
_normalized_process_facts(id::Symbol, process::Consumption, components::NamedTuple) =
    _consumption_facts(id, process, components)
_normalized_process_facts(id::Symbol, process::Mortality, components::NamedTuple) =
    _mortality_facts(id, process, components)
_normalized_process_facts(id::Symbol, process::Remineralization, components::NamedTuple) =
    _remineralization_facts(id, process, components)

function _validate_process(id::Symbol, process, components::NamedTuple)
    process isa AbstractProcess || throw(
        ArgumentError("process :$id must be an AbstractProcess; got $(typeof(process))"),
    )
    component_names = keys(components)
    missing = unique(Tuple(
        reference for reference in _process_component_references(process)
        if !(reference in component_names)
    ))
    isempty(missing) || throw(
        ArgumentError("process :$id references unknown components $missing"),
    )
    return NamedProcess(id, process, _normalized_process_facts(id, process, components))
end

function _canonical_processes(processes::NamedTuple, components::NamedTuple)
    names = sort!(collect(keys(processes)); by=String)
    names_tuple = Tuple(names)
    process_values = Tuple(
        _validate_process(name, getproperty(processes, name), components) for name in names
    )
    return NamedTuple{names_tuple}(process_values)
end

function _collect_driver_identities!(identities::Vector{Symbol}, factor::AbstractFactor)
    for input in factor_inputs(factor)
        input isa FactorDriver || continue
        input.identity in identities || push!(identities, input.identity)
    end
    for child in values(_validated_factor_children(factor))
        _collect_driver_identities!(identities, child)
    end
    return nothing
end

function _canonical_driver_identities(processes::NamedTuple)
    identities = Symbol[]
    for process in values(processes), factor in values(factors(process))
        _collect_driver_identities!(identities, factor)
    end
    sort!(identities; by=String)
    return Tuple(identities)
end

function _validate_binding_slot_names(node)
    slots = parameter_slots(_parameter_slot_source(node))
    bindings = if applicable(authored_parameter_bindings, node)
        authored_parameter_bindings(node)
    else
        isempty(slots) || throw(ArgumentError(
            "slot-owning node $(typeof(node)) must implement `authored_parameter_bindings`",
        ))
        NamedTuple()
    end
    slot_names = Tuple(slot.name for slot in slots)
    unknown = Tuple(name for name in keys(bindings) if !(name in slot_names))
    isempty(unknown) || throw(ArgumentError(
        "inline bindings contain unknown slots $unknown; expected $slot_names",
    ))
    return bindings
end

function _binding_value(bindings::NamedTuple, slot::ParameterSlot, qualifier)
    explicit = hasproperty(bindings, slot.name)
    value = explicit ? getproperty(bindings, slot.name) : slot.name
    value isa Symbol && return value, explicit

    value isa NamedTuple || throw(ArgumentError(
        "binding :$(slot.name) must be a parameter Symbol or one-level qualifier NamedTuple",
    ))
    isnothing(slot.qualify) && throw(ArgumentError(
        "binding :$(slot.name) uses a qualifier map but the slot is unqualified",
    ))
    qualifier isa Qualifier && qualifier.axis === slot.qualify || throw(ArgumentError(
        "binding :$(slot.name) requires qualifier :$(slot.qualify)",
    ))
    qualifier_value = qualifier.value
    hasproperty(value, qualifier_value) || throw(ArgumentError(
        "binding :$(slot.name) has no entry for qualifier :$qualifier_value",
    ))
    parameter = getproperty(value, qualifier_value)
    parameter isa Symbol || throw(ArgumentError(
        "binding :$(slot.name) qualifier :$qualifier_value must map to a parameter Symbol",
    ))
    return parameter, true
end

function _declared_parameter_uses(processes::NamedTuple)
    uses = Any[]
    for named in values(processes)
        _visit_process_slots!(uses, named)
    end
    keys = Tuple(_binding_key(use) for use in uses)
    length(unique(keys)) == length(keys) || throw(
        ArgumentError("normalized processes declare duplicate parameter binding keys"),
    )
    return Tuple(uses)
end

function _storage_rank(axes)
    axes === nothing && return nothing
    axes isa Symbol && return 1
    return length(axes)
end

const RESERVED_PARAMETER_KEYS = (:x, :y, :z, :t)

function _normalize_parameter_definitions(definitions)
    isnothing(definitions) && return nothing, (), Set{Symbol}()
    definitions isa NamedTuple || throw(
        ArgumentError("model parameters must be a NamedTuple of Parameter values"),
    )
    all(parameter -> parameter isa Parameter, values(definitions)) || throw(
        ArgumentError("model parameters must contain only Parameter values"),
    )

    definition_names = Tuple(keys(definitions))
    definition_set = Set(definition_names)
    for name in definition_names
        name in RESERVED_PARAMETER_KEYS && throw(
            ArgumentError("model parameters declare reserved parameter name :$name"),
        )
    end

    dependency_names = Set{Symbol}()
    derived_names = Symbol[]
    for (name, parameter) in pairs(definitions)
        default = parameter.default
        default isa DerivedDefault || continue
        push!(derived_names, name)
        for dependency in default.deps
            dependency in definition_set || throw(
                ArgumentError(
                    "parameter :$name default depends on undeclared parameter :$dependency",
                ),
            )
            push!(dependency_names, dependency)
        end
    end

    derived_set = Set(derived_names)
    resolved = Set{Symbol}()
    order = Symbol[]
    pending = copy(derived_names)
    while !isempty(pending)
        remaining = Symbol[]
        progressed = false
        for name in pending
            dependencies = Tuple(
                dep for dep in getproperty(definitions, name).default.deps if dep in derived_set
            )
            if all(dep -> dep in resolved, dependencies)
                push!(order, name)
                push!(resolved, name)
                progressed = true
            else
                push!(remaining, name)
            end
        end
        progressed || throw(
            ArgumentError(
                "derived parameter defaults contain a dependency cycle among: " *
                join(string.(Tuple(remaining)), ", "),
            ),
        )
        pending = remaining
    end
    return definitions, Tuple(order), dependency_names
end

function _normalize_parameter_bindings(
    processes::NamedTuple, definitions, dependency_names::Set{Symbol}
)
    uses = _declared_parameter_uses(processes)
    if isnothing(definitions)
        bindings = Tuple(
            ParameterBinding(
                use.process, use.path, use.slot, use.qualifier,
                use.axes, use.parameter, nothing,
            )
            for use in uses
        )
        return bindings
    end

    definition_names = Set(keys(definitions))
    for use in uses
        use.parameter in definition_names || throw(ArgumentError(
            "inline binding for process :$(use.process) path $(use.path) slot :$(use.slot) " *
            "qualifier $(use.qualifier) names undeclared parameter :$(use.parameter)",
        ))
    end

    by_parameter = Dict{Symbol,Vector{Any}}()
    for use in uses
        push!(get!(by_parameter, use.parameter, Any[]), use)
    end
    for (parameter, parameter_uses) in by_parameter
        if length(parameter_uses) > 1 && any(use -> !use.explicit, parameter_uses)
            locations = Tuple(
                (process=use.process, path=use.path, slot=use.slot, qualifier=use.qualifier)
                for use in parameter_uses
            )
            throw(ArgumentError(
                "parameter :$parameter is implicitly bound by multiple parameter slots $locations; " *
                "bind every shared use explicitly",
            ))
        end
        required_ranks = unique(Tuple(length(use.axes) for use in parameter_uses))
        length(required_ranks) == 1 || throw(ArgumentError(
            "parameter :$parameter is bound to slots with incompatible dimensionality $(Tuple(required_ranks))",
        ))
    end

    for (name, parameter) in pairs(definitions)
        parameter_uses = get(by_parameter, name, Any[])
        isempty(parameter_uses) && !(name in dependency_names) && throw(ArgumentError(
            "parameter :$name is neither bound to a process slot nor used by a derived default",
        ))
        isempty(parameter_uses) && continue
        required_rank = only(unique(Tuple(length(use.axes) for use in parameter_uses)))
        storage_rank = _storage_rank(parameter.spec.axes)
        isnothing(storage_rank) || storage_rank == required_rank || throw(ArgumentError(
            "parameter :$name storage axes imply rank $storage_rank but its bound slots require rank $required_rank",
        ))
    end

    return Tuple(
        ParameterBinding(
            use.process,
            use.path,
            use.slot,
            use.qualifier,
            use.axes,
            use.parameter,
            getproperty(definitions, use.parameter).spec.axes,
        )
        for use in uses
    )
end

"""Normalize process identity and resolve inline parameter bindings.

Process instances are canonicalized by stable process ID, so declaration order does
not change normalized scientific identity. Component ordering is preserved because it
still participates in concrete tracer realization. Local formulation slots bind directly
to stable model parameter names during normalization.
"""
function normalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Population,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Population or Pool values"))
    normalized_processes = _canonical_processes(
        definition.processes, definition.components
    )
    parameters, derived_order, dependency_names = _normalize_parameter_definitions(
        definition.parameters
    )
    bindings = _normalize_parameter_bindings(
        normalized_processes, parameters, dependency_names
    )
    lookup = Dict(_binding_key(binding) => binding for binding in bindings)
    return NormalizedModelDefinition(
        definition.components,
        normalized_processes,
        parameters,
        derived_order,
        _canonical_driver_identities(normalized_processes),
        bindings,
        lookup,
    )
end

function _axis_components(
    process::NamedProcess, binding::ParameterBinding, axis::Symbol
)
    qualifier = binding.qualifier
    if qualifier isa Qualifier && qualifier.axis === axis
        return (qualifier.value,)
    end
    process_participants = participants(process)
    hasproperty(process_participants, axis) || throw(
        ArgumentError(
            "parameter applicability axis :$axis is not a participant role of process :$(process_id(process))",
        ),
    )
    return getproperty(process_participants, axis)
end

function _axis_classes(layout::ModelLayout, components::Tuple)
    classes = Symbol[]
    for component in components
        hasproperty(layout.component_classes, component) || throw(
            ArgumentError("parameter applicability references unrealized component :$component"),
        )
        append!(classes, component_classes(layout, component))
    end
    return Tuple(classes)
end

"""Resolve each semantic parameter axis onto ecological component class identities.

The result is setup-time applicability metadata. Population state multiplicity does not
change parameter-axis length: applicability follows ecological classes rather than physical
state tracers. Scalar qualified slots retain participant applicability when their qualifier
names a process participant role; non-participant qualifiers such as currency do not create
ecological axes.
"""
function _applicability_axes(process::NamedProcess, binding::ParameterBinding)
    isempty(binding.axes) || return binding.axes
    qualifier = binding.qualifier
    qualifier isa Qualifier || return ()
    process_participants = participants(process)
    return hasproperty(process_participants, qualifier.axis) ? (qualifier.axis,) : ()
end

function resolve_parameter_applicability(
    definition::NormalizedModelDefinition, layout::ModelLayout
)
    return map(definition.parameter_bindings) do binding
        process = getproperty(definition.processes, binding.process)
        applicability_axes = _applicability_axes(process, binding)
        axis_components = map(
            axis -> _axis_components(process, binding, axis), applicability_axes
        )
        axis_classes = map(components -> _axis_classes(layout, components), axis_components)
        ParameterApplicability(binding, axis_components, axis_classes)
    end
end
