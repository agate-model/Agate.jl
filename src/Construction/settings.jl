using ..ModelFamilies: AbstractModelFamily, ModelSetting, setting_definitions

_setting_domain_valid(value, ::Val{:any}) = true
_setting_domain_valid(value, ::Val{:finite}) =
    value isa Real && !(value isa Bool) && isfinite(value)
_setting_domain_valid(value, ::Val{:nonnegative}) =
    _setting_domain_valid(value, Val(:finite)) && value >= zero(value)
_setting_domain_valid(value, ::Val{:positive}) =
    _setting_domain_valid(value, Val(:finite)) && value > zero(value)
_setting_domain_valid(value, ::Val{:unit_interval}) =
    _setting_domain_valid(value, Val(:finite)) && zero(value) <= value <= one(value)

function _typed_setting(value, ::Type{T}) where {T<:Real}
    value isa Real && !(value isa Bool) && return convert(T, value)
    return value
end

"""Resolve and validate family-level scientific settings for one construction."""
function resolve_model_settings(
    family::AbstractModelFamily, overrides::NamedTuple, ::Type{T}
) where {T<:Real}
    definitions = setting_definitions(family)
    unknown = Tuple(name for name in keys(overrides) if !hasproperty(definitions, name))
    isempty(unknown) || throw(
        ArgumentError("unknown model setting override(s): $(join(string.(unknown), ", "))"),
    )

    names = keys(definitions)
    values = ntuple(length(names)) do i
        name = names[i]
        definition = getproperty(definitions, name)
        definition isa ModelSetting || throw(
            ArgumentError("setting_definitions must contain only ModelSetting values"),
        )
        value = hasproperty(overrides, name) ? getproperty(overrides, name) : definition.default
        value = _typed_setting(value, T)
        _setting_domain_valid(value, Val(definition.domain)) || throw(
            ArgumentError(
                "model setting :$name must satisfy domain :$(definition.domain); got $(repr(value))"
            ),
        )
        value
    end
    return NamedTuple{names}(values)
end
