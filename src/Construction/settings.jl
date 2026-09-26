using ..ModelFamilies: AbstractModelFamily, ModelSetting, setting_definitions
using ..Processes: parameter_domain_valid

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
        parameter_domain_valid(value, definition.domain) || throw(
            ArgumentError(
                "model setting :$name must satisfy domain :$(definition.domain); got $(repr(value))"
            ),
        )
        value
    end
    return NamedTuple{names}(values)
end
