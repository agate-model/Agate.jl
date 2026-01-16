# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Directory / pretty printing
# ----------------------------------------------------------------------------

"""Return a lightweight directory of parameters for `factory`."""
function parameter_directory(factory)
    reg = parameter_registry(factory)
    return map(reg.specs) do s
        default_form = isnothing(s.provider) ? :required : typeof(s.provider)
        (name=s.name,
         shape=s.shape,
         missing_policy=s.missing_policy,
         value_kind=s.value_kind,
         doc=s.doc,
         default=default_form)
    end
end

@inline _provider_string(p) = isnothing(p) ? "REQUIRED" : string(typeof(p))

function show(io::IO, ::MIME"text/plain", reg::ParamRegistry)
    println(io, "Agate.ParamRegistry with $(length(reg.specs)) parameters")
    for s in reg.specs
        println(io, "  * ", s.name, " (", s.shape, ", ", s.value_kind, ") [", s.missing_policy, "]")
        println(io, "      default: ", _provider_string(s.provider))
        if !isempty(s.doc)
            for line in split(s.doc, '\n')
                println(io, "      ", line)
            end
        end
    end
end
