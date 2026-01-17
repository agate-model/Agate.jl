# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Directory / pretty printing
# ----------------------------------------------------------------------------

"""Return a lightweight directory of parameters for `factory`."""
function parameter_directory(factory)
    reg = parameter_registry(factory)
    return map(reg.specs) do s
        default_form = if isnothing(s.provider)
            (s.shape === :matrix && (s.name in (:palatability_matrix, :assimilation_matrix))) ? :derived : :required
        else
            typeof(s.provider)
        end

        (name=s.name,
         shape=s.shape,
         value_kind=s.value_kind,
         doc=s.doc,
         default=default_form)
    end
end

@inline function _provider_string(name::Symbol, shape::Symbol, p)
    if isnothing(p)
        return (shape === :matrix && (name in (:palatability_matrix, :assimilation_matrix))) ? "DERIVED" : "REQUIRED"
    end
    return string(typeof(p))
end

function show(io::IO, ::MIME"text/plain", reg::ParamRegistry)
    println(io, "Agate.ParamRegistry with ", length(reg.specs), " parameters")
    for s in reg.specs
        println(io, "  * ", s.name, " (", s.shape, ", ", s.value_kind, ")")
        println(io, "      default: ", _provider_string(s.name, s.shape, s.provider))
        if !isempty(s.doc)
            for line in split(s.doc, '\n')
                println(io, "      ", line)
            end
        end
    end
end
