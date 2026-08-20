using Agate.Compilation:
    model_fluxes,
    group_fluxes,
    compile_model_tendencies
using Agate.Configuration: realize_components
using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Agate.ModelFamilies: default_components
using Agate.Processes: ModelDefinition, driver_identities, normalize_model

const NIPIZD_PROCESS_TRACER_ORDER = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
const NIPIZD_PROCESS_ARGS = (
    0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, 0.01, 0.02, 100.0
)

function nipizd_process_compilation(::Type{T}=Float64) where {T<:Real}
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    normalized = normalize_model(ModelDefinition(family))
    layout = realize_components(default_components(family); scalar_type=T)
    context = Agate.Configuration.parse_community(
        T, default_nipizd_community(); biogeochem_tracers=(:N, :D)
    )
    fluxes = model_fluxes(normalized, layout, context)
    grouped = group_fluxes(
        fluxes; target_order=NIPIZD_PROCESS_TRACER_ORDER
    )
    compiled = compile_model_tendencies(
        normalized, layout, context; target_order=NIPIZD_PROCESS_TRACER_ORDER
    )
    drivers = driver_identities(normalized)
    return (; normalized, layout, context, fluxes, grouped, compiled, drivers)
end

function full_process_bgc()
    return Agate.Models.NiPiZD.construct(;
        parameters=(
            detritus_remineralization=1.4e-6,
            mortality_export_fraction=0.23,
            linear_mortality=(Z_1=1.1e-6, Z_2=1.3e-6, P_1=0.7e-6, P_2=0.9e-6),
            quadratic_mortality=(Z_1=1.7e-6, Z_2=2.1e-6),
            maximum_growth_rate=(P_1=2.2e-5, P_2=1.7e-5),
            nutrient_half_saturation=(P_1=0.21, P_2=0.34),
            alpha=(P_1=2.4e-6, P_2=1.8e-6),
            maximum_predation_rate=(Z_1=1.1e-4, Z_2=0.8e-4),
            holling_half_saturation=(Z_1=0.12, Z_2=0.18),
        ),
        palatability_matrix=[0.7 0.2; 0.4 0.9],
        assimilation_matrix=[0.3 0.35; 0.4 0.45],
    )
end

@testset "Complete NiPiZD process flux compiler" begin
    compilation = nipizd_process_compilation()
    fluxes = compilation.fluxes
    grouped = compilation.grouped
    compiled = compilation.compiled

    @test compilation.drivers == (:PAR,)
    @test length(fluxes) == 36
    @test map(length, grouped) == (N=9, D=11, Z_1=4, Z_2=4, P_1=4, P_2=4)
    @test keys(compiled) == NIPIZD_PROCESS_TRACER_ORDER
    @test all(equation -> equation isa CompiledEquation, values(compiled))
    @test all(equation -> isbitstype(typeof(equation.f)), values(compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms), values(compiled)
    )

    bgc = full_process_bgc()
    tracer_index = Agate.Runtime.build_tracer_index(
        compilation.context,
        NIPIZD_PROCESS_TRACER_ORDER,
        compilation.drivers;
        n_biogeochem_tracers=2,
    )
    factory = define_tracer_functions(
        bgc.parameters,
        compiled;
        auxiliary_fields=compilation.drivers,
        tracer_index,
    )
    generated_bgc = factory(bgc.parameters)
    @test !any(type -> type === Any, fieldtypes(typeof(generated_bgc.tracer_functions)))

    generated = map(
        target -> generated_bgc(Val(target), NIPIZD_PROCESS_ARGS...),
        NIPIZD_PROCESS_TRACER_ORDER,
    )
    constructed = map(
        target -> bgc(Val(target), NIPIZD_PROCESS_ARGS...),
        NIPIZD_PROCESS_TRACER_ORDER,
    )
    @test all(process_compiler_isapprox.(generated, constructed))
    @test isapprox(sum(generated), 0; atol=10 * eps(sum(abs, generated)))
    @test process_compiler_isapprox(
        @inferred(compiled.P_1(bgc, NIPIZD_PROCESS_ARGS...)), constructed[5]
    )
end
