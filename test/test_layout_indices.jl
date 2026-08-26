using Agate
using Test

using Agate.Configuration:
    Population, DiameterRangeSpecification, realize_model_layout, normalize_diameters

@testset "Realized runtime input indices" begin
    layout = default_nipizd_layout(Float64; auxiliary_fields=(:PAR,))

    @test layout.tracer_indices == (N=1, D=2, Z_1=3, Z_2=4, P_1=5, P_2=6)
    @test layout.auxiliary_indices == (PAR=7,)
    @test_throws ArgumentError Agate.Compilation.input_operand(layout, :temperature)
end

@testset "Diameter input normalization" begin
    components = (Z=Population(:carbon), P=Population(:carbon))
    layout = realize_model_layout(
        components;
        population_groups=(Z=(:Z,), P=(:P,)),
        group_diameters=(
            Z=DiameterRangeSpecification(2, 20.0, 100.0, :linear_splitting),
            P=(n=3, min_esd=2.0, max_esd=10.0, splitting=:log_splitting),
        ),
    )

    @test layout.diameters[1:2] == (20.0, 100.0)
    @test collect(layout.diameters[3:5]) ≈ [2.0, sqrt(20.0), 10.0]

    bad_path = "population group :P diameters"
    for invalid in (
        Float64[], [1.0, 0.0], [1.0, Inf], [true],
        (n=0, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=true, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=0.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=2.0, max_esd=1.0, splitting=:log_splitting),
        (n=2, min_esd=1.0, max_esd=2.0, splitting=:unsupported),
    )
        err = try
            realize_model_layout(
                (P=Population(:carbon),);
                population_groups=(P=(:P,),), group_diameters=(P=invalid,),
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin(bad_path, sprint(showerror, err))
    end

    @test normalize_diameters([1.0, 2.0]).n == 2
end
