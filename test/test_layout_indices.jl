using Agate
using Test

using Agate.Configuration:
    Plankton, DiameterRangeSpecification, realize_model_layout, component_diameters

@testset "Size-structure input normalization" begin
    components = (
        Z=Plankton(; states=:carbon, reference_state=:carbon),
        P=Plankton(; states=:carbon, reference_state=:carbon),
    )
    layout = realize_model_layout(
        components;
        plankton_pfts=(Z=(:Z,), P=(:P,)),
        pft_size_structures=(
            Z=DiameterRangeSpecification(2, 20.0, 100.0, :linear_splitting),
            P=(n=3, min_esd=2.0, max_esd=10.0, splitting=:log_splitting),
        ),
    )

    @test component_diameters(layout, :Z) == (20.0, 100.0)
    @test collect(component_diameters(layout, :P)) ≈ [2.0, sqrt(20.0), 10.0]

    bad_path = "plankton PFT :P size_structure"
    for invalid in (
        Float64[], [1.0, 0.0], [1.0, Inf], [true],
        (n=0, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=true, min_esd=1.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=0.0, max_esd=2.0, splitting=:log_splitting),
        (n=2, min_esd=2.0, max_esd=1.0, splitting=:log_splitting),
        (n=2, min_esd=1.0, max_esd=2.0, splitting=:unsupported),
    )
        message = argument_error_message(() -> realize_model_layout(
            (P=Plankton(; states=:carbon, reference_state=:carbon),);
            plankton_pfts=(P=(:P,),), pft_size_structures=(P=invalid,),
        ))
        @test occursin(bad_path, message)
    end

end
