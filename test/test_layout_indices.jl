using Agate
using Test

using Agate.Configuration:
    Plankton, realize_model_layout, component_entities, component_diameters

@testset "Size-structure input normalization" begin
    components = (
        Z=Plankton(; states=:carbon, reference_state=:carbon),
        P=Plankton(; states=:carbon, reference_state=:carbon),
    )
    layout = realize_model_layout(
        components;
        plankton_pfts=(Z=(:Z,), P=(:P,)),
        pft_size_structures=(
            P=(n=3, min_esd=2.0, max_esd=10.0, spacing=:log),
            Z=(n=2, min_esd=20.0, max_esd=100.0, spacing=:linear),
        ),
    )

    @test component_entities(layout, :Z) == (:Z_1, :Z_2)
    @test component_entities(layout, :P) == (:P_1, :P_2, :P_3)
    @test component_diameters(layout, :Z) == (20.0, 100.0)
    @test collect(component_diameters(layout, :P)) ≈ [2.0, sqrt(20.0), 10.0]


    unsized = realize_model_layout(
        (P=Plankton(; states=:carbon, reference_state=:carbon),);
        plankton_pfts=(P=(:small, :large),),
        pft_size_structures=(small=nothing, large=nothing),
    )
    @test component_entities(unsized, :P) == (:small, :large)
    @test component_diameters(unsized, :P) === nothing

    one_size_class = realize_model_layout(
        (P=Plankton(; states=:carbon, reference_state=:carbon),);
        plankton_pfts=(P=(:P,),),
        pft_size_structures=(P=(n=1, min_esd=5.0, max_esd=5.0, spacing=:log),),
    )
    @test component_entities(one_size_class, :P) == (:P_1,)
    @test component_diameters(one_size_class, :P) == (5.0,)

    mixed = realize_model_layout(
        (P=Plankton(; states=:carbon, reference_state=:carbon),);
        plankton_pfts=(P=(:plain, :sized),),
        pft_size_structures=(plain=nothing, sized=[5.0]),
    )
    @test component_entities(mixed, :P) == (:plain, :sized_1)
    @test component_diameters(mixed, :P) == (nothing, 5.0)

    bad_path = "plankton PFT :P size_structure"
    for invalid in (
        Float64[], [1.0, 0.0], [1.0, Inf], [true],
        (n=0, min_esd=1.0, max_esd=2.0, spacing=:log),
        (n=true, min_esd=1.0, max_esd=2.0, spacing=:log),
        (n=2, min_esd=0.0, max_esd=2.0, spacing=:log),
        (n=2, min_esd=2.0, max_esd=1.0, spacing=:log),
        (n=2, min_esd=1.0, max_esd=2.0, spacing=:unsupported),
    )
        message = argument_error_message(() -> realize_model_layout(
            (P=Plankton(; states=:carbon, reference_state=:carbon),);
            plankton_pfts=(P=(:P,),), pft_size_structures=(P=invalid,),
        ))
        @test occursin(bad_path, message)
    end

end
