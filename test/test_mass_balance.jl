using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField
using Agate.Constructors

const year = years = 365day

@testset "mass_balance" begin
    @testset "size_structured_NPZD box model" begin
        model = construct_size_structured_NPZD()

        bgc_model = Biogeochemistry(
            model(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        box_model = BoxModel(; biogeochemistry=bgc_model)
        set!(box_model; N=7, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

        function estimate_mass(box_model)
            return box_model.fields.N.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] +
                   box_model.fields.P2.data[1, 1, 1] +
                   box_model.fields.Z1.data[1, 1, 1] +
                   box_model.fields.Z2.data[1, 1, 1] +
                   box_model.fields.D.data[1, 1, 1]
        end

        initial_nitrogen_mass = estimate_mass(box_model)

        for n in 1:1000
            time_step!(box_model, 0.1)
        end

        final_nitrogen_mass = estimate_mass(box_model)
        @test isapprox(initial_nitrogen_mass, final_nitrogen_mass)
    end
    @testset "thunder_egg_1 box model" begin
        model = construct_thunder_egg_1()

        bgc_model = Biogeochemistry(
            model(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        box_model = BoxModel(; biogeochemistry=bgc_model)
        set!(box_model; DIN=7, PO4=3, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, DOC=0, POC=0.0)

        function estimate_carbon_mass(box_model)
            return box_model.fields.DIC.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] +
                   box_model.fields.P2.data[1, 1, 1] +
                   box_model.fields.Z1.data[1, 1, 1] +
                   box_model.fields.Z2.data[1, 1, 1] +
                   box_model.fields.POC.data[1, 1, 1] +
                   box_model.fields.DOC.data[1, 1, 1]
        end

        function estimate_nitrogen_mass(box_model)
            return box_model.fields.DIN.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.P2.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.Z1.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.Z2.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.POC.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.DOC.data[1, 1, 1] * model().nitrogen_to_carbon
        end

        function estimate_phosphorus_mass(box_model)
            return box_model.fields.PO4.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.P2.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.Z1.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.Z2.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.POC.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.DOC.data[1, 1, 1] * model().phosphorus_to_carbon
        end

        carbon_mass_records = Float64[]
        nitrogen_mass_records = Float64[]
        phosphorus_mass_records = Float64[]
        initial_carbon_mass = estimate_carbon_mass(box_model)
        initial_nitrogen_mass = estimate_nitrogen_mass(box_model)
        initial_phosphorus_mass = estimate_phosphorus_mass(box_model)
        push!(carbon_mass_records, initial_carbon_mass)
        push!(nitrogen_mass_records, initial_nitrogen_mass)
        push!(phosphorus_mass_records, initial_phosphorus_mass)

        for n in 1:1000
            time_step!(box_model, 1.0)
            if n % 100 == 0
                push!(carbon_mass_records, estimate_carbon_mass(box_model))
                push!(nitrogen_mass_records, estimate_nitrogen_mass(box_model))
                push!(phosphorus_mass_records, estimate_phosphorus_mass(box_model))
            end
        end

        final_carbon_mass = estimate_carbon_mass(box_model)
        final_nitrogen_mass = estimate_nitrogen_mass(box_model)
        final_phosphorus_mass = estimate_phosphorus_mass(box_model)
        push!(carbon_mass_records, final_carbon_mass)
        push!(nitrogen_mass_records, final_nitrogen_mass)
        push!(phosphorus_mass_records, final_phosphorus_mass)

        atol = 1e-8

        # Ensure there is no systematic mass decrease beyond numerical noise
        @test all(diff(carbon_mass_records) .≥ -atol)
        @test all(diff(nitrogen_mass_records) .≥ -atol)
        @test all(diff(phosphorus_mass_records) .≥ -atol)
    end
end
