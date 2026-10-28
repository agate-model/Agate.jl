import Pkg
Pkg.activate(@__DIR__)
Pkg.add(Pkg.PackageSpec(; url="https://github.com/agate-model/Onyx.jl", rev="c22c460b862d2656c9f429ccbe790ddb165d67be"))
Pkg.develop(Pkg.PackageSpec(; path=joinpath(@__DIR__, "..")))
Pkg.instantiate()
