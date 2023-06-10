# LSP indexing solution
# https://github.com/julia-vscode/julia-vscode/issues/800#issuecomment-650085983
if isdefined(@__MODULE__, :LanguageServer)
    include("../src/IncompressibleNavierStokes.jl")
    using .IncompressibleNavierStokes
end

using Aqua
using CairoMakie
using IncompressibleNavierStokes
using LinearAlgebra
using Statistics
using Test

@testset "IncompressibleNavierStokes" begin
    @testset "Aqua" begin
        Aqua.test_all(
            IncompressibleNavierStokes;
            ambiguities = false,
            project_toml_formatting = false, # https://github.com/JuliaTesting/Aqua.jl/issues/72
        )
    end
end
