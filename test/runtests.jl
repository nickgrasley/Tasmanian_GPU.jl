using Tasmanian_GPU
using Test

# write your own tests here
@testset "Basics" begin

    dim = 2
    out = 1
    depth = 5
    tsg = Tasmanian.TasmanianSG(dim, out, depth)
    @test Tasmanian.getVersion() == VersionNumber("8.1")
    Tasmanian.makeLocalPolynomialGrid!(tsg)
    @test Tasmanian.isLocalPolynomial(tsg)
    @test !Tasmanian.isGlobal(tsg)
    @test Tasmanian.getNumDimensions(tsg) == dim
    @test Tasmanian.getNumOutputs(tsg) == out
    @test size(Tasmanian.getPoints(tsg)) == (dim, 145)
end

@testset "run basic" begin
    t = Tasmanian.run()
    @test isa(t,TasmanianSG)
end

@testset "run examples" begin
    Tasmanian.ex1()
    Tasmanian.ex2()
    Tasmanian.ex3()
end

@testset "BasciIO" begin
    include("testBasicIO.jl")
    include("testCommon.jl")
    checkCopySubgrid()
end
