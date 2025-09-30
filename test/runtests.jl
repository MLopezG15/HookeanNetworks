using Test
using HookeanNetworks
N=16
testvert,testedg,testkvec,testkint=TriangLattice(N,"Ord",1,2,20.0)
@testset "Hookean Networks Tests" begin
    @test typeof(testvert) <: Matrix{Float64}
    @test size(testvert)== (2,(N+1)^2)
    @test typeof(testedg)<: Vector{Tuple{Int64,Int64}}
    @test typeof(testkvec) <: Matrix{Float64}
end