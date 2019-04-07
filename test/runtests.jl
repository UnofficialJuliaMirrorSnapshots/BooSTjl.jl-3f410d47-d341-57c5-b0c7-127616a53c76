using BooSTjl
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
function test()
    x::Matrix{Float64} = [randn(10) randn(10)]
    y::Vector{Float64} = randn(10)

    t1 = BooST(x,y, M = 5,display = false)
    t2 = predictBooST(t1,x)
    t3 = BooSTMore(x,y,t1,M = 5, display = false)
    t4 = estimate_derivatives(t1,x,1)
    return t1, t2, t3, t4

end

(t1,t2,t3,t4) = test()

@test typeof(t1) == Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}}
@test typeof(t2) == Vector{Float64}
@test typeof(t3) == Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}}
@test typeof(t4) == Vector{Float64}
