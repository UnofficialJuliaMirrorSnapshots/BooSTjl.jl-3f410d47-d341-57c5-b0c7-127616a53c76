```@docs
BooST(x::Matrix{Float64},y::Vector{Float64}; v::Float64 = 0.2, p::Float64 = 2/3, d_max::Int64 = 4, gamma::Vector{Float64} = collect(0.5:0.01:5), M::Int64 = 300, display::Bool = true, node_obs::Int64 = 1)
estimate_derivatives(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},x::Matrix{Float64},variable::Int64)
predictBooST(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},newx::Matrix{Float64})
BooSTMore(x::Matrix{Float64},y::Vector{Float64},object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}}; M::Int64=50, display::Bool=true)
```
