<a id='BooSTjl.BooST-Tuple{Array{Float64,2},Array{Float64,1}}' href='#BooSTjl.BooST-Tuple{Array{Float64,2},Array{Float64,1}}'>#</a>
**`BooSTjl.BooST`** &mdash; *Method*.



```
BooST(x::Matrix{Float64},y::Vector{Float64}; v::Float64 = 0.2, p::Float64 = 2/3, d_max::Int64 = 4, gamma::Vector{Float64} = collect(0.5:0.01:5), M::Int64 = 300, display::Bool = true, node_obs::Int64 = 1)
```

**Arguments**

  * **`x`** Matrix of characteristics
  * **`y`** Response variable vector
  * **`v`** Step size (shrinkage)
  * **`p`** Proportion of variables tested in each new split (default 2/3)
  * **`d_max`** number of splits
  * **`gamma`** Transiction function intensity. Bigger numbers makes the transition less smoth. The default is a sequence of values (0.5:5) to be randomized in each new node. Multiple values may be supplied in a vector to increase the model randomness
  * **`M`** Number of Boosting iterations to add
  * **`display`** if `true`, prints iteration counter
  * **`node_obs`** Equivalent to the minimum number of observations in a termina node for a discrete tree

Estimates the BooST and returns a tupple with the following objects below. The entire tupple is required to call the functions `BooSTMore`, `predictBooST` and `estimate_derivatives`.

**Outputs (Tupple order)**

save_tree, phi, brmse, ybar, params, nvar, save_rho

  * **`model`** Dict with information of each tree
  * **`phi`** Fitted values
  * **`brmse`** RMSE in each iteration
  * **`ybar`** Mean of y
  * **`params`** Parameters from call
  * **`nvar`** Number of characteristics
  * **`rho`** Vector with the parameter phi estimated in each iteration


<a target='_blank' href='https://github.com/gabrielrvsc/BooSTjl/blob/3ffe5b6cd888b5e37dbf01debcbd39074b49b71f/src/export.jl#L108-L137' class='documenter-source'>source</a><br>

<a id='BooSTjl.estimate_derivatives-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2},Int64}' href='#BooSTjl.estimate_derivatives-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2},Int64}'>#</a>
**`BooSTjl.estimate_derivatives`** &mdash; *Method*.



```
estimate_derivatives(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},x::Matrix{Float64},variable::Int64)
```

**Arguments**

  * **`object`** Output object from `BooST` function
  * **`x`** Matrix of characteristics
  * **`variable`** Column index of the desired variable

Returns the derivative of y with respect to `x[:,variable]`.


<a target='_blank' href='https://github.com/gabrielrvsc/BooSTjl/blob/3ffe5b6cd888b5e37dbf01debcbd39074b49b71f/src/export.jl#L3-L13' class='documenter-source'>source</a><br>

<a id='BooSTjl.predictBooST-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2}}' href='#BooSTjl.predictBooST-Tuple{Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}},Array{Float64,2}}'>#</a>
**`BooSTjl.predictBooST`** &mdash; *Method*.



```
predictBooST(object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}},newx::Matrix{Float64})
```

**Arguments**

  * **`object`** Output object from `BooST` function
  * **`newx`** Matrix of characteristics

Returns predicted values from a BooST model. If `newx = x` returns the fitted values.


<a target='_blank' href='https://github.com/gabrielrvsc/BooSTjl/blob/3ffe5b6cd888b5e37dbf01debcbd39074b49b71f/src/export.jl#L80-L89' class='documenter-source'>source</a><br>

<a id='BooSTjl.BooSTMore-Tuple{Array{Float64,2},Array{Float64,1},Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}}}' href='#BooSTjl.BooSTMore-Tuple{Array{Float64,2},Array{Float64,1},Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64,Array{Float64,1}}}'>#</a>
**`BooSTjl.BooSTMore`** &mdash; *Method*.



```
BooSTMore(x::Matrix{Float64},y::Vector{Float64},object::Tuple{Dict{Float64,Array{Float64,2}},Array{Float64,1},Array{Float64,1},Float64,Tuple{Float64,Float64,Int64,Array{Float64,1},Int64,Float64},Int64, Vector{Float64}}; M::Int64=50, display::Bool=true)
```

**Arguments**

  * **`x`** Matrix of characteristics
  * **`y`** Response variable vector
  * **`object`** Output object from `BooST` function
  * **`M`** Number of Boosting iterations to add
  * **`display`** if `true`, prints iteration counter

Estimates more trees for a previously estimated BooST and returns the updated BooST.


<a target='_blank' href='https://github.com/gabrielrvsc/BooSTjl/blob/3ffe5b6cd888b5e37dbf01debcbd39074b49b71f/src/export.jl#L31-L43' class='documenter-source'>source</a><br>

