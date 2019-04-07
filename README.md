
[![Build Status](https://travis-ci.com/gabrielrvsc/BooSTjl.jl.svg?branch=master)](https://travis-ci.com/gabrielrvsc/BooSTjl.jl)

# BooSTjl

Set of functions to estimate the BooST model 

## Installation

```
Pkg.clone("https://github.com/gabrielrvsc/BooSTjl.jl")
```

## List of Functions

- `BooST` estimates the BooST
- `BooSTMore` estimates more trees for a previously estimated BooST
- `predictBooST` Estimates predicted values for a given matrix of characteristics
- `estimate_derivatives` estimates the derivative of y with respect to `x[:,variable]`