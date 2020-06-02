# PoissonFE

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ew-git.github.io/PoissonFE.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ew-git.github.io/PoissonFE.jl/dev)
[![Build Status](https://travis-ci.com/ew-git/PoissonFE.jl.svg?branch=master)](https://travis-ci.com/ew-git/PoissonFE.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ew-git/PoissonFE.jl?svg=true)](https://ci.appveyor.com/project/ew-git/PoissonFE-jl)
[![Coverage](https://codecov.io/gh/ew-git/PoissonFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ew-git/PoissonFE.jl)
[![Coverage](https://coveralls.io/repos/github/ew-git/PoissonFE.jl/badge.svg?branch=master)](https://coveralls.io/github/ew-git/PoissonFE.jl?branch=master)

A Julia package to estimate Poisson regression with fixed effects. 
The fixed effects are not estimated, similar to Stata's command `xtpoisson y x, fe`. 
Consequently, estimation should be faster than [GLM.jl](https://github.com/JuliaStats/GLM.jl) for problems with many (thousands) of fixed effects. 
The package only supports one variable defining the fixed effects. 

To be implemented:
  * Standard errors from [Wooldridge (1999)](https://doi.org/10.1016/S0304-4076%2898%2900033-5)

## Installation

```julia
]add https://github.com/ew-git/PoissonFE.jl
```

## Example

The outcomes `y` and the predictors `X` should have the same Float type. 
The array `X` should be an _n_ x _k_ matrix where _n_ is the number of observations 
and _k_ is the number of parameters, excluding the fixed effects. 

```julia
y = [1.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.0, 1.0, 0.0, 0.0, 6.0,
2.0, 3.0, 0.0, 1.0, 6.0, 0.0, 7.0, 0.0, 21.0]
id = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
X = [0.398106 0.496961; -0.612026 -0.224875; 0.34112 -1.11714;
-1.12936 -0.394995; 1.43302 1.54983; 1.9804 -0.743514; -0.367221 -2.33171;
-1.04413 0.812245; 0.56972 -0.501311; -0.135055 -0.510887;
2.40162 -1.21536; -0.03924 -0.0225586; 0.689739 0.701239;
0.0280022 -0.587482; -0.743273 -0.606728; 0.188792 1.09664;
-1.80496 -0.24751; 1.46555 -0.159902; 0.153253 -0.625778;
2.17261 0.900435]
m = PoissonFEModel(y, X, id)
result = fit(m) # The coefficients of x1 and x2, [0.9244995462032737, 0.8693710731458504]
```
