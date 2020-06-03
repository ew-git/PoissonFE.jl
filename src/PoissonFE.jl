module PoissonFE

import StatsFuns: logsumexp
using Optim
using LinearAlgebra

include("fit.jl")

export PoissonFEModel, fit

end
