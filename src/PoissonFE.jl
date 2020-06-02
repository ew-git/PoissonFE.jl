module PoissonFE

import StatsFuns: logsumexp
using Optim

include("fit.jl")

export PoissonFEModel, fit

end
