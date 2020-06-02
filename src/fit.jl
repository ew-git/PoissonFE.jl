

pois_log_like_i = function(istart, iend, η, y)
    logdenom = logsumexp(η[istart:iend]) # from StatsFuns
    third_term = zero(eltype(η))
    for t = istart:iend
        third_term += y[t] * (η[t] - logdenom)
    end
    return third_term
end

pois_log_like = function(params, x, y, ids)
    # params is an array of β
    # x is a 2d array of independent variables
    # y is an array of outcomes, should have same type as x for performance
    # ids is a vector of indices of the last id by group. In other words,
    #   if the raw ids are [1,1,1,2,2,2] then ids is [3, 6]
    likelihood = zero(eltype(x))
    η = x*params
    istart = 1
    for iend in ids
        # maybe just the relevant views of eta and y?
        likelihood += pois_log_like_i(istart, iend, η, y)
        istart = iend + 1
    end
    return likelihood
end

"""
    struct PoissonFEModel{T}

Struct to hold input arrays for fitting a Poisson fixed effects regression
"""
struct PoissonFEModel{T}
    # length of y == size(x, 1)
    y::Array{T}
    x::Array{T}
    id # array-like  for findlast
end

"""
    fit(model::PoissonFEModel{T})

Returns the estimated coefficients of a `PoissonFEModel` object.
"""
function fit(model::PoissonFEModel{T}) where T <: AbstractFloat
    id_indices = findlast.(isequal.(unique(model.id)), [model.id])
    llike_optim(params) = -pois_log_like(params, model.x, model.y, id_indices)
    result = optimize(llike_optim, ones(eltype(model.x), size(model.x, 2)), LBFGS(); autodiff = :forward)
    return result.minimizer
end

# llike_optim(params) = -pois_log_like(params, x_array, y_array, id_indices)
# @time result2 = optimize(llike_optim, [1.5, 1.5], LBFGS(); autodiff = :forward);