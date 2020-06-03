

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
    initial_params = ones(eltype(model.x), size(model.x, 2))
    llike_optim2 = TwiceDifferentiable(llike_optim, initial_params; autodiff=:forward)
    result = optimize(llike_optim2, initial_params)
    if !Optim.converged(result)
        @warn "Optimization failed to converge."
    end
    coefs = result.minimizer
    se = poisfe_se(model.y, model.x, model.id, coefs)
    # The following "non-robust" standard errors don't exactly match GLM.jl,
    # but they're close:
    # inv(Optim.hessian!(llike_optim2, result.minimizer))
    return coefs, se
end

"""
    poisfe_se(y, id, X, qcmle_coefs)

Returns robust standard errors.
"""
poisfe_se = function(y, X, id, qcmle_coefs)
    # No "hats" on variables are used for brevity.
    # We are already dealing with estimated quantities.
    β = qcmle_coefs # Estimated parameters, p.79
    K = length(β) # Number of parameters, p.79
    all_ids = unique(id)
    N = length(all_ids) # Number of groups, p.83

    # Initialize arrays described later
    𝐊 = zeros(Float64, K, K)
    𝐀 = zeros(Float64, K, K)
    𝐁 = zeros(Float64, K, K)
    𝚲 = Array{Array{Float64}}(undef, N)
    𝐖 = Array{Array{Float64}}(undef, N)
    𝐮 = Array{Array{Float64}}(undef, N)

    # Particular functional form (Poisson), p.79
    μ = function(𝐱ᵢₜ, β)
        return exp(𝐱ᵢₜ ⋅ β)
    end

    for i=1:N
        this_index = searchsorted(id, all_ids[i])
        T = length(this_index)
        yᵢ = y[this_index] # Tx1 vector of outcomes, p.79
        𝐱ᵢ = X[this_index, :] # TxK matrix of indep vars, p.79
        nᵢ = sum(yᵢ) # Scalar sum of outcomes, p.79

        𝛍ᵢ = similar(yᵢ) # Tx1 vector of conditional means, p.82
        for t in 1:T
            𝛍ᵢ[t] = μ(𝐱ᵢ[t, :], β)
        end
        Σ𝛍ᵢ = sum(𝛍ᵢ) # Scalar sum of conditional means, p.82

        𝐩 = 𝛍ᵢ / Σ𝛍ᵢ # Tx1 vector of "choice probabilities," p.79 (2.5)
        # 𝚲ᵢ (TxK) is the derivative of 𝐩 wrt β; p.86
        𝚲ᵢ = similar(𝐱ᵢ)
        for k in 1:K
            last_term = 0
            for s in 1:T
                last_term += 𝐱ᵢ[s, k] * 𝐩[s]
            end
            for t in 1:T
                𝚲ᵢ[t, k] = 𝐩[t] * (𝐱ᵢ[t, k] - last_term)
            end
        end

        𝐖ᵢ = Diagonal((1 ./ 𝐩))  # TxT, p.82
        𝐮ᵢ = yᵢ - nᵢ * 𝐩 # Tx1, p.82

        # Each group's 𝚲ᵢ, 𝐖ᵢ, and 𝐮ᵢ are stored in an array of arrays for
        # later computation of a hypothesis test, p.86
        𝚲[i] = 𝚲ᵢ
        𝐖[i] = 𝐖ᵢ
        𝐮[i] = 𝐮ᵢ

        𝐊 += 1/N * 𝚲ᵢ' * (nᵢ * 𝚲ᵢ) # p.86 (3.16)
        𝐀 += 1/N * nᵢ * (𝚲ᵢ' * 𝐖ᵢ * 𝚲ᵢ) # p.83 (3.9)
        𝐁 += 1/N * 𝚲ᵢ' * 𝐖ᵢ * 𝐮ᵢ * 𝐮ᵢ' * 𝐖ᵢ * 𝚲ᵢ # p.83 (3.10)
    end

    𝐀⁻¹ = inv(𝐀)
    se_rob = broadcast(sqrt, diag((𝐀⁻¹ * 𝐁 * 𝐀⁻¹) / N )) # p.83, (3.11)

    # Test of the conditional mean specification (3.1) with the null
    # hypothesis (3.14), p.85.
    # NOT YET IMPLEMENTED
    # 𝐫 = Array{Float64, 2}(undef, N, K) # p.86, (3.17)
    # for i=1:N
    #     𝐫[i,:] = 𝐮[i]' * (𝚲[i] - 𝐖[i] * 𝚲[i] * 𝐀⁻¹ * 𝐊')
    # end
    # ssr = sum((ones(N) - 𝐫 * (𝐫 \ ones(N))) .^ 2) # p.86 (3.18)
    # p_value = 1 - chisqcdf(K, N - ssr) # p.86
    return se_rob
end