using Statistics
using Optim

function YJ(y::Vector, λ)::Vector
## Purpose: Performs the Yeo-Johnson power transformation
## Input: a vector and λ
## Output: a vector
    for i in 1:length(y)
        if y[i] ≥ 0
            if λ ≠ 0
                y[i] = ((y[i] + 1)^λ - 1)/λ
            else
                y[i] = log(y[i] + 1)
            end
        else
            if λ ≠ 2
                y[i] = -((-y[i] + 1)^(2 - λ) - 1)/(2 - λ)
            else
                y[i] = -log(-y[i] + 1)
            end
        end
    end
    return y
end

function yeojohnson(y::Vector; min = -2, max = 2, λ = 0.1, opt = true)::Vector
## Purpose: Calls the YJ function depending on whether lambda should be optimised or not
## Input: A Vector. Optional inputs are the min and max search range for the optimiser, an optional λ parameter, and a flag for whether to optimise or use the λ input provided
## Output: a Vector of Yeo-Johnson transformed values
    if opt == true
        λ = λoptimum(y, min, max)
        return YJ(y, λ) 
    else
        return YJ(y, λ)
    end
end

function LogLike(y, λ)::Float32
## Purpose: Computes the log-likelihood of the Yeo-Johnson power transformation
## Input: a Vector and λ parameter
## Output: a Float
## Source: Algorithm from SciPy.stats `yeojohnson_lif` function in _morestats.py
    N = length(y)
    σ̂² = var(YJ(y, λ))
    LL = -N/2 * log(σ̂²) + (λ - 1) * sum(sgn.(y) .* log.(abs.(y) .+ 1))
    return LL
end

function sgn(x)
## Purpose: Calculates the "sign" of a number (1 if positive, -1 if negative, zero otherwise)
## Input: a number
## Output: 1, -1, or 0
    if x ≠ 0.0
        return x/abs(x)
    else
        return 0
    end
end

function λoptimum(y::Vector, min, max)::Float32
## Purpose: Computes the value of λ that minimises the log-likelihood of the Yeo-Johnson power transformation
## Input: a Vector and minimum and maximum values for the search range
## Output: a Value for λ̂, the minimiser
    LL(λ) = LogLike(y, λ)
    λ̂ = optimize(LL, min, max).minimizer
    return λ̂
end

