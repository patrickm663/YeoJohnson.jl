using Statistics
using Optim

function YJ(y3::Vector, λ3)::Vector
## Purpose: Performs the Yeo-Johnson power transformation
## Input: a vector and λ
## Output: a vector
    size = length(y3)
    yTransformed = zeros(size)
    for i in 1:size
        if y3[i] ≥ 0
            if λ3 ≠ 0
                yTransformed[i] = ((y3[i] + 1)^λ3 - 1)/λ3
            else
                yTransformed[i] = log(y3[i] + 1)
            end
        else
            if λ3 ≠ 2
                yTransformed[i] = -((-y3[i] + 1)^(2 - λ3) - 1)/(2 - λ3)
            else
                yTransformed[i] = -log(-y3[i] + 1)
            end
        end
    end
    return yTransformed
end

function yeojohnson(y2::Vector; min = -2, max = 2, λ2 = 0.1, opt = true)::Vector
## Purpose: Calls the YJ function depending on whether lambda should be optimised or not
## Input: A Vector. Optional inputs are the min and max search range for the optimiser, an optional λ parameter, and a flag for whether to optimise or use the λ input provided
## Output: a Vector of Yeo-Johnson transformed values
    if opt == true
        λ2 = λoptimum(y2, min, max)
        return YJ(y2, λ2) 
    else
        return YJ(y2, λ2)
    end
end

function LogLike(y1, λ1)::Float32
## Purpose: Computes the log-likelihood of the Yeo-Johnson power transformation
## Input: a Vector and λ parameter
## Output: a Float
## Source: Algorithm from SciPy.stats `yeojohnson_lif` function in _morestats.py
    N = length(y1)
    σ̂² = var(YJ(y1, λ1))
    LL = -N/2 * log(σ̂²) + (λ1 - 1) * sum(sgn.(y1) .* log.(abs.(y1) .+ 1))
    return -LL
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

