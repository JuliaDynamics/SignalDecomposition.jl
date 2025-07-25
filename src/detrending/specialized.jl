# Hodrick–Prescott filter obtained from
# https://github.com/sdBrinkmann/HPFilter.jl,
export HodrickPrescott
using SparseArrays, LinearAlgebra

"""
    HodrickPrescott(; λ = 1600, iter = 1)

Decompose timeseries `s` into a **sum** `x + r` where `x` is the trend
and `r` the residual according to the
[Hodrick-Prescott](https://en.wikipedia.org/wiki/Hodrick%E2%80%93Prescott_filter) filter.
The residual is called the "cyclic" component in this context.
The keyword `iter` controls the boosted (or iterative) version of the filter,
based on Peter Phillips and Zhentao Shi (2019): "Boosting the Hodrick-Prescott Filter".
"""
@kwdef struct HodrickPrescott <: Decomposition
    λ::Int = 1600
    iter::Int = 1
end

function decompose(t, x::AbstractVector, method::HodrickPrescott)
    (; λ, iter) = method
    n = length(x)
    m = 2
    @assert n > m
    I = Diagonal(ones(n))
    D = spdiagm(
        0 => fill(1, n-m),
        -1 => fill(-2, n-m),
        -2 => fill(1, n-m)
    )
    @inbounds D = D[1:n,1:n-m]
    S = (I + λ * D * D')
    function solve(S,x,iter)
        f = S \ x
        if iter > 1
            return solve(S,x-f,iter-1)
        else
            return x - f
        end
    end
    solution = solve(S, x, iter)
    trend = x - solve(S, x, iter)
    res = -solution
    return trend, res
end
