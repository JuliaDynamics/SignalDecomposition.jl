export PolynomialDetrending
import Polynomials

"""
    PolynomialDetrending(degree::Int = 1) <: Decomposition

Decompose timeseries `s` into a **sum** `x + r` where `x` is the trend
and `r` the residual. The trend is a fitted polynomial guaranteed to have given `degree`.
For `degree = 1` this is linear detrending (ordinary least squares).
"""
@kwdef struct PolynomialDetrending <: Decomposition
    degree::Int = 1
end

function decompose(t, s, method::PolynomialDetrending)
    p = Polynomials.fit(t, s, method.degree)
    trend = p.(t)
    return trend, trend .- s
end

"""
    NoDecomposition <: Decomposition

Decompose timeseries `s` into `s` and zeros.
Dummy method provided simply so that one can loop over various decomposition
methods and no decomposition at all while still using the `decompose` function.
"""
struct NoDecomposition <: Decomposition end

function decompose(t, s, ::NoDecomposition)
    return s, zeros(eltype(s), length(s))
end
