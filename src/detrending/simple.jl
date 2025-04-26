import Polynomials

"""
    PolyNomialDetrending(degree::Int = 2)

Decompose timeseries `s` into a **sum** `x + r` where `x` is the trend
and `r` the residual, utilizing a polynomial guaranteed to have given `degree`.
For `degree = 1` this is linear detrending (ordinary least squares).
"""
@kwdef struct PonynomialDetrending <: Decomposition
    degree::Int = 1
end

function decompose(t, s, method::PonynomialDetrending)
    p = Polynomials.fit(t, s, method.degree)
    trend = p.(t)
    return trend, trend .- s
end

@kwdef struct MovingAverage
    window
end

struct NoDecomposition <: Decomposition end

function decompose(t, s, ::NoDecomposition)
    return s, similar(s)
end
