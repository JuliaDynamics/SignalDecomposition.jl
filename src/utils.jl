function isequispaced(t)
    d1 = t[2] - t[1]
    for i in 3:length(t)
        d = t[i] - t[i-1]
        d == d1 || return false
    end
    return true
end

"""
    findnearest(val, A)
Return the index of `A` which has value nearest to `val`.
"""
function findnearest(val::Real, A::AbstractVector{<:Real})
    i = 1
    d = abs(val - A[i])
    @inbounds for j in 1:length(A)
        dd = abs(val - A[j])
        if dd < d
            i = j
            d = dd
        end
    end
    return i
end


export nrmse, rmse
function mse(x, y)
    m = length(x)
    @assert m == length(y)
    @inbounds mse = sum(abs2(x[i] - y[i]) for i in 1:m) / m
    return mse
end

using Statistics: mean

"""
    rmse(x, y) → e
Return the root mean square error `e` of the "fit" `y` into data `x`.
"""
rmse(x, y) = sqrt(mse(x, y))

"""
    nrmse(x, y) → e
Return the normalized root mean square error of the "fit" `y` into data `x`.
This number is the relative error of `y` to `x` versus `mean(x)` to `x`, i.e.
if `e < 1` the fit `y` is better than using `mean(x)` as a fit.
"""
function nrmse(x, y)
    m = length(x)
    mean_out = mean(x)
    _mse     = mse(x, y)
    @inbounds msemean = sum(abs2(x[i] - mean_out) for i in 1:m) / m
    nrmse   = sqrt(_mse / msemean)
    return nrmse
end

export nrmse, rmse
