using LPVSpectral, Statistics, Random
export Sinusoidal

"""
    Sinusoidal(fs, x=true)
Decompose a timeseries `s` into a **sum** `x + r`, where `x`
are sinusoidal components with the given frequencies `fs` that minimize
coefficients ``A_i, \\phi_i`` of the expression
```math
s - \\bar{s} \\approx \\sum_i A_i \\cos(2\\pi f_i t + \\phi_i)
```
with ``\\bar{s}`` the mean.
There is arbitrarity of which part of the signal `x, r`
gets the mean value of `s`, because it is deducted for a better fit.
The argument `x=true` attributes it to `x`, use `false` for `r`.

This method uses a new least-squares algorithm in frequency domain using the package
[LPVSpectral.jl](https://github.com/baggepinnen/LPVSpectral.jl), see[^Bagge2017].

This method works for non-equispaced `t` axis (and also normal), is generally very accurate
(if choosen frequencies are not too close), but has performance scaling of
O(N^2.4) instead of O(n log(n)) of [`Fourier`](@ref).

[^Bagge2018]: F. Bagge Carlson et al., [Linear Parameter-Varying Spectral Decomposition](https://lup.lub.lu.se/search/publication/ac32368e-e199-44ff-b76a-36668ac7d595), 2017.
"""
struct Sinusoidal{F<:AbstractArray{<:AbstractFloat}} <: Decomposition
    fs::F
    x::Bool
end
Sinusoidal(fs) = Sinusoidal(fs, true)

function decompose(t, s, method::Sinusoidal)
    fs = method.fs
    μ = mean(s)
    χ, = ls_spectral(s .- μ, t, fs)
    As = abs.(χ) ./ sqrt(2*length(fs))
    φs = angle.(χ)
    x = zero(s)
    for i in 1:length(fs)
        @. x += As[i]*cos(2π*fs[i]*t + φs[i])
    end
    if method.x
        return x .+ μ, s .- x .- μ
    else
        return x, s .- x
    end
end
