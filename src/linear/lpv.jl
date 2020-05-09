using LPVSpectral, Statistics, Random
export Sinusoidal

"""
    Sinusoidal(fs)
Decompose a timeseries `s` into a **sum** `x + r`, where `x`
are sinusoidal components with the given frequencies `fs` that minimize
coefficients ``A, \\phi`` of the expression
```math
s \\approx A_0 + \\sum_i A_i \\cos(2\\pi f_i t + \\phi_i)
```
with ``\\bar{s}`` the mean.

This method uses a new least-squares algorithm in frequency domain using the package
[LPVSpectral.jl](https://github.com/baggepinnen/LPVSpectral.jl), see[^Bagge2017].
It works for non-equispaced `t` axis (and also normal), is generally very accurate
(if choosen frequencies are not too close), but has performance scaling of
O(N^2.4) instead of O(n log(n)) of [`Fourier`](@ref).

Because it can work with arbitrary signal length the method always estimates the
zero-frequency Fourier component, and attributes it to `x`.
The fitted coefficients ``A, \\phi`` are available as fields `.A` and `.φ` of the struct
(first entry is zero-frequency component, i.e. the mean with respect to the sinusoidals).

[^Bagge2017]: F. Bagge Carlson et al., [Linear Parameter-Varying Spectral Decomposition](https://lup.lub.lu.se/search/publication/ac32368e-e199-44ff-b76a-36668ac7d595).
"""
struct Sinusoidal{F<:AbstractArray{<:AbstractFloat}} <: Decomposition
    fs::F
    A::Vector{Float64}
    φ::Vector{Float64}
end
Sinusoidal(fs) = Sinusoidal(fs, zeros(length(fs)+1), zeros(length(fs)+1))

function decompose(t, s, method::Sinusoidal)
    @assert 0 ∉ method.fs
    fs = [0, method.fs...]
    χ, = ls_spectral(s, t, fs)
    As = abs.(χ) ./ sqrt(2*length(fs))
    φs = angle.(χ)
    method.A .= As
    method.φ .= φs
    x = fill(As[1], length(s))
    for i in 2:length(fs)
        @. x += As[i]*cos(2π*fs[i]*t + φs[i])
    end
    return x, s .- x
end
