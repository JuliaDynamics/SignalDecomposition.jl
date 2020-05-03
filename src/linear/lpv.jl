using LPVSpectral, Statistics, PyPlot, Random
export Sinusoidal

"""
    Sinusoidal(fs)
Decompose a timeseries `s` into a **sum** `x + r`, where `x`
are sinusoidal components with the given frequencies `fs` that minimize
coefficients ``A_i, \\phi_i`` of the expression
```math
s - \\bar{s} \\approx \\sum_i A_i \\cos(f_i t + \\phi_i)
```
with ``\\bar{s}`` the mean.
There is arbitrarity of which part of the signal `x, r`
gets the mean value (we just attribute it to `x`).

This method uses a new least-squares algorithm in frequency domain using the package
[LPVSpectral.jl](https://github.com/baggepinnen/LPVSpectral.jl), see[^Bagge2018] (Ch. 10)
for more details on the algorithm.

This method works for non-equispaced `t` axis (and also normal) and is only
2x slower than [`Fourier`](@ref) while not having the limiatations of `Fourier`.

[^Bagge2018]: Bagge Carlson, Fredrik, [PhD Thesis 2018, Lund University](https://lup.lub.lu.se/search/ws/files/55399444/phdthesis.pdf).
"""
struct Sinusoidal{F}
    fs::F
end

function decompose(t, s, method::Sinusoidal)
    fs = method.fs
    ωs = fs ./ 2π
    μ = mean(s)
    χ, = ls_spectral(s .- μ, t, ωs)
    As = abs.(χ) ./ sqrt(2*length(ωs))
    φs = angle.(χ)
    x = zero(s)
    for i in 1:length(fs)
        @. x += As[i]*cos(fs[i]*t + φs[i])
    end
    return x .+ μ, s .- x .- μ
end
