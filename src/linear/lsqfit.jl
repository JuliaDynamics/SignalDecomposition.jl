export SinusoidalFit
import LsqFit

"""
    SinusoidalFit(s, fs) <: Decomposition
Decompose a timeseries `s` into a **sum** `x + r` (with `x` the periodic component),
by fitting sinuisoidals with given
frequencies `fs` to the signal `s` using the package `LsqFit`. Specifically, fit
```math
s - \\bar{s} \\approx \\sum_i A_i \\sin(2\\pi f_i t + \\phi_i)
```
with ``\\bar{s}`` the mean. The fit happens on the amplitudes
and phases ``A_i, \\phi_i``. After the `decomposition` you can find these in the struct's
fields `amps, phases`. The fit is done on `s` versus `t`, so be sure that you have
transformed `t` appropriately (e.g. if `t` is "days" but your frequencies are multiples of
years, then you should give `t/365.26`).

    SinusoidalFit(fs, φ0s, A0s [, ub, lb])

The quality of the fit depends dramatically on the initial guesses for the phases
and amplitudes, `φ0s, A0s`. This second constructor
gives full control over initial phases and amplitudes, as well as
upper and lower bounds on the amplitudes `ub, lb` (also vectors of `length(fs)`).
In the first constructor `φ0 = 0, A0 = abs(-(extrema(s)...))/2, ub = Inf, lb = -Inf` for
all frequencies. (the bounds for the phases are always ± π)

**Notice**: `LsqFit` performs poorly for fitting sinusoidals and `Fourier`
should be preferred over this method if the signal given is in
multiples of the expected periods (and of course `t` is equally spaced).
"""
mutable struct SinusoidalFit{T<:AbstractFloat} <: Decomposition
    fs::Vector{T}
    φ0s::Vector{T}
    A0s::Vector{T}
    ub::Vector{T}
    lb::Vector{T}
    amps::Vector{T}
    phases::Vector{T}
end

function SinusoidalFit(s::AbstractVector{T}, fs) where {T}
    φ0s = zero(fs)
    A0s = fill(abs(-(extrema(s)...))/2, length(fs))
    ub = fill(T(Inf), length(fs))
    lb = fill(T(-Inf), length(fs))
    SinusoidalFit{T}(fs, φ0s, A0s, ub, lb, zero(fs), zero(fs))
end

SinusoidalFit(fs, φ0s, A0s, ub=fill(Inf, length(fs)), lb=fill(-Inf, length(fs))) =
SinusoidalFit{eltype(fs)}(fs, φ0s, A0s, ub, lb, zero(fs), zero(fs))

function decompose(t, s, method::SinusoidalFit{T}) where {T}
    fs = method.fs
    L = length(fs)
    function sinuisoidal(t, p)
        r = zeros(T, length(t))
        for i in 1:L
            @. r += p[i]*sin(2π * (fs[i]*t + p[i+L]))
        end
        return r
    end

    p0 = vcat(method.A0s, method.φ0s)
    ub = vcat(method.ub, fill( T(π), L))
    lb = vcat(method.lb, fill(-T(π), L))

    fit = LsqFit.curve_fit(sinuisoidal, t, s .- mean(s), p0; lower=lb, upper=ub)
    method.amps .= LsqFit.coef(fit)[1:L]
    method.phases .= LsqFit.coef(fit)[L+1:2L]
    periodic = sinuisoidal(t, LsqFit.coef(fit)) .+ mean(s)
    return periodic, s .- periodic
end
