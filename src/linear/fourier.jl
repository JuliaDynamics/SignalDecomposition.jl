using FFTW
export Fourier, FrequencySplit

# TODO: Allow width to the removed frequencies, ie. remove
# a Gaussian around the specified frequency
# This is possibly done via a combination of the
# DSP.Filters.Bandstop and probably some windowing function from DSP.

"""
    Fourier([s, ] frequencies, x=true) <: Decomposition

Decompose a timeseries `s` into a **sum** `x + r`, by identifying specific `frequencies`
at the Fourier space and removing them from the signal.
`x` is the removed periodic component while `r` is the residual.
If a given frequency is not exactly matching the Fourier frequencies, the closest one is
removed.

**Important**: periods/frequencies are defined with respect to the `t` axis length,
the actual `t` values are not used in this method. So, frequency `1/12` (a period of `12`)
means `12` data points (whose actual value depends on `t`).

If you provide `s` the method plans the forward and inverse
Fourier transforms (so that it is efficient to re-use it for `s` of same type and length).

This method works well when a periodic signal P is superimposed on fluctuations S,
and you have a good educated guess of what frequencies compose P.
This method works well if the given signal has length multiple of the periods given.

There is arbitrarity of which part of the signal `x, r`
gets the mean value of `s`, because it is deducted for a better fit.
The argument `x=true` attributes it to `x`, use `false` for `r`.
"""
struct Fourier{F, I} <: Decomposition
    fs::Vector{Float64}
    forward::F
    inverse::I
    x::Bool
end
Fourier(fs, x::Bool=true) = Fourier(fs, nothing, nothing, x)
function Fourier(s::AbstractVector, fs::AbstractVector, x::Bool=true)
    forward = plan_rfft(s)
    inverse = plan_irfft(forward*s, length(s))
    return Fourier(fs, forward, inverse, x)
end

function decompose(t, s, method::Fourier)
    isequispaced(t) || error("Input time axis must be equispaced for method Fourier.")
    if length(s) % round(Int, maximum(1/f for f in method.fs)) ≠ 0
        @warn "The signal length is not a multiple of largest period."
    end

    m = mean(s)
    𝓕 = isnothing(method.forward) ? rfft(s .- m) : method.forward*(s .- m)
    fs = rfftfreq(length(s))
    for f in method.fs
        i = findnearest(f, fs)
        𝓕[i] = 0.0
    end
    inv_𝓕 = isnothing(method.inverse) ? irfft(𝓕, length(s)) : method.inverse*𝓕
    residual = inv_𝓕
    periodic = s .- residual
    if !method.x
        periodic .-= m
        residual .+= m
    end
    return periodic, residual
end


"""
    FrequencySplit([s, ] f::Real) <: Decomposition
Similar to the [`Fourier`](@ref) method, but now the "residual" signal is the part
of `s` with frequencies higher than `f`, while the "seasonal" part has frequencies `≤ f`.
"""
struct FrequencySplit{F, I} <: Decomposition
    f::Float64
    forward::F
    inverse::I
end
FrequencySplit(fs) = FrequencySplit(fs, nothing, nothing)
function FrequencySplit(s, fs)
    forward = plan_rfft(s)
    inverse = plan_irfft(forward*s, length(s))
    return FrequencySplit(fs, forward, inverse)
end

function decompose(t, s, method::FrequencySplit)
    isequispaced(t) || error("Input time axis must be equispaced for method FrequencySplit.")

    m = mean(s)
    𝓕 = isnothing(method.forward) ? rfft(s .- m) : method.forward*(s .- m)
    fs = rfftfreq(length(s))
    i = findlast(f -> f ≤ method.f, fs)
    𝓕[1:i+1] .= 0.0
    inv_𝓕 = isnothing(method.inverse) ? irfft(𝓕, length(s)) : method.inverse*𝓕
    residual = inv_𝓕
    periodic = s .- residual
    return periodic, residual
end
