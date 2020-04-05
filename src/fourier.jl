using FFTW
export Fourier

# TODO: Allow width to the removed frequencies, ie. remove
# a Gaussian around the specified frequency
# This is possibly done via a combination of the
# DSP.Filters.Bandstop and probably some windowing function from DSP.

"""
    Fourier([s, ] frequencies) <: Decomposition
Simple decomposition method that identifies specific frequencies
at the Fourier space and removes them from the signal. The residual is the inverse of the
remaining Fourier signal, while the periodic part is just the removed components.
If a given frequency is not exactly matching the Fourier frequencies, the closest one is
removed.

If you provide `s` the method plans the forward and inverse
fourier transforms (so that it is efficient to re-use it for `s` of same type and length).

**Important**: periods/frequencies are defined with respect to the time axis length,
the actual time axis `t` is not used in this method. So, frequency `1/12` (a period of `12`)
means `12` time points (whose actual time value depends on `t`).

This method works well when a periodic signal P is superimposed on fluctuations S,
and you have a good educated guess of what frequencies compose P.
This method works if the given signal has length multiple of the largest period given.
"""
struct Fourier{F, I} <: Decomposition
    fs::Vector{Float64}
    forward::F
    inverse::I
end
Fourier(fs) = Fourier(fs, nothing, nothing)
function Fourier(s, fs)
    forward = plan_rfft(s)
    inverse = plan_irfft(forward*s, length(s))
    return Fourier(fs, forward, inverse)
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
    return periodic, residual
end
