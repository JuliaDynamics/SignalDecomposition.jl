using FFTW
export Fourier

# TODO: Allow width to the removed frequencies, ie. remove
# a Gaussian around the specified frequency
# This is possibly done via a combination of the
# DSP.Filters.Bandstop and probably some windowing function from DSP.

"""
    Fourier(frequencies) <: Decomposition
Simple decomposition method that identifies specific frequencies
at the Fourier space and removes them from the signal. The residual is the inverse of the
remaining Fourier signal, while the periodic part is just the removed components.
If a given frequency is exactly matching the Fourier frequencies, the closest one is
removed.

**Important**: periods/frequencies are defined with respect to the time axis length,
the actual time axis `t` is not used in this method. So, frequency `1/12` (a period of `12`)
means `12` time points (whose actual time value depends on `t`).

This method works well when a periodic signal P is superimposed on fluctuations S,
and you have a good educated guess of what frequencies compose P.
This method works if the given signal has length multiple of the largest period given.
"""
struct Fourier{F} <: Decomposition
    fs::F
end

function decompose(t, s, method::Fourier)
    isequispaced(t) || error("Input time axis must be equispaced")
    if length(s) % round(Int, maximum(1/f for f in method.fs)) ≠ 0
        @warn "The signal length is not a multiple of largest period."
    end

    m = mean(s)
    𝓕 = rfft(s .- m)
    fs = rfftfreq(length(s))
    for f in method.fs
        i = findnearest(f, fs)
        𝓕[i] = 0.0
    end
    inv_𝓕 = irfft(𝓕, length(s))
    residual = inv_𝓕
    periodic = s .- residual
    return periodic, residual
end
