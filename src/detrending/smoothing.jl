export MovingAverageSmoothing

"""
    MovingAverageSmoothing(window::Int = 10) <: Decomposition

Decompose timeseries `s` into a **sum** `x + r` where `x` is the smoothened singla
and `r` the residual. The smoothing is done via a basic moving average using a fixed
rectangualr window of size `window`.
At the end and start of teh timeseries only half of the window can be used for averaging.
"""
@kwdef struct MovingAverageSmoothing
    window::Int = 10
end

function decompose(t, s, mas::MovingAverageSmoothing)
    window = mas.window
    n = length(s)
    smooth = zeros(eltype(s), n)
    for i in 1:n
        half = (window - 1) ÷ 2
        start_idx = max(1, i - half)
        end_idx = min(n, i + half)
        smooth[i] = mean(@view s[start_idx:end_idx])
    end
    return smooth, s .- smooth
end