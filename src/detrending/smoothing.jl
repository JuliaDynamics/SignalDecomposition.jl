export MovingAverageSmoothing, LoessSmoothing
import Loess

"""
    MovingAverageSmoothing(window::Int = 10) <: Decomposition

Decompose timeseries `s` into a **sum** `x + r` where `x` is the smoothened signal
and `r` the residual. The smoothing is done via a basic moving average using a fixed
rectangualr window of size `window`.
At the end and start of teh timeseries only half of the window can be used for averaging.
"""
@kwdef struct MovingAverageSmoothing <: Decomposition
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

"""
    LoessSmoothing(; span, degree, cell) <: Decomposition

Decompose timeseries `s` into a **sum** `x + r` where `x` is the smoothened signal
and `r` the residual. The smoothing is done via locally estimated scatterplot smoothing
(loess) using Loess.jl and the provided keywords.

- `span`: The degree of smoothing, typically in [0,1]. Smaller values result in smaller
  local context in fitting.
- `degree`: Polynomial degree.
- `cell`: Control parameter for bucket size. Internal interpolation nodes will be
added to the K-D tree until the number of bucket element is below `n * cell * span`.
"""
@kwdef struct LoessSmoothing{S, C} <: Decomposition
    span::S = 0.75
    degree::Int = 2
    cell::C = 0.2
end

function decompose(t, s, method::LoessSmoothing)
    model = Loess.loess(t, s, span = method.span)
    trend = Loess.predict(model, t)
    return trend, s .- trend
end