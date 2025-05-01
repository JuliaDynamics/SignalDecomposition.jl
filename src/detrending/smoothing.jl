export MovingAverageSmoothing

"""
    MovingAverageSmoothing(window::Int = 10) <: Decomposition

TODO.
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