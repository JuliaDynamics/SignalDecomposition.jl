using Dates
export TimeAnomaly

daymonth(t) = day(t), month(t)

"""
    TimeAnomaly()
Decompose timeseries `s` into its temporal average `x` and the anomalies `r`
so that `s = x + r`. Each unique day+month combination in `t` is identified, and the
values of `s` for each combination are averaged
(no special handling for February is and leap years).
As a result, the time vector `t` must be `<:AbstractVector{<:TimeType}`.

This method is very common in climate science, referred to as simply "anomalies".
"""
struct TimeAnomaly end

decompose(t, s, method::TimeAnomaly) = error("`t` must be `<:AbstractVector{<:TimeType}`")

function decompose(t::AbstractVector{<:TimeType}, s, method::TimeAnomaly)
    @assert length(t) == length(s)
    udates = unique!(sort!(daymonth.(t)))
    counters = Dict(u => [0, zero(eltype(s))] for u in udates)
    x = copy(s)
    # First run: calculate climatology averages
    @inbounds for i in 1:length(t)
        date = daymonth(t[i])
        counters[date] .+= (1, s[i])
    end
    # Second run: make x
    for (key, val) in counters
        @inbounds val[2] = val[2]/val[1] # take the average
    end
    @inbounds for i in 1:length(t)
        date = daymonth(t[i])
        x[i] = counters[date][2]
    end
    return x, s .- x
end
