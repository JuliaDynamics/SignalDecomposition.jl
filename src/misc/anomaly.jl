using Dates
export TimeAnomaly

daymonth(t) = day(t), month(t)

"""
    TimeAnomaly()
Decompose timeseries `s` into its temporal average `x` and the anomalies `r`
so that `s = x + r`. Each unique day+month combination in `t` is identified, and the
values of `s` for each year that has this day+month combination are averaged.
As a result, the time vector `t` must be `<:AbstractVector{<:TimeType}`.

This method is very common in climate science, referred to as simply "anomalies".
"""
struct TimeAnomaly end

decompose(t, s, method::TimeAnomaly) = error("`t` must be `<:AbstractVector{<:TimeType}`")

function decompose(t::AbstractVector{<:TimeType}, s, method::TimeAnomaly)
    @assert length(t) == length(s)
    udates = unique!(sort!(daymonth.(t)))
    c = Dict(u => 0 for u in udates)
    v = Dict(u => zero(eltype(s)) for u in udates)
    x = copy(s)
    # First run: calculate climatology averages
    @inbounds for i in 1:length(t)
        date = daymonth(t[i])
        c[date] += 1
        v[date] += s[i]
    end
    # Second run: make x
    for (key, n) in c
        @inbounds v[key] = v[key]/n # take the average
    end
    @inbounds for i in 1:length(t)
        date = daymonth(t[i])
        x[i] = v[date]
    end
    return x, s .- x
end
