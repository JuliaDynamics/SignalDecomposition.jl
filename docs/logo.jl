using SignalDecomposition, Dates, Random, Plots
using ColorSchemes, Colors

pyplot()

Random.seed!(41516)
y = Date(2001):Day(1):Date(2025)
dy = dayofyear.(y)
cy =  @. 4 + 7.2cos(2π*dy/365.26) + 5.6cos(4π*dy/365.26 + 3π/5)
r0 = randn(length(dy))/2
sy = cy .+ r0

t = collect(1:length(y)) ./ 365.26 # true time in years
x, r = decompose(t, sy, Sinusoidal([1.0, 2.0]))

is = collect(0:0.02:1.0)
s1 = x .+ 6
s2 = r .- 10

cg1 = ColorScheme([parse(RGB, "black"), parse(RGB, "#745bb1")])
cg2 = ColorScheme([parse(RGB, "black"), parse(RGB, "#1aa0a4")])

append!(is, ones(100))
p1 = plot(t, sy, label = "input", color = "black")
annotate!(p1, [(1, 22.5, text("SignalDecomposition.jl", 48; family = "Bahnschrift"))])
ylims!(p1, -12, 25)
xlims!(0, 2)
plot!(p1, t, s1, label = "decomposed", alpha = 1.0)
plot!(p1; size = (1200, 600), showaxis = false, grid = false, legend=false)

# %%
anim = @animate for i in is
    p1 = plot(t, sy, label = "input", color = "black")
    x1 = @. i*s1 + (1-i)*sy
    x2 = @. i*s2 + (1-i)*sy
    plot!(p1, t, x1, label = "decomposed", alpha = 1.0, color = get(cg1, i))
    plot!(p1, t, x2, color = get(cg2, i))
    xlims!(0, 2)
    annotate!(p1, [(1, 22.5, text("SignalDecomposition.jl", 48; family = "Bahnschrift"))])
    ylims!(p1, -12, 25)
    plot!(p1; size = (1200, 600), showaxis = false, grid = false, legend=false)
end
gif(anim, "signaldecomposition.gif", fps=30)
