# Examples
Only a few examples are shown here. Every method has an example (and plotting code) in the `test` folder!

## Nonlinear
```@example docs
using SignalDecomposition, DynamicalSystems, Random, Plots, Statistics

he = Systems.henon()
tr = trajectory(he, 10000; Ttr = 100)
Random.seed!(151521)
z = tr[:, 1]
s = z .+ randn(10001)*0.1*std(z)
m = 5
w = 1 # theiler window
metric = Euclidean()
k = 30
Q = [2, 2, 2, 3, 3, 3, 3]
x, r = decompose(s, ManifoldProjection(m, Q, k))
```

```@example docs
p1 = plot(s, label = "input")
plot!(p1, z, color = :black, ls = :dash, label = "real")
plot!(p1, x, alpha = 0.5, label = "output")
xlims!(p1, 0, 50)
p1
```

Alright, this doesn't seem much of a difference to be honest.
One sees a big difference once going into the state space and looking at the attractor:

```@example docs
p2 = scatter(s[1:end-1], s[2:end], ms = 1, label = "input", msw = 0)
scatter!(p2, z[1:end-1], z[2:end], ms = 1, label = "real", color = :black, msw = 0)
scatter!(p2, x[1:end-1], x[2:end], ms = 1, label = "output", alpha = 0.5, msw = 0)
p2
```

## Time and Sinusoidal
```@example docs
using SignalDecomposition, Dates, Random, Plots
Random.seed!(41516)
y = Date(2001):Day(1):Date(2025)
dy = dayofyear.(y)
cy =  @. 4 + 7.2cos(2π*dy/365.26) + 5.6cos(4π*dy/365.26 + 3π/5)
r0 = randn(length(dy))/2
sy = cy .+ r0

x, r = decompose(y, sy, TimeAnomaly())

t = collect(1:length(y)) ./ 365.26 # true time in years

x2, r2 = decompose(t, sy, Sinusoidal([1.0, 2.0]))

p3 = plot(t, sy, label = "input")
plot!(p3, t, cy, label = "true periodic", color = :black, ls = :dash)
plot!(p3, t, x, label = "TimeAnomaly", alpha = 1.0, color = :red)
plot!(p3, t, x2, label = "Sinusoidal", alpha = 0.5, color = :green)
xlabel!(p3, "years")
xlims!(p3, 0, 1) # zoom in
```

Allthough not immediatelly obvious from the figure, `Sinusoidal` performs better:
```@example docs
rmse(cy, x), rmse(cy, x2)
```
