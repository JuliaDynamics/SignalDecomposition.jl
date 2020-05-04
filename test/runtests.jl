cd(@__DIR__)
using SignalDecomposition
using DelimitedFiles, Test, Random, Statistics
using DynamicalSystemsBase

ds = Systems.lorenz()
tr = trajectory(ds, 500.0; dt = 0.05, Ttr = 100)
lorenzx = tr[:, 1]/std(tr[:, 1])

tr = trajectory(ds, 20; dt = 0.002, Ttr = 100)
lorenzx_slow = tr[:, 1]/std(tr[:, 1])

ds = Systems.roessler()
tr = trajectory(ds, 500.0, dt = 0.05, Ttr = 10)
roesslerz = tr[:, 3]/std(tr[:, 3])

L = length(lorenzx)

noise = randn(Random.MersenneTwister(12441), L)

# trivial periodic component:
te = range(0, 24π, length = length(lorenzx))
periodicf(t) = @. 4 + 7.2cos(t) + 5.6cos(2t + 3π/5)
periodic = periodicf(te)
twopiidx = findfirst(x -> x  ≥ 2π, te)
tperiods = [twopiidx, twopiidx/2]

# t = cumsum(rand(1000)/2)

# %% Run tests
include("periodic_test.jl")
include("lpv_test.jl")
include("product_test.jl")
include("nonlinear_test.jl")
