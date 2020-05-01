cd(@__DIR__)
using SignalDecomposition
using DelimitedFiles, Test, Random, Statistics

# %% Generate timeseries
using DynamicalSystemsBase

ds = Systems.lorenz()
tr = trajectory(ds, 500.0; dt = 0.05, Ttr = 100)
lorenzx = tr[:, 1]/std(tr[:, 1])

tr = trajectory(ds, 20; dt = 0.002, Ttr = 100)
lorenzx_slow = tr[:, 1]/std(tr[:, 1])

ds = Systems.roessler()
tr = trajectory(ds, 500.0, dt = 0.05, Ttr = 10)
roesslerz = tr[:, 3]/std(tr[:, 3])

noise = randn(Random.MersenneTwister(12441), length(lorenzx))

# trivial periodic component:
t = range(0, 24π, length = length(lorenzx))
periodic = @. 4 + 7.2cos(t) + 5.6cos(2t + 3π/5)
twopiidx = findfirst(x -> x  ≥ 2π, t)
tperiods = [twopiidx, twopiidx/2]

# %% Run tests
include("periodic_test.jl")
include("product_test.jl")
include("nonlinear_test.jl")
