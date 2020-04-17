cd(@__DIR__)
using SignalDecomposition
using DelimitedFiles, Test, Random

# load timeseries
lorenzx = vec(readdlm("lorenzx.txt"))
roeslerz = vec(readdlm("roeslerz.txt"))
lorenzx_slow = vec(readdlm("lorenzx_slow.txt"))

noise = randn(Random.MersenneTwister(12441), length(lorenzx))
# trivial periodic component:
t = range(0, 24π, length = length(lorenzx))
periodic = @. 4 + 7.2cos(t) + 5.6cos(2t + 3π/5)
twopiidx = findfirst(x -> x  ≥ 2π, t)
tperiods = [twopiidx, twopiidx/2]

include("periodic_test.jl")
include("product_test.jl")
