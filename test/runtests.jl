using SignalDecomposition
cd(@__DIR__)
using DelimitedFiles, Test, Random
Random.seed!(14141)

lorenzx = vec(readdlm("lorenzx.txt"))
roeslerz = vec(readdlm("roeslerz.txt"))

noise = randn(length(lorenzx))
residuals = (lorenzx, roeslerz, noise)
resnames = ("lorenz", "roessler", "gaussian")

# trivial structure component:
t = range(0, 24π, length = length(lorenzx))
structure = @. 7.2cos(t) + 5.6cos(2t + 3π/4)
twopiidx = findfirst(x -> x  ≥ 2π, t)
periods = [twopiidx, twopiidx/2]

m1 = Fourier(1 ./ periods)
m2 = FrequencySplit(maximum(1 ./ periods))

methods = (m1, m2)
m = m2
close("all")
for (name, re) in zip(resnames, residuals)
    s = structure .+ re
    for m in methods
        println(name*" with "*string(nameof(typeof(m))))
        x, r = decompose(s, m)

        @show nrmse(x, structure)
        @test nrmse(x, structure) < 0.1
        @show nrmse(r, re)
        @test nrmse(r, re) < 0.5
        # nice plots
        # figure()
        # title(string(m))
        # ax1 = subplot(211)
        # plot(structure; alpha = 0.75, label = "original")
        # plot(x; alpha = 0.75, ls = "dashed", label = "decomposed")
        # ylabel("structure"); legend()
        # subplot(212; sharex = ax1)
        # plot(re; alpha = 0.75, label = "original")
        # plot(r; alpha = 0.75, ls = "dashed", label = "decomposed")
        # ylabel("noise"); legend()
    end
end
