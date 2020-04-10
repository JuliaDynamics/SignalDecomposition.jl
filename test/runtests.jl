using SignalDecomposition
cd(@__DIR__)
using DelimitedFiles, Test, Random
Random.seed!(14141)

lorenzx = vec(readdlm("lorenzx.txt"))
roeslerz = vec(readdlm("roeslerz.txt"))
noise = randn(length(lorenzx))
residuals = (lorenzx, roeslerz, noise)

# trivial structure component:
t = range(0, 24π, length = length(lorenzx))
structure = @. 7.2cos(t) + 5.6cos(2t + 3π/4)
twopiidx = findfirst(x -> x  ≥ 2π, t)
periods = [twopiidx, twopiidx/2]

method = Fourier(1 ./ periods)

for re in residuals
    s = structure .+ re
    x, r = decompose(s, method)

    @test nrmse(x, structure) < 0.1
    @test nrmse(r, re) < 0.5
    # nice plots
    # figure()
    # title(string(method))
    # ax1 = subplot(211)
    # plot(structure; alpha = 0.75, label = "original")
    # plot(x; alpha = 0.75, ls = "dashed", label = "decomposed")
    # ylabel("structure"); legend()
    # subplot(212; sharex = ax1)
    # plot(re; alpha = 0.75, label = "original")
    # plot(r; alpha = 0.75, ls = "dashed", label = "decomposed")
    # ylabel("noise"); legend()
end
