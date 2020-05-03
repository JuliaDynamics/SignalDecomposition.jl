using SignalDecomposition, Random, PyPlot, Statistics
m = Sinusoidal([1, 2])
Random.seed!(151512)
tu = cumsum(rand(L)/4) # non-equispaced timeaxis
pu = periodicf(tu)

@testset "Sinusoidal" begin
    for (name, re) in zip(("cos(3t)", "gaussian"), (cos.(3tu), noise))
        s = pu .+ re .- mean(re)
        x, r = decompose(tu, s, m)
        errper = nrmse(pu, x)
        @test errper < 0.1
        errres = nrmse(re, r)
        @test nrmse(r, re) < 0.1
        # println("  "*name*" errper=$errper, errres=$errres")

        # figure()
        # ax1 = subplot(211)
        # plot(pu; alpha = 0.75, label = "x component")
        # plot(x; alpha = 0.75, ls = "dashed", label = "err=$errper")
        # ylabel("x"); legend()
        # subplot(212; sharex = ax1)
        # plot(re; alpha = 0.75, label = "r component")
        # plot(r; alpha = 0.75, ls = "dashed", label = "err=$errres")
        # ylabel("r"); legend()
    end
end
