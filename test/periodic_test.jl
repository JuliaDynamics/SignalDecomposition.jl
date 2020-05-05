m1 = Fourier(periodic, 1 ./ tperiods)
m2 = FrequencySplit(maximum(1 ./ tperiods))
m3 = Sinusoidal([1, 2] ./ 2π)
mthods = (m1, m2, m3)

for m in mthods
    mstring = string(nameof(typeof(m)))
    @testset "Standard periodic, $mstring" begin
    for (name, re) in zip(("lorenz", "roessler", "gaussian"), (lorenzx, roesslerz, noise))
        s = periodic .+ re
        t = m isa Sinusoidal ? te : (1:length(s))
        x, r = decompose(t, s, m)
        errper = nrmse(periodic, x)
        @test errper < 0.1
        errres = nrmse(re, r)
        @test nrmse(r, re) < 0.5

        errori = nrmse(s, x .+ r)
        @test errori < 1e-15

        # println("  "*name*" errper=$errper, reserr=$reserr ")

        # figure()
        # title(string(m))
        # ax1 = subplot(211)
        # plot(periodic; alpha = 0.75, label = "original")
        # plot(x; alpha = 0.75, ls = "dashed", label = "err=$errper")
        # ylabel("x"); legend()
        # subplot(212; sharex = ax1)
        # plot(re; alpha = 0.75, label = "original")
        # plot(r; alpha = 0.75, ls = "dashed", label = "err=errres")
        # ylabel("r"); legend()
    end
    end
end
