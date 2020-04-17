m1 = Fourier(1 ./ tperiods)
m2 = FrequencySplit(maximum(1 ./ tperiods))
mthods = (m1, m2)

for m in mthods
    # println("Method: "string(nameof(typeof(m)))))
    for (name, re) in zip(("lorenz", "roessler", "gaussian"), (lorenzx, roeslerz, noise))
        s = periodic .+ re
        x, r = decompose(s, m)
        errper = nrmse(periodic, x)
        @test errper < 0.1
        errres = nrmse(re, r)
        @test nrmse(r, re) < 0.5
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
