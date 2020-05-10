@testset "ExtremelySimpleNL" begin
# input
s1 = lorenzx + 0.1noise # input
s2 = roesslerz + 0.1noise # input
k = 5 # backward dim
ℓ = 5 # forward dim
ε = 0.3 # range
w = 1 # theiler window


for (name, s) in zip(("lorenz", "roessler"), (s1, s2))

    τ = estimate_delay(s, "mi_min") #  5 # delaytime
    method = ExtremelySimpleNL(k, ℓ, τ, w, ε)

    x, r = decompose(s, method)

    err = nrmse(x, s)
    @test err < 0.2

    # figure()
    # title(name)
    # plot(s; alpha = 0.75, label = "original")
    # plot(x; alpha = 0.75, ls = "dashed", label = "noiseless, err=$err")
    # legend()
end
end

@testset "ManifoldProjection" begin
    @testset "Henon" begin
        he = Systems.henon()
        tr = trajectory(he, 10000; Ttr = 100)
        Random.seed!(151521)
        z = tr[:, 1]
        s = z .+ randn(10001)*0.1*std(z)
        # s = lorenzx + 0.1noise
        m = 5
        w = 1 # theiler window
        metric = Euclidean()
        searchtype = 30
        Q = 2
        Q = [2, 2, 2, 3, 3, 3, 3]

        x, r = decompose(s, ManifoldProjection(m, Q, searchtype))

        err = nrmse(x, tr[:, 1])
        @test err < 0.1

        # figure()
        # plot(s, label = "input")
        # plot(z, color = "k", ls = "dashed", label = "real")
        # plot(x, alpha = 0.5, label = "output")

        # figure()
        # scatter(s[1:end-1], s[2:end], s = 0.1, label = "input")
        # scatter(z[1:end-1], z[2:end], s = 0.1, label = "real", color = "k")
        # scatter(x[1:end-1], x[2:end], s = 0.1, label = "output", alpha = 0.5)
        # legend()
    end

    @testset "lorenz" begin
        lo = Systems.lorenz()
        tr = trajectory(lo, 1000; Ttr = 100, dt = 0.1)
        Random.seed!(151521)
        z = tr[:, 1]
        s = z .+ randn(10001)*0.2*std(z)
        m = 13
        w = 1 # theiler window
        metric = Euclidean()
        k = 30
        Q = 2
        Q = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]

        x, r = decompose(s, ManifoldProjection(m, Q, 30))
        err = nrmse(x, z)
        @test err < 0.2

        # figure()
        # plot(s, label = "input")
        # plot(z, color = "k", ls = "dashed", label = "real")
        # plot(x, alpha = 0.5, label = "output")
        # estimate_delay(z, "mi_min")
        # xlim(0, 200)
        #
        # τ = 1
        # d = 2
        # figure()
        # plot(columns(embed(s, d, τ))..., marker = "o", ms = 1, lw = 0.1, label = "input")
        # plot(columns(embed(z, d, τ))..., marker = "o", ms = 1, lw = 0.1, label = "real", color = "k")
        # plot(columns(embed(x, d, τ))..., marker = "o", ms = 1, lw = 0.1, label = "real", alpha = 0.5)
    end
end
