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
    @test err < 0.1

    # figure()
    # title(name)
    # plot(s; alpha = 0.75, label = "original")
    # plot(x; alpha = 0.75, ls = "dashed", label = "noiseless, err=$err")
    # legend()
end
