# input
s = lorenzx + 0.1randn(length(lorenzx)) # input
k = 5 # backward dim
ℓ = 5 # forward dim
τ = 1 # delaytime
w = 1 # theiler window
ε = 0.3 # range

method = ExtremelySimpleNL(k, ℓ, τ, w, ε)

x, r = decompose(s, method)

@test nrmse(x, s) < 0.1
