#=
Extremely simple method from paper:

https://journals.aps.org/pre/abstract/10.1103/PhysRevE.47.2401

where a delay embedding is done, and the corrected timeseries is in fact just
an average of nearest neigbhors in the embedded space
=#
using NearestNeighbors
using DelayEmbeddings

# input
x = s
s = lorenzx + 0.1randn(length(lorenzx))
k = 5 # forward dim
ℓ = 5 # backward dim
τ = 1 # delaytime
w = 1 # theiler window
metric = NearestNeighbors.Chebyshev()
ε = 0.1 # range

delays = [-i*τ for i in k:-1:1]
push!(delays, 0)
append!(delays, [i*τ for i in 1:ℓ])
embedding = GeneralizedEmbedding(Tuple(delays))
X = genembed(x, Tuple(delays))
tree = KDTree(X.data, metric)
c = copy(x) # corrected timeseries

for n in τrange(embedding)
𝒰n = inrange(tree, X[n], neighborhood with ε and theiler window)[1] # idxs only
K = length(𝒰n)
c[n] = (1/K) * sum(x[j] for j in 𝒰n)
end

return c, x .- c
