# TODO: make it into a sctruct and extent `decompose`

# Method name: locally linear projections, nonlinear noise reduction.

# https://journals.aps.org/pre/abstract/10.1103/PhysRevE.53.R4326
# and references therein

# Algorithm "IV" of paper "On noise reduction methods for chaotic data".
# Further used in the following two refs:...

using Neighborhood
using DelayEmbeddings
using LinearAlgebra # for eigenvalues
using StaticArrays

# TODO: add @inbounds everywhere once it is done

# inputs
he = Systems.henon()
tr = trajectory(he, 10000; Ttr = 100)
s = tr[:, 1] .+ randn(10001)/15

# s = lorenzx + 0.1noise
m = 5
w = 0 # theiler window
metric = Euclidean()
r = 1000.0 # r is large says the paper. "penaltizing factor"
searchtype = NeighborNumber(30)
Q = 2 # the difference of `m` to 2 times the capacity dimension of the chaotic set is Q
# Q is "the dimension of the local manifold projected on"
# So Q is the dimension of the expected manifold dimension of the attractor of the
# noiseless timeseries

# start algorithm
@assert Q < m
V = embed(s, m+1, 1) # embedded space, which is "bold x" in paper
Vcorrected = deepcopy(V)
x = s # paper uses x throughout, simpler this way
noiseless = copy(s)

# Initialize a bunch of stuff
tree = searchstructure(KDTree, V.data, metric)
R = ones(eltype(s), m+1)
R[1] = R[end] = r # something "large"
Γ = C = zeros(eltype(s), m+1, m+1)
𝒟 = copy(C)
 # notice: C and Γ are the same entity, only for clarity of source I use 2 symbols
ξ = zeros(eltype(s), m+1)
θ = zeros(eltype(s), m+1) # correction vector

𝒰 = bulkisearch(tree, V.data, searchtype, Theiler(w)) # all neighbors for all points

for n in 1:length(V) # start process for every point in the embedded space
    𝒰n = 𝒰[n]
    K = length(𝒰n)
    # ensure found indices

    for i in 0:m; ξ[i+1] = (1/K)*sum(x[k+i] for k in 𝒰n); end

    for j in 1:m+1, i in 1:m+1
        C[i, j] = @fastmath (1/K)*sum(x[k+i-1]*x[k+j-1] - ξ[i]*ξ[j] for k in 𝒰n)
        Γ[i, j] = @fastmath R[i]*C[i, j]*R[j] # make Γ
    end
    eig = eigvecs(Γ) # these are sorted by eigenvalue. Important for the algorithm.
    for j in 1:m+1, i in 1:m+1
        𝒟[i, j] = sum(eig[i, q]*eig[j, q] for q in 1:Q)
    end

    # corrections vectors:
    θ = @SVector [
        (1/R[i]) * sum(𝒟[i, j] * R[j] * (ξ[j] - x[n+j-1]) for j in 1:m+1)
        for i in 1:m+1
    ]

    Vcorrected.data[n] = V.data[n] + θ
end

# The paper states that "every point in the timeseries" participates in "exactly"
# m+1 delay vectors. It is TERRIBLE that something so wrong can be said so casually.
# This is of course false, the start and end of the timeseries do not satisfy this

# okay, now use corrected points to estimate new points, by taking the average.

# This is incorrect here:

for n in (m+1):length(V)-m
    noiseless[n] = 1/(m+1)*sum(Vcorrected[n+j, j+1] for j in 0:m)
end

# Edge cases
for n in 1:m
    noiseless[n] = 1/(n)*sum(Vcorrected[n+j-1, j] for j in 1:n)
end
for (z, n) in enumerate((length(V)-m):-1:(length(V)-2m))
    noiseless[n] = 1/(z)*sum(Vcorrected[n-j, m+1-j] for j in 0:(z-1))
end

# return noiseless, x .- noiseless
scatter(s[1:end-1], s[2:end], s = 1)
scatter(Vcorrected[1:end-1, 1], Vcorrected[2:end, 2], s = 1)
scatter(noiseless[1:end-1], noiseless[2:end], s = 1)
