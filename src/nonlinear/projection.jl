# TODO

# Method name: locally projective nonlinear noise reduction.

# https://journals.aps.org/pre/abstract/10.1103/PhysRevE.53.R4326
# and references therein
using NearestNeighbors
using DelayEmbeddings
using LinearAlgebra # for eigenvalues

# inputs
s = lorenzx + 0.1randn(length(lorenzx))
m = 40
w = 1 # theiler window
metric = NearestNeighbors.Euclidean()
ntype = 50 # neighborhood type
r = 1000.0 # r is large says the paper. "penaltizing factor"
Q = m÷2 # the difference of `m` to 2 times the capacity dimension of the chaotic set
# underlying the dynamics of input timeseries

# start algorithm
@assert Q < m
x = s # the paper uses the symbol `y` for the input timeseries
X = embed(x, m+1, 1) # capital y is paper's bold y
tree = KDTree(Y.data, metric)
R = ones(eltype(y), m+1)
R[1] = R[end] = r # something "large"
Γ = C = zeros(eltype(η), m+1, m+1) # notice: C and Γ are the same entity
𝒟 = zeros(eltype(x), Q, Q)
ξ = zeros(eltype(x), m+1)

# notice: the paper uses backwards embedding, but in DelayEmbeddings we traditionally
# use forward embedding by default. Thus e.g. equations (1) and (2) change a little
# but we are doing the same operations in the end.

for n in 1:length(y) # start process for every point in the embedded space

𝒰n = knn(tree, X[n], neighborhood)
K = length(𝒰n)

for i in 0:m; ξ[i+1] = (1/K)*sum(x[k+i] for k in 𝒰n); end

for j in 1:m+1, i in 1:m+1
    C[i, j] = @fastmath (1/K)*sum(x[k+i]*x[k+j] - ξ[i]*ξ[j] for k in 𝒰n)
    Γ[i, j] = @fastmath R[i]*C[i, j]*R[j] # make Γ
end
eig = eigenvec(Γ) # these are sorted by eigenvalue
for j in 1:m+1, i in 1:m+1
    𝒟[i, j] = sum(eig[i, q]*eig[j, q] for q in 1:Q)
end

# corrections vectors:
for i in 1:m+1
    θ[n][i] = (1/R[i]) * sum(𝒟[i, j] * R[j] * (ξ[j] - x[n+j-1]) for j in 1:m+1)
end


# TODO: add @inbounds everywhere once it is done
