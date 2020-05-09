using Neighborhood
using DelayEmbeddings
using LinearAlgebra # for eigenvalues
using StaticArrays

# TODO: Would be nice to be able to generalize this to have `τ` as delay time
# instead of forcing `τ=1`.

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

`k` can be either `Int` or a `SearchType` from [Neighborhood.jl](https://julianeighbors.github.io/Neighborhood.jl/stable/#Search-types-1).
If `Int`, the `k` nearest neighbors are choosen as the neighborhood 𝓤 of each point.
The paper contains an involved process for determining optimal `k::Int`, see eq.(5.4).
`w` is just the Theiler window, while `r` is the value of the edge entries of vector R,
and probably has not much impact (the rest of the entries are 1).

In the paper too big correction vectors were rescaled to the average magnitude of corrections,
using as a criterion the distribution of their size. This is not implemented here (as it
is not clear *exactly* what it means computationally, what is "too big"?)
Contributing it is welcomed if you know how...

See also [^Schreiber1996] for an application of the same algorithm in real ECG data.

[^Schreiber1996]: Schreiber & Kaplan (1996). Nonlinear noise reduction for electrocardiograms. [Chaos, 6(1), 87–92](https://doi.org/10.1063/1.166148)

[^Grassberger1993]: Grassberger et al., (1993). On noise reduction methods for chaotic data. Chaos 3(2), 127–141](https://doi.org/10.1063/1.165979)
"""
struct ManifoldProjection{ST <: SearchType} <: Decomposition
    m::Int
    Qs::Vector{Int}
    searchtype::ST
    w::Int
    r::Float64
end
ManifoldProjection(m, Q::Int, args...) = ManifoldProjection(m, [Q], args...)
ManifoldProjection(m, Q, k::Int, args...) = ManifoldProjection(m, Q, NeighborNumber(k), args...)
ManifoldProjection(m, Q, st::ST, w = 0, r = 1000.0) where {ST<:SearchType} =
ManifoldProjection{ST}(m, Q, st, w, r)

function SignalDecomposition.decompose(t, s, method::ManifoldProjection)
    m = method.m
    @assert maximum(method.Qs) < m
    x = copy(s) # x is "input timeseries", iteratively corrected
    ge = GeneralizedEmbedding(tuple(0:m...))

    # Function barrier to make `m` type parameter:
    iteratively_project!(x, method, ge)
    return x, s .- x
end

# TODO: Replace all 0:m and 1:m+1 with 1:D ?
# TODO: add inbounds everywhere

function iteratively_project!(x, method, ge::GeneralizedEmbedding{D}) where {D}
    m = D-1
    # Initialize a bunch of stuff
    r = copy(x)
    R = ones(eltype(x), m+1)
    R[1] = R[end] = method.r # something "large"
    Γ = C = zeros(eltype(x), m+1, m+1)
    𝒟 = copy(C)
    ξ = zeros(eltype(x), m+1)
    θ = zeros(eltype(x), m+1)
    Vc = genembed(x, ge)   # corrected embedded set

    # Main algorithm loop, repeated for Qs
    for Q in method.Qs
        V = genembed(x, ge) # embedded set, which is "bold x" in paper
        tree = searchstructure(KDTree, V.data, Euclidean())
        𝒰 = bulkisearch(tree, V.data, method.st, Theiler(method.w); sortds = false)

        for n in 1:length(V) # start process for every point in the embedded space
            𝒰n = 𝒰[n]
            K = length(𝒰n)
            iK = 1/K

            for i in 0:m; ξ[i+1] = iK*sum(x[k+i] for k in 𝒰n); end

            for j in 1:m+1, i in 1:m+1
                ξij = ξ[i]*ξ[j]
                C[i, j] = @fastmath iK*sum(x[k+i-1]*x[k+j-1] - ξij for k in 𝒰n)
                Γ[i, j] = @fastmath R[i]*C[i, j]*R[j] # make Γ
            end
            eig = eigvecs(Γ) # these are sorted by eigenvalue. Important for the algorithm.
            for j in 1:m+1, i in 1:m+1
                𝒟[i, j] = sum(eig[i, q]*eig[j, q] for q in 1:Q)
            end

            # corrections vectors:
            for i in 1:m+1
                θ[i] = (1/R[i]) * sum(𝒟[i, j] * R[j] * (ξ[j] - x[n+j-1]) for j in 1:m+1)
            end

            Vc.data[n] = V.data[n] + SVector{D}(θ)
        end
        # okay, now use corrected points to estimate new points, by taking the average.
        for n in (m+1):length(V)-m
            r[n] = 1/(m+1)*sum(Vc[n-j, j+1] for j in 0:m)
        end
        # Edge cases
        for n in 1:m
            r[n] = sum(Vc[n-j+1, j] for j in 1:n)/n
        end
        for (z, n) in enumerate((length(V)-m):-1:(length(V)-2m))
            r[n] = sum(Vc[n-j, m+1-j] for j in 0:(z-1))/z
        end
        x .= r # set corrected vector to input vector.
    end
end

# return noiseless, x .- noiseless
scatter(s[1:end-1], s[2:end], s = 1)
scatter(Vcorrected[1:end-1, 1], Vcorrected[2:end, 2], s = 1)
scatter(noiseless[1:end-1], noiseless[2:end], s = 1)
