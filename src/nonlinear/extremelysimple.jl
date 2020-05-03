#=
Extremely simple method from paper:

https://journals.aps.org/pre/abstract/10.1103/PhysRevE.47.2401

where a delay embedding is done, and the corrected timeseries is in fact just
an average of nearest neigbhors in the embedded space
=#
using Neighborhood
using DelayEmbeddings
export ExtremelySimpleNL

"""
    ExtremelySimpleNL(k::Int, ℓ::Int, w::Int, τ::Int, ε::Real) <: Decomposition
Quite literally the "extremely simple nonlinear noise-reduction method"[^Schreiber1993].
It decomposes `s` into the **sum** `x + r` with `x` being the "noiseless" timeseries.
This is the average position of neighbors in the delay embedded space.

This method works well if your timeseries is composed by the
addition of a structured component (which follows deterministic and
stationary dynamics which the embedding should approximate) and some noise.

The cited paper has some info on choosing optimal `ε`.

Arguments:
* `k` amount of past delay
* `ℓ` amount of forward delay
* `w` Theiler window
* `τ` delay time
* `ε` radius of the neighborhood in the embedded space

[^Schreiber1993]: [Schreiber, (1993) Extremely simple nonlinear noise-reduction method. Physical Review E, 47(4)](https://doi.org/10.1103/PhysRevE.47.2401)
"""
struct ExtremelySimpleNL <: Decomposition
    k::Int
    ℓ::Int
    w::Int
    τ::Int
    ε::Float64
end

function decompose(t, s, method::ExtremelySimpleNL)
    k, ℓ, w, τ, ε = getproperty.(Ref(method), (:k, :ℓ, :w, :τ, :ε))
    delays = [i*τ for i in -k:1:ℓ]
    embedding = GeneralizedEmbedding(Tuple(delays))
    X = genembed(s, Tuple(delays))
    tree = searchstructure(KDTree, X.data, Chebyshev())
    c = copy(s) # corrected timeseries
    theiler = Theiler(w)

    vec_of_idxs = bulkisearch(tree, X.data, WithinRange(ε), theiler)
    zc = 0

    for (i, n) in enumerate(τrange(s, embedding))
        𝒰n = vec_of_idxs[i]
        K = length(𝒰n)
        if K == 0
            zc += 1
            continue
        end
        c[n] = (1/K) * sum(s[j+k] for j in 𝒰n)
    end
    zc > 0 && @warn "$(zc) points were not corrected because no neighbors where found."
    return c, s .- c
end
