module SignalDecomposition

export decompose, Decomposition

"Supertype of all decomposition methods."
abstract type Decomposition end

"""
    decompose([t, ] s, method::Decomposition) → x, r
Decompose an 1D input signal or timeseries `s(t)` into components, `x, r`,
using the given `method`. `t` defaults to `1:length(s)`.

What are `x` and `r` really depend on your point of view and your application.
They can be structure `x` and noise `r` (i.e. noise reduction). They can be
seasonal/periodic `x` and residual component `r`. They can even be multiplier `x` and
input `r`.
"""
decompose(s::AbstractVector, method::Decomposition; kwargs...) =
decompose(1:length(s), s, method; kwargs...)

using Statistics
include("utils.jl")
include("linear/fourier.jl")
include("linear/lpv.jl")
include("product/matrixinversion.jl")
include("nonlinear/extremelysimple.jl")
# include("nonlinear/projection.jl")

end # module
