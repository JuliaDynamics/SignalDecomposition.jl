module SignalDecomposition

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end SignalDecomposition

export decompose, Decomposition

"Supertype of all decomposition methods."
abstract type Decomposition end

"""
    decompose([t, ] s, method::Decomposition) → x, r

Decompose an 1D input signal or timeseries `s(t)` into components `x, r`
using the given `method`. `t` defaults to `eachindex(s)`.

What are `x` and `r` are, and how they combine to give `s`, depends on `method`.
See the online documentation for all subtypes of `Decomposition`.
"""
decompose(s::AbstractVector, method::Decomposition; kwargs...) =
decompose(eachindex(s), s, method; kwargs...)

using Statistics
include("utils.jl")
include("linear/fourier.jl")
include("linear/lpv.jl")
include("product/matrixinversion.jl")
include("nonlinear/extremelysimple.jl")
include("nonlinear/projection.jl")
include("misc/anomaly.jl")
include("detrending/simple.jl")
include("detrending/smoothing.jl")



end # module
