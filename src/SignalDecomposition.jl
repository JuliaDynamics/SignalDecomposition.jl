module SignalDecomposition

export decompose, Decomposition

"Supertype of all decomposition methods."
abstract type Decomposition end

"""
    decompose([t, ] s, method::Decomposition) → x, r
Decompose an 1D input signal or timeseries `s(t)` into two components,
the structure `x` and noise `r`. Depending on your application, "structure" could be
a seasonal/periodic uninteresting component, and "noise" a residual interesting component.
The decomposition is done via dispatch on the `method` (see the methods for more).

Most of the time `x` also includes the mean value of `s`.
"""
decompose(s::AbstractVector, method::Decomposition; kwargs...) =
decompose(1:length(s), s, method; kwargs...)

using Statistics
include("utils.jl")
include("fourier.jl")
include("sinuisoidal.jl")
include("deconvolution.jl")
include("ssa.jl")
include("delay.jl")

end # module
