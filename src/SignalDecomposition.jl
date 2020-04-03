module SignalDecomposition

export decompose, Decomposition

"Abstract type of all decomposition methods."
abstract type Decomposition end

"""
    decompose([t, ] s, method::Decomposition) → p, r
Decompose an 1D input signal or timeseries `s(t)` into its periodic and residual
components `p, r` using the given `method`. `t` defaults to `1:length(s)`.
Other common names for `p, r` are "seasonal" and "trend".

The decomposition is under addition, i.e. `s = p .+ r`, unless
explicitly noted otherwise by the `method`.
By convention, the periodic component also includes the mean value of `s`.
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
