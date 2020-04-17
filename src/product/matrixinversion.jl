using BandedMatrices, LinearAlgebra
export ProductInversion

"""
    ProductInversion(r, μ) <: Decomposition
Decompose a timeseries `s` into a **product** `x * r`, given that you have a good
estimate of one of the two factors `r` (the "input") and you need `x` (the "multiplier")
but you can't do simply `x =  s ./ r` because `r` contains zeros.

This method works well when the characteristic timescales of `x` are comparable,
or larger than those of `r`.

The second argument `μ` is a regularization parameter. In short, we estimate `r`
by minimizing a cost with two components: that `x` is close to `s/r` and that `x`
is smooth. `μ` is the multiplier of the smoothness cost.

You can give a vector as `μ`. The process will be repeated for all `μ` and the
[`rmse`](@ref) between `s` and the estimated `x * r` will be computed. The
`x` that gives the least error will be returned finally.
Use the low level `SignalDecomposition.matrix_invert(s, r, μ::Real) → x, err`
to get the error values.
"""
struct ProductInversion{R, M} <: Decomposition
    r::R
    μ::M
    x::R
    xdummy::R
end
ProductInversion(r, μ) = ProductInversion(r, μ, copy(r), copy(r))

function decompose(t, s, method::ProductInversion{R, <: Real}) where R
    x, err = matrix_invert(s, method.r, method.μ)
    return x, r
end

function decompose(t, s, method::ProductInversion{R, <: AbstractArray}) where R
    isequispaced(t) || error("Input time axis must be equispaced for method ProductInversion.")
    leasterr = Inf
    bestx = copy(s)
    bestμ = method.μ[1]
    for μ in method.μ
        x, err = matrix_invert(s, method.r, μ, method.x, method.xdummy)
        if err < leasterr
            bestx .= x
            leasterr = err
            bestμ = μ
        end
    end
    return bestx, method.r
end

function matrix_invert(s, r, μ::Real, x = copy(r), xdummy = copy(r))
    L = length(s)
    diags = Dict(
        0 => [6μ + r[i]^2 for i in 1:L],
        1 => fill(-4μ, L-1),
        -1 => fill(-4μ, L-1),
        2 => fill(μ, L-2),
        -2 => fill(μ, L-2),
    )
    # take care of border elements
    diags[1][1] = -8μ
    diags[2][1:2] .= 2μ
    diags[-1][end] = -8μ
    diags[-2][end-1:end] .= 2μ
    B = BandedMatrix(diags...)
    @. xdummy = s * r
    ldiv!(x, B, xdummy)
    @. xdummy = x * r
    err = rmse(s, xdummy)
    return x, err
end
