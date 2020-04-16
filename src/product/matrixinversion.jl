using BandedMatrices
export ProductInversion

"""
    ProductInversion(r, μ) <: Decomposition
Decompose a timeseries `s` into a **product** `x * r`, given that you have a good
estimate of one of the two factors `r` (the "input") and you need `x` (the "multiplier")
but you can't do simply `x =  s ./ r` because `r` contains zeros.

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
end

function decompose(t, s, method::ProductInversion{R, <: Real}) where R
    x, err = matrix_invert(s, method.r, method.μ)
    return x, r
end

function decompose(t, s, method::ProductInversion{R, <: AbstractArray}) where R
    leasterr = Inf
    bestx = copy(s)
    for μ in method.μ
        x, err = matrix_invert(s, method.r, method.μ)
        if err < leasterr
            bestx .= x
            leasterr = err
        end
    end
    return bestx, r
end

function matrix_invert(s, r, μ::Real)
    L = length(s)
    diags = Dict(i => zeros(eltype(s), L - abs(i)) for i in -2:2)
    for i in 1:L-2
        diags[0][i+2] = 6μ + s[i]^2
        diags[1][i+1] = -4μ
        diags[2][i] = μ
        diags[-1][i+1] = -4μ
        diags[-2][i] = μ
    end
    # take care of border elements
    for i in (L-1, L); diags[0][i] = 6μ + s[i]^2; end
    diags[1][1] = -8μ
    diags[1][end] = -4μ
    diags[2][1:2] .= 2μ
    diags[-1][1] = -4μ
    diags[-1][end] = -8μ
    diags[-2][end-1:end] .= 2μ
    B = BandedMatrix(diags...)
    # TODO: we can optimize this by not generating a new `x` and `x .* r` every time
    x = B \ (s .* r)
    err = rmse(s, x .* r)
    return x, err
end
