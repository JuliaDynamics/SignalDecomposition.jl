function isequispaced(t)
    d1 = t[2] - t[1]
    for i in 3:length(t)
        d = t[i] - t[i-1]
        d == d1 || return false
    end
    return true
end

"""
    findnearest(val, A)
Return the index of `A` which has value nearest to `val`.
"""
function findnearest(val::Real, A::AbstractVector{<:Real})
    i = 1
    d = abs(val - A[i])
    @inbounds for j in 1:length(A)
        dd = abs(val - A[j])
        if dd < d
            i = j
            d = dd
        end
    end
    return i
end
