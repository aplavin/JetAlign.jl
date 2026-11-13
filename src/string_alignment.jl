"""
    difference_runs(s1, s2) -> (runs1, runs2)

Split each string into consecutive runs `(; text, differs)`, where `differs` marks characters
that have no match in the other string (per [`align_strings`](@ref)).
"""
function difference_runs(s1, s2)
    alignment = align_strings(s1, s2)
    map((2, 3)) do k
        chars = filtermap(a -> isnothing(a[k]) ? nothing : (a[k], a[1] !== :match), alignment)
        starts = filter(i -> i == firstindex(chars) || chars[i][2] != chars[i-1][2], eachindex(chars))
        stops = filter(i -> i == lastindex(chars) || chars[i][2] != chars[i+1][2], eachindex(chars))
        map(starts, stops) do a, b
            (; text=String(first.(chars[a:b])), differs=chars[a][2])
        end
    end
end

"""
    align_strings(s1, s2)

Minimal edit-distance alignment of two strings, as a vector of `(op, char1, char2)` with `op` one of
`:match`, `:substitute`, `:delete` (char only in `s1`), `:insert` (char only in `s2`).
"""
function align_strings(s1, s2)
    chars1, chars2 = collect(s1), collect(s2)
    m, n = length(chars1), length(chars2)

    # edit distances and the last operation of each optimal prefix alignment
    dp = Matrix{Int}(undef, m + 1, n + 1)
    ops = Matrix{Symbol}(undef, m + 1, n + 1)
    dp[:, 1] .= 0:m; ops[:, 1] .= :delete
    dp[1, :] .= 0:n; ops[1, :] .= :insert

    for i in 1:m, j in 1:n
        if chars1[i] == chars2[j]
            dp[i + 1, j + 1] = dp[i, j]
            ops[i + 1, j + 1] = :match
        else
            min_cost, min_idx = findmin((dp[i, j], dp[i, j + 1], dp[i + 1, j]))
            dp[i + 1, j + 1] = min_cost + 1
            ops[i + 1, j + 1] = (:substitute, :delete, :insert)[min_idx]
        end
    end

    alignment = Tuple{Symbol, Union{Char,Nothing}, Union{Char,Nothing}}[]
    i, j = m + 1, n + 1
    while i > 1 || j > 1
        op = ops[i, j]
        elt = if op === :match || op === :substitute
            i, j = i - 1, j - 1
            (op, chars1[i], chars2[j])
        elseif op === :delete
            i -= 1
            (op, chars1[i], nothing)
        else
            j -= 1
            (op, nothing, chars2[j])
        end
        pushfirst!(alignment, elt)
    end
    alignment
end
