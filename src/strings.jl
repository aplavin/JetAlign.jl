function highlight_string_differences(strs; colors)
    length(strs) == 2 || error("Currently supports exactly 2 strings")
    
    str1, str2 = strs
    color1, color2 = colors
    
    # Simple alignment using edit distance approach
    alignment = align_strings(str1, str2)
    
    result1_parts = []
    result2_parts = []
    
    for (op, char1, char2) in alignment
        if op == :match
            push!(result1_parts, rich(string(char1)))
            push!(result2_parts, rich(string(char2)))
        elseif op == :substitute  
            push!(result1_parts, rich(string(char1); color=color1))
            push!(result2_parts, rich(string(char2); color=color2))
        elseif op == :delete
            push!(result1_parts, rich(string(char1); color=color1))
        elseif op == :insert
            push!(result2_parts, rich(string(char2); color=color2))
        end
    end
    
    [rich(result1_parts...), rich(result2_parts...)]
end

function align_strings(s1, s2)
    chars1, chars2 = collect(s1), collect(s2)
    m, n = length(chars1), length(chars2)
    
    # Dynamic programming table for edit distance with operations tracking
    dp = Matrix{Int}(undef, m + 1, n + 1)
    ops = Matrix{Symbol}(undef, m + 1, n + 1)
    
    # Initialize base cases
    @p 0:m map((dp[_ + 1, 1] = _; ops[_ + 1, 1] = :delete))
    @p 0:n map((dp[1, _ + 1] = _; ops[1, _ + 1] = :insert))
    
    # Fill the DP table
    for i in 1:m, j in 1:n
        if chars1[i] == chars2[j]
            dp[i + 1, j + 1] = dp[i, j]
            ops[i + 1, j + 1] = :match
        else
            costs = (
                dp[i, j] + 1,     # substitute
                dp[i, j + 1] + 1, # delete  
                dp[i + 1, j] + 1  # insert
            )
            min_cost, min_idx = findmin(costs)
            dp[i + 1, j + 1] = min_cost
            ops[i + 1, j + 1] = (:substitute, :delete, :insert)[min_idx]
        end
    end
    
    # Backtrack to get alignment
    alignment = []
    i, j = m + 1, n + 1
    
    while i > 1 || j > 1
        op = ops[i, j]
        elt = if op == :match
            i, j = i - 1, j - 1
            (:match, chars1[i], chars2[j])
        elseif op == :substitute
            i, j = i - 1, j - 1
            (:substitute, chars1[i], chars2[j])
        elseif op == :delete
            i = i - 1
            (:delete, chars1[i], nothing)
        elseif op == :insert
            j = j - 1
            (:insert, nothing, chars2[j])
        end
        pushfirst!(alignment, elt)
    end
    
    alignment
end
