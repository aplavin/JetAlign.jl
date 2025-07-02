round_(x::Number, precision::Number) = round(Integer, x / precision) * precision

function highlight_string_differences(strs; colors)
	if allequal(length, strs)
		equalchars = @p let
			strs
			map() do str
				collect(str) .== collect(first(strs))
			end
			stack
			eachslice(dims=1)
			map(all)
		end
		map(strs, colors) do str, color
			parts = map(collect(str), equalchars) do ch, isequal
				kws = isequal ? (;) : (;color)
				rich(string(ch); kws...)
			end
			rich(parts...)
		end
	else
		neqchars_start = @p let
			0:minimum(length, strs)
			filter() do i
				allequal(s -> (@view s[begin:i]), strs)
			end
			maximum
		end
		neqchars_end = @p let
			0:minimum(length, strs)
			filter() do i
				allequal(s -> (@view s[end-i+1:end]), strs)
			end
			minimum
		end
		map(strs, colors) do str, color
			rich(
				str[1:neqchars_start],
				rich(str[neqchars_start+1:end-neqchars_end]; color),
				str[end-neqchars_end+1:end],
			)
		end
	end
end
