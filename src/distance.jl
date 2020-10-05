  
#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $0))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

function hamming_distance(s1::AbstractString, s2::AbstractString)::Integer
    if ! isequal(length(s1), length(s2))
        throw(error("Cannot compute Hamming Distance on strings of unequal length."))
    end
    
    distance = 0
        
    for i in 1:length(s1)
        if s1[i] ≠ s2[i]
            distance += 1
        end
    end
    
    return distance
end

function code_distance(C::Array{String,1})::Integer
	distances = Integer[]
	
	for c in C, c′ in C
		if c ≠ c′
			push!(distances, hamming_distance(c, c′))
		end
	end

	return minimum(distances)
end

function code_distance(C::Array{Array{Int64, 1}, 1})::Integer
	string_code = String[]
	
	for c in C
		push!(string_code, join(c))
	end
	
	return code_distance(string_code)
end
