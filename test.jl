include("src/messages.jl")

function hamming_sphere(M::AbstractArray{T}, w::AbstractArray, e::Integer) where T
	if e < 0
	    throw(error("e (the ball \"radius\") must be a non-negative number."))
	end
	
	if eltype(w) isa Symbol # ensure element types of the input are interpreted as symbols
	else
		w = __deepsym(w)
	end

	sphere = [] # initialise the ball with no words
	for v in __deepsym(M)
		if hamming_distance(w, v) == e # if the words in the set of messages are within radius e of the word w, add them to the list of words
		    push!(sphere, v)
		end
	end
	
	return sphere # return said list of words
end

println(get_codewords(["a", "b", "c"], 3, 3, 1))
println()
println(hamming_sphere(get_codewords(["a", "b", "c"], 3, 3, 1), ["a", "b", "c"], 2))
