using CSV, DataFrames, Luxor, Colors, Cairo

include(joinpath(dirname(dirname(@__FILE__)), "src", "messages.jl"))
include(joinpath(dirname(dirname(@__FILE__)), "src", "distance.jl"))
include(joinpath(dirname(dirname(@__FILE__)), "src", "utils.jl"))
include(joinpath(dirname(dirname(@__FILE__)), "src", "primes.jl"))

function integer_search(stop_at::Integer)::Array{Array{Number, 1}}
    A = []
    upper_bound = 10
    increment_bound = 10
    processed = []
    i = 0

    while true
        for q in 1:upper_bound, n in 1:upper_bound, d in 1:upper_bound
            # skip configurations that have already been processed
            # arelessthan(upper_bound - increment_bound, q, n, d) && continue
            (q, n, d) ∈ processed && continue
            hb = hamming_bound(q, n, d, no_round)
            push!(processed, (q, n, d))
            # i += 1; println("$i:\t$q, $n, $d")
            push!(A, [q, n, d, hb])
            isequal(length(A), stop_at) && return A
        end
        upper_bound += increment_bound
    end
end

function make_integer_df(stop_at::Integer)
    A = integer_search(stop_at)
    D = DataFrame(; q = Number[], n = Number[], d = Number[], hamming_bound = Number[])

    for i in A
        push!(D, i)
    end

    return D
end

function draw_turtle(stop_at::Integer; specific::Bool = false)
    rigidness = "non_rigid"
    if specific
        rigidness = "rigid"
    end

    df = make_integer_df(stop_at)

    golden_angle = 137.5
    angle = golden_angle
    # angle = 42

    # Drawing(600, 400, joinpath(dirname(dirname(@__FILE__)), "other", "turtle_drawing.pdf"))
    Drawing(
        10000,
        10000,
        joinpath(
            homedir(),
            "projects",
            "CodingTheory.jl",
            "other",
            "$(rigidness)",
            "turtle_drawing_$(stop_at)_$(angle).pdf",
        ),
    )
    origin()
    # background("midnightblue")

    🐢 = Turtle()
    Pencolor(🐢, "black")
    Penwidth(🐢, 1.5)
    n = 5

    # distance shouldn't be larger than the block length; filter trivial distances of one; we filter out combinations when all are equal; filter out codes that are generalised hamming codes; filter out even distances
    for row in eachrow(df)
        if specific
            if row.d < row.n &&
                !isone(row.d) &&
                !allequal(row.q, row.n, row.d) &&
                !iszero(mod(row.d, 2)) &&
                !isone(row.hamming_bound) &&
                !ishammingbound(BigInt(row.q), BigInt(row.n), BigInt(row.d))
                continue
            end
        else
            if row.d < row.n && !isone(row.d)
                continue
            end
        end

        if isinteger(row.hamming_bound)
            Turn(🐢, angle)
            # HueShift(🐢)
        end
        Forward(🐢, n)
        # n += 0.75
    end

    # fontsize(20)
    # Message(🐢, "finished")
    Luxor.finish()

    return nothing
end
