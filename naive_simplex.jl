using LinearAlgebra

function solve(c::Vector, A::Matrix, b::Vector, basicSolution::Vector)

    # useful to not call size(A) over and over again
    num_variables = size(A, 2)
    num_constraints = size(A, 1)

    # The user should pass in a non-degenerate
    # BFS to start with, so we can figure out the 
    # basic indicies on our own. Later we'll
    # have to manually update these in case we
    # arrive at a degenerate solution, where 
    # some basic variables will be zero-valued
    basic_indicies = Vector{Int8}(undef, 0)
    for i in 1:length(basicSolution)
        if basicSolution[i] != 0.0
            append!(basic_indicies, i)
        end
    end

    # 300 is an arbitrary number of iterations to terminate
    # at, in case it cycles (which shouldn't happen because
    # I always select the smallest index j that has a negative 
    # reduced cost, but hey anything can happen, who knows)
    for iter in 1:300

        println("")
        println("--- Iteration ", iter, " ---")
        println("Beginning at vertex ", basicSolution, " with cost ", dot(c, basicSolution))



        println("Found basic indicies: ", basic_indicies)

        basic_matrix = Matrix{Float64}(undef, num_constraints, num_constraints)
        for i in 1:length(basic_indicies)
            basic_matrix[:, i] = A[:, basic_indicies[i]]
        end

        #=
            Step 1, figure out what direction to move in
        =#

        # We want to move in a cost-reducing manner, let's compute reduced costs!   
        # we're looking for a direction vector d⃗ to move along. dⱼ will equal 1, the index of 
        # the variable that will enter the basis, and by allowing it to enter the basis,
        # the remaining basic variables must change to allow for that entry. Their values
        # must conform to Ad = 0 ⟹ Bdᵦ + Aⱼ, so dᵦ = -B⁻¹Aⱼ

        # Let's compute all the reduced costs!
        index_entering_basis = -1
        direction = Vector{Float64}(undef, num_variables)

        for j in 1:length(c)
            # We only care about non-basic variables
            if j in basic_indicies
                continue
            end

            # Initialize the direction vector
            # The non-basics will be 0, except for j,
            # and we'll need to calculate dᵦ
            d = zeros(Float64, num_variables, 1)
            d[j] = 1


            basic_directions = (inv(basic_matrix) * (-A[:, j]))
            for i in 1:length(basic_indicies)
                d[basic_indicies[i]] = basic_directions[i]
            end

            # Okay, now we have d, time to compute the reduced cost
            reduced_cost = dot(c, d)
            
            # We're not doing anything fancy, first cost-minimizing direction
            # is good enough for us! This will hopefully prevent cycling, too
            if reduced_cost < 0
                println("Found negative reduced cost c̄_", j, " = ", reduced_cost)
                index_entering_basis = j
                direction = d
                break
            end
        end

        if index_entering_basis == -1
            # no reduced cost was negative,
            # so we are unable to find a lower-cost
            # solution, so we're optimal!
            println("")
            printstyled("Found optimal solution: " * string(basicSolution) * "\n", color = :green)
            println("")
            return
        end

        println("Index ", index_entering_basis, " will enter the basis")

        #=
            Step 2, figure out how far we can move in that direction
        =#

        # We can move in the direction `direction` until
        # a basic variable hits zero value. I.e
        # basicSolution[i] + θ⋅direction[i] = 0

        θstar = Inf
        l = -1

        for i in 1:num_variables
            if !(i in basic_indicies) || direction[i] >= 0
                continue
            end
            newθ = -basicSolution[i] / direction[i]
            if newθ < θstar
                θstar = newθ
                l = i
            end
        end


        println("Index ", l, " will exit the basis, with θ* = ", θstar)

        if θstar == Inf
            println("Could move forever, optimal cost is -∞")
            return
        end

        #= 
            Step 3, move and repeat!
        =#
        basicSolution += θstar .* direction
        println("Moved to new vertex ", basicSolution)

        # We also need to manually update the basic indicies, in case of degeneracy
        for i in 1:length(basic_indicies)
            if basic_indicies[i] == l
                basic_indicies[i] = index_entering_basis
                break
            end
        end
    end
end

solve(
    [-10, -12, -12, 0, 0, 0], # Cost vector (c)
    [   1 2 2 1 0 0 # Constraint matrix (A)
        2 1 2 0 1 0
        2 2 1 0 0 1 ],
    [
        20 # b, as in Ax = b
        20
        20
    ],
    [0, 0, 0, 20, 20, 20] # An initial basic feasible solution
)