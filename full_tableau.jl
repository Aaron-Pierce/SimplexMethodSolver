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

    # We'll construct the initial simplex tableau here,
    # and will update it at the end of each iteration.
    # The simpled tableau consists of B^-1 A, which is 
    # useful to be able to pick off columns for anything
    # that needs u or d_b. Reduced cost, minimal l, pretty
    # much everything needs B^-1 A_something.
    # There's also an extra column, the first column,
    # which is B^-1 b, the values of your basic variables.
    # And an extra row, the reduced costs of each variable.
    # Finally, the top left cell, (1, 1), the intersection of the new row and col
    # keeps -c_b ⋅ b, which represents the current negative cost. 

    tableau = Matrix{Float64}(undef, num_constraints + 1, num_variables + 1)


    # These are annotated with a double underscore
    # to denote that they won't be used in the main
    # loop, only while constructing the initial tableau.
    __basic_matrix = Matrix{Float64}(undef, num_constraints, num_constraints)
    for i in 1:length(basic_indicies)
        __basic_matrix[:, i] = A[:, basic_indicies[i]]
    end
    __b_inv = inv(__basic_matrix)

    # We'll first populate the B^-1 A part.
    for i in 1:num_variables
        tableau[:, i+1] = vcat(0, __b_inv * A[:, i])
    end

    # Then the first column, B^-1 b
    tableau[:, 1] = vcat(0, __b_inv * b)

    
    # Then the first row, c_j - c_b ⋅ B^-1 A_j for all j
    for j in 1:num_variables
        tableau[1, 1 + j] = c[j] - dot(c[basic_indicies], __b_inv * A[:, j])
    end

    # And finally, the negative overall cost, c_b ⋅ B^-1b
    tableau[1, 1] = dot(c[basic_indicies], __b_inv * b)

    # Phew, bunch of work, but it'll be easy to update later!

    println("Initial tableau:")
    display(tableau)
    println()

    

    # 300 is an arbitrary number of iterations to terminate
    # at, in case it cycles (which shouldn't happen because
    # I always select the smallest index j that has a negative 
    # reduced cost, but hey anything can happen, who knows)
    for iter in 1:4

        println("")
        println("--- Iteration ", iter, " ---")
        println("Beginning at vertex ", basicSolution, " with cost ", dot(c, basicSolution))
        println("and basic indicies: ", basic_indicies)

        #=
            Step 1, figure out what direction to move in
        =#


        # We want to move in a cost-reducing manner.
        # Our tableau very kindly keeps track of 
        # All of the reduced costs already.
        # All we need to do is find a negative one!
        index_entering_basis = -1 # this is j
        direction = Vector{Float64}(undef, num_variables) # this is d⃗

        for j in 1:length(c)
            # We only care about non-basic variables
            if j in basic_indicies
                continue
            end

            # We'll just pull the reduced cost straight
            # off the table, selecting the first
            # one that is negative, as an anti-
            # cycling strategy.

            reduced_cost = tableau[1, 1 + j]

            if reduced_cost < 0
                println("Found negative reduced cost c̄_", j, " = ", reduced_cost)
                index_entering_basis = j

                # Now we need to compute the direction
                # vector we move along to take advantage
                # of this negative reduced cost.
                # Thankfully, we can easily find -d_b by picking off
                # B^-1 A_j, which is the (j+1)'st tableau col

                d = zeros(num_variables)

                d[basic_indicies] = -tableau[2:(num_constraints+1), j+1]
                d[j] = 1

                println("created direction ", d)

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
            println("Final tableau:")
            display(tableau)
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
        # The tableau will really help us out.
        # Each column holds the value of -d_b,
        # we just need to find the smallest magnitude
        # positive element of -d_b

        θstar = Inf
        l = -1

        u = tableau[2:(num_constraints+1), 1 + index_entering_basis]

        for i in 1:size(u, 1) # size(u) == num_constraints
            if u[i] > 0
                newθ = basicSolution[basic_indicies[i]] / u[i]
                if newθ < θstar
                   θstar = newθ
                   l = i
                end
            end
        end


        println("Index B(", l, ") = ", basic_indicies[l], " will exit the basis, with θ* = ", θstar)

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
        # when we wouldn't be able to just read off the zero entries
        basic_indicies[l] = index_entering_basis

        # And finally, we do the row reductions.
        # B^-1 B̄ = I, except for the lth column, which is valued B^-1 A_j
        # Your goal is to figure out how to make Q B^-1 B̄ = I, for real, no exceptions.
        # How do you do that? Row reduce B^-1A_j into the l'th identity column.
        # Thankfully, we already know what B^-1A_j is from the tableau, so this is easy.

        # All of the columns in tableau, even the extra ones, rely on B^-1 in some way,
        # (the reduced cost row is a little more complicated) But the value of that table
        # is B^-1 times something. B̄^-1 = QB^-1. Matrix multiplication associates,
        # so hitting the table with Q is the same as hitting B̄^-1 by the right hand
        # side of whatever is necessary to build that side of the table, so applying
        # the row reductions directly to the tableau still works out.

        # This is (m+1)x(m+1) to match the dimension of the 
        # rows of the tableau. 
        Q = 1 * Matrix{Float64}(I, num_constraints+1, num_constraints+1)

        # Knock out all rows using the l'th
        for i in 1:(num_constraints+1)
            if i != l+1 # l + 1, not l, because of the extra row
                Q[i,:] -= Q[l+1,:] * (tableau[i, index_entering_basis+1] / u[l]) 
            end
        end

        Q[l+1, :] /= u[l]

        tableau = Q * tableau

        println("Built new tableau: ")
        display(tableau)
        println()
    end
end

solve(
    [-10, -12, -12, 0, 0, 0], # Cost vector (c)
    [1 2 2 1 0 0 # Constraint matrix (A)
        2 1 2 0 1 0
        2 2 1 0 0 1],
    [
        20 # b, as in Ax = b
        20
        20
    ],
    [0, 0, 0, 20, 20, 20] # An initial basic feasible solution
)