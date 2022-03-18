using LinearAlgebra

function find_initial_solution(problem_cost::Vector, problem_A::Matrix, problem_b::Vector)
    # We'll solve another LP problem to find a solution
    # To construct that problem, we add m extra decision variables, one per constraint,
    # Which just have the values of b, which gives an easy solution
    # to this newly constructed problem. We format our cost function
    # as the sum of all of these variables, so that an optimal solution
    # occurs when they're all zero, so the basic variables that make
    # up this solution will also be a basic feasible solution for our
    # original problem, and then we can proceed to simplex on it

    num_constraints = size(problem_A, 1)
    num_variables = size(problem_A, 2) + num_constraints # We add one new variable per constraint

    basic_indicies = collect(size(problem_A, 2)+1:num_variables) # All of the auxiliary variables

    c = zeros(num_variables)
    c[basic_indicies] .= 1

    println("created cost function: ")
    display(c)
    println()

    A = copy(problem_A)
    b = copy(problem_b)

    for i in 1:num_constraints
        if b[i] < 0
            b[i] *= -1
            A[i, :] .*= -1
        end
    end

    A = hcat(A, Matrix{Float64}(I, num_constraints, num_constraints))

    # Construct the tableau
    tableau = Matrix{Float64}(undef, num_constraints + 1, num_variables + 1)

    # The cost in the top left, should be the (negative) sum of b
    tableau[1, 1] = -dot(c[basic_indicies], b)

    # The values of basic variables, again just b
    tableau[2:(num_constraints+1)] = b

    # The reduced costs, c_j - B^-1A_j
    for j in 1:num_variables
        tableau[1, j+1] = c[j] - dot(c[basic_indicies], A[:, j]) # B is the identity, we can skip it
    end

    # The entire A matrix, because B^-1 is the identity
    tableau[2:num_constraints+1, 2:num_variables+1] = A

    println("Constructed initial auxiliary tableau:")
    display(tableau)
    println()

    println("Now just run the simplex method on it")

    currentSolution = zeros(num_variables)
    currentSolution[size(problem_A, 2)+1:num_variables] = b


    for iter in 1:300

        if tableau[1, 1] == 0
            println("Terminated at cost 0")


            redundant_constraints = []
            problem_basic_indicies = copy(basic_indicies)
            for ind in 1:size(basic_indicies, 1)
                if basic_indicies[ind] > num_constraints
                    # An artificial variable is in the basis, we need to swap it for a 'real' variable
                    # We need to search for a non-zero pivot element to do this

                    index_entering = -1
                    print("Checking: ")
                    for j in 1:size(problem_A, 2) # number of 'real' variables
                        if tableau[ind+1, j+1] != 0
                            index_entering = j
                            break
                        end
                    end

                    if index_entering != -1
                        println("We can swap artificial " * string(basic_indicies[ind]) * " for real " * string(index_entering))
                        Q = 1 * Matrix{Float64}(I, num_constraints + 1, num_constraints + 1)

                        # Knock out all rows using the l'th
                        for i in 1:(num_constraints+1)
                            if i != ind + 1 # l + 1, not l, because of the extra row
                                Q[i, :] -= Q[ind+1, :] * (tableau[i, index_entering+1] / tableau[ind+1, index_entering+1])
                            end
                        end

                        Q[ind+1, :] /= tableau[ind+1, index_entering+1]

                        tableau = Q * tableau

                        for k in 1:size(problem_basic_indicies)
                            if problem_basic_indicies[k] == ind
                                problem_basic_indicies[k] = index_entering
                            end
                        end
                    else
                        println("Found redundant constraint: " * string(ind))
                        append!(redundant_constraints, ind)
                    end

                end
            end


            problem_tableau = tableau[setdiff(collect(1:num_constraints+1), redundant_constraints .+ 1), 1:size(problem_A, 2)+1]
            problem_basic_indicies = deleteat!(problem_basic_indicies, redundant_constraints)

            # we need to update reduced costs for the real problem now 

            for j in 1:size(problem_A, 2)
                problem_tableau[1, j+1] = problem_cost[j] - dot(problem_cost[problem_basic_indicies], problem_tableau[2:end, j+1])
            end

            println("Final tableau: ")
            display(problem_tableau)
            println()
            return problem_tableau, problem_basic_indicies
        end

        println("")
        println("--- Iteration ", iter, " ---")
        println("Beginning at vertex ", currentSolution, " with cost ", dot(c, currentSolution))
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

        for j in 1:num_variables
            # We only care about non-basic variables
            if j in basic_indicies
                continue
            end

            # We'll just pull the reduced cost straight
            # off the table, selecting the first
            # one that is negative, as an anti-
            # cycling strategy.

            reduced_cost = tableau[1, 1+j]

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

            if tableau[1, 1] != 0
                println("Auxiliary problem terminated with cost " * tableau[1, 1])
                if tableau[1, 1] > 0
                    println("The initial problem was infeasible")
                else
                    println("Optimal cost was not positive, why wasn't this caught earlier in the loop?")
                    return
                end
            end

            println("")
            printstyled("Found optimal solution: " * string(currentSolution) * "\n", color = :green)
            println("")
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

        u = tableau[2:(num_constraints+1), 1+index_entering_basis]

        for i in 1:size(u, 1) # size(u) == num_constraints
            if u[i] > 0
                newθ = currentSolution[basic_indicies[i]] / u[i]
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
        currentSolution += θstar .* direction
        println("Moved to new vertex ", currentSolution)

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
        Q = 1 * Matrix{Float64}(I, num_constraints + 1, num_constraints + 1)

        # Knock out all rows using the l'th
        for i in 1:(num_constraints+1)
            if i != l + 1 # l + 1, not l, because of the extra row
                Q[i, :] -= Q[l+1, :] * (tableau[i, index_entering_basis+1] / u[l])
            end
        end

        Q[l+1, :] /= u[l]

        tableau = Q * tableau

        println("Built new tableau: ")
        display(tableau)
        println()
    end

end

function solve(c::Vector, A::Matrix, b::Vector)
    tableau, basic_indicies = find_initial_solution(c, A, b)


    println("\n\n-------------------------")
    println("Found an initial BFS with basic indicies: " * string(basic_indicies) * " and initial tableau ")
    display(tableau)
    println()

    num_variables = size(tableau, 2) - 1
    num_constraints = size(tableau, 1) - 1
    current_solution = zeros(num_variables)
    current_solution[basic_indicies] = tableau[basic_indicies .+ 1, 1]

    for iter in 1:300

        println("")
        println("--- Iteration ", iter, " ---")
        println("Beginning at vertex ", current_solution, " with cost ", dot(c, current_solution))
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

            reduced_cost = tableau[1, 1+j]

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
            printstyled("Found optimal solution: " * string(current_solution) * "\n", color = :green)
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
        # current_solution[i] + θ⋅direction[i] = 0
        # The tableau will really help us out.
        # Each column holds the value of -d_b,
        # we just need to find the smallest magnitude
        # positive element of -d_b

        θstar = Inf
        l = -1

        u = tableau[2:(num_constraints+1), 1+index_entering_basis]

        for i in 1:size(u, 1) # size(u) == num_constraints
            if u[i] > 0
                newθ = current_solution[basic_indicies[i]] / u[i]
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
        current_solution += θstar .* direction
        println("Moved to new vertex ", current_solution)

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
        Q = 1 * Matrix{Float64}(I, num_constraints + 1, num_constraints + 1)

        # Knock out all rows using the l'th
        for i in 1:(num_constraints+1)
            if i != l + 1 # l + 1, not l, because of the extra row
                Q[i, :] -= Q[l+1, :] * (tableau[i, index_entering_basis+1] / u[l])
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
    [0, 0, 0, -2, -3, 1, 12],   # Cost vector (c)
    [1 0 0 -2 -9 1 9 # Constraint matrix (A)
        0 1 0 1/3 1 -1/3 -2
        0 1 0 1/3 1 -1/3 -2
        0 1 0 1/3 1 -1/3 -2
        0 0 1 2 3 -1 -12],
    [
        0 # b, as in Ax = b
        0
        0
        0
        2
    ],
)