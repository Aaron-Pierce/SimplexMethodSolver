# Linear Programming Solver

Solves LP problems via the simplex method.

~~Ideally I'll implement the naive, revised, and full tableau implementations.~~ 
~~For now~~ I've ~~only~~ written the naive (in [naive_simplex.jl](./naive_simplex.jl)),
revised (in [revised_simplex.jl](./revised_simplex.jl)), and full tableau methods ([full_tableau.jl](./full_tableau.jl)).

Each file will have plenty of comments to understand the approach - they pretty closely follow [the textbook](http://athenasc.com/linoptbook.html)'s way of thinking.
[The wikipedia article](https://en.wikipedia.org/wiki/Simplex_algorithm) has the same idea, but uses a rather different notation and approach (at least for the naive implementation, the tableau is pretty similar but uses a different notation to describe the tableau's contents)

`naive_simplex.jl`, `revised_simplex.jl`, and `full_tableau.jl` all require the user to
find a non-degenerate initial basic feasible solution to seed the algorithm with.
`full_tableau.jl` will let you supply a degenerate BFS so long as you give it 
the basic indicies, but even finding one BFS by hand can be a pain.
`finds_initial_solution.jl` runs the full tableau implementation, but
will do the auxillary problem of finding an initial solution for you.
This is very nice because it can identify infeasible problems
and redundant constraints right at the outset, so this is the
one to use if you actually need to solve an LP problem.

# Example

(From `full_tableau.jl`)
```julia

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

```

Which outputs:

```
Initial tableau:
4×7 Matrix{Float64}:
  0.0  -10.0  -12.0  -12.0  0.0  0.0  0.0
 20.0    1.0    2.0    2.0  1.0  0.0  0.0
 20.0    2.0    1.0    2.0  0.0  1.0  0.0
 20.0    2.0    2.0    1.0  0.0  0.0  1.0

--- Iteration 1 ---
Beginning at vertex [0, 0, 0, 20, 20, 20] with cost 0
and basic indicies: Int8[4, 5, 6]
Found negative reduced cost c̄_1 = -10.0
created direction [1.0, 0.0, 0.0, -1.0, -2.0, -2.0]
Index 1 will enter the basis
Index B(2) = 5 will exit the basis, with θ* = 10.0
Moved to new vertex [10.0, 0.0, 0.0, 10.0, 0.0, 0.0]
Built new tableau: 
4×7 Matrix{Float64}:
 100.0  0.0  -7.0  -2.0  0.0   5.0  0.0
  10.0  0.0   1.5   1.0  1.0  -0.5  0.0
  10.0  1.0   0.5   1.0  0.0   0.5  0.0
   0.0  0.0   1.0  -1.0  0.0  -1.0  1.0

--- Iteration 2 ---
Beginning at vertex [10.0, 0.0, 0.0, 10.0, 0.0, 0.0] with cost -100.0
and basic indicies: Int8[4, 1, 6]
Found negative reduced cost c̄_2 = -7.0
created direction [-0.5, 1.0, 0.0, -1.5, 0.0, -1.0]
Index 2 will enter the basis
Index B(3) = 6 will exit the basis, with θ* = 0.0
Moved to new vertex [10.0, 0.0, 0.0, 10.0, 0.0, 0.0]
Built new tableau: 
4×7 Matrix{Float64}:
 100.0  0.0  0.0  -9.0  0.0  -2.0   7.0
  10.0  0.0  0.0   2.5  1.0   1.0  -1.5
  10.0  1.0  0.0   1.5  0.0   1.0  -0.5
   0.0  0.0  1.0  -1.0  0.0  -1.0   1.0

--- Iteration 3 ---
Beginning at vertex [10.0, 0.0, 0.0, 10.0, 0.0, 0.0] with cost -100.0
and basic indicies: Int8[4, 1, 2]
Found negative reduced cost c̄_3 = -9.0
created direction [-1.5, 1.0, 1.0, -2.5, 0.0, 0.0]
Index 3 will enter the basis
Index B(1) = 4 will exit the basis, with θ* = 4.0
Moved to new vertex [4.0, 4.0, 4.0, 0.0, 0.0, 0.0]
Built new tableau: 
4×7 Matrix{Float64}:
 136.0  0.0  0.0  2.22045e-16   3.6   1.6   1.6
   4.0  0.0  0.0  1.0           0.4   0.4  -0.6
   4.0  1.0  0.0  0.0          -0.6   0.4   0.4
   4.0  0.0  1.0  0.0           0.4  -0.6   0.4

--- Iteration 4 ---
Beginning at vertex [4.0, 4.0, 4.0, 0.0, 0.0, 0.0] with cost -136.0
and basic indicies: Int8[3, 1, 2]

Found optimal solution: [4.0, 4.0, 4.0, 0.0, 0.0, 0.0]

Final tableau:
4×7 Matrix{Float64}:
 136.0  0.0  0.0  2.22045e-16   3.6   1.6   1.6
   4.0  0.0  0.0  1.0           0.4   0.4  -0.6
   4.0  1.0  0.0  0.0          -0.6   0.4   0.4
   4.0  0.0  1.0  0.0           0.4  -0.6   0.4

```