# Linear Programming Solver

Solves LP problems via the simplex method.

Ideally I'll implement the naive, revised, and full tableau implementations.
For now I've only written the naive (in [naive_simplex.jl](./naive_simplex.jl)) 

Each file will have plenty of comments to understand the approach - they pretty closely follow [the textbook](http://athenasc.com/linoptbook.html)'s way of thinking.
[The wikipedia article](https://en.wikipedia.org/wiki/Simplex_algorithm) has the same idea, but uses rather different notation and approach (at least for the naive implementation, the tableau is pretty similar)

# Example

(For `naive_simplex.jl`)

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
--- Iteration 1 ---
Beginning at vertex [0, 0, 0, 20, 20, 20] with cost 0
Found basic indicies: Int8[4, 5, 6]
Found negative reduced cost c̄_1 = -10.0
Index 1 will enter the basis
Index 5 will exit the basis, with θ* = 10.0
Moved to new vertex [10.0; 0.0; 0.0; 10.0; 0.0; 0.0;;]

--- Iteration 2 ---
Beginning at vertex [10.0; 0.0; 0.0; 10.0; 0.0; 0.0;;] with cost -100.0
Found basic indicies: Int8[4, 1, 6]
Found negative reduced cost c̄_2 = -7.0
Index 2 will enter the basis
Index 6 will exit the basis, with θ* = 0.0
Moved to new vertex [10.0; 0.0; 0.0; 10.0; 0.0; 0.0;;]

--- Iteration 3 ---
Beginning at vertex [10.0; 0.0; 0.0; 10.0; 0.0; 0.0;;] with cost -100.0
Found basic indicies: Int8[4, 1, 2]
Found negative reduced cost c̄_3 = -9.0
Index 3 will enter the basis
Index 4 will exit the basis, with θ* = 4.0
Moved to new vertex [4.0; 4.0; 4.0; 0.0; 0.0; 0.0;;]

--- Iteration 4 ---
Beginning at vertex [4.0; 4.0; 4.0; 0.0; 0.0; 0.0;;] with cost -136.0
Found basic indicies: Int8[3, 1, 2]

Found optimal solution: [4.0; 4.0; 4.0; 0.0; 0.0; 0.0;;]
```