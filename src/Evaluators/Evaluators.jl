
"""
    AbstractNLPEvaluator

AbstractNLPEvaluator implements the bridge between the
problem formulation (see `ExaPF.AbstractFormulation`) and the optimization
solver. Once the problem formulation bridged, the evaluator allows
to evaluate:
- the objective;
- the gradient of the objective;
- the constraints;
- the Jacobian of the constraints;
- the Jacobian-vector and transpose-Jacobian vector products of the constraints;
- the Hessian of the objective;
- the Hessian of the Lagrangian.

"""
abstract type AbstractNLPEvaluator end

"""
    optimize!(optimizer, nlp::AbstractNLPEvaluator, x0)

Use optimization routine implemented in `optimizer` to optimize
the optimal power flow problem specified in the evaluator `nlp`.
Initial point is specified by `x0`.

Return the solution as a named tuple, with fields
- `status::MOI.TerminationStatus`: Solver's termination status, as specified by MOI
- `minimum::Float64`: final objective
- `minimizer::AbstractVector`: final solution vector, with same ordering as the `Variables` specified in `nlp`.


    optimize!(optimizer, nlp::AbstractNLPEvaluator)

Wrap previous `optimize!` function and pass as initial guess `x0`
the initial value specified when calling `initial(nlp)`.

## Examples

```julia
nlp = ExaPF.ReducedSpaceEvaluator(datafile)
optimizer = Ipopt.Optimizer()
solution = ExaPF.optimize!(optimizer, nlp)

```

## Notes
By default, the optimization routine solves a minimization
problem.

"""
function optimize! end

abstract type AbstractNLPAttribute end

"""
    Variables <: AbstractNLPAttribute end

Attribute corresponding to the optimization variables attached
to a given [`AbstractNLPEvaluator`](@ref).
"""
struct Variables <: AbstractNLPAttribute end

"""
    Constraints <: AbstractNLPAttribute end

Attribute corresponding to the constraints  attached
to a given [`AbstractNLPEvaluator`](@ref).
"""
struct Constraints <: AbstractNLPAttribute end

"""
    n_variables(nlp::AbstractNLPEvaluator)
Get the number of variables in the problem.
"""
function n_variables end

"""
    n_constraints(nlp::AbstractNLPEvaluator)
Get the number of constraints in the problem.
"""
function n_constraints end

# Callbacks
"""
    objective(nlp::AbstractNLPEvaluator, u)::Float64

Evaluate the objective at given variable `u`.
"""
function objective end

"""
    gradient!(nlp::AbstractNLPEvaluator, g, u)

Evaluate the gradient of the objective, at given variable `u`.
Store the result inplace in the vector `g`.

## Note
The vector `g` should have the same dimension as `u`.

"""
function gradient! end

"""
    constraint!(nlp::AbstractNLPEvaluator, cons, u)

Evaluate the constraints of the problem at given variable `u`. Store
the result inplace, in the vector `cons`.

## Note
The vector `cons` should have the same dimension as the result
returned by `n_constraints(nlp)`.

"""
function constraint! end

"""
    jacobian_structure(nlp::AbstractNLPEvaluator)

Return the sparsity pattern of the Jacobian matrix.
Return two vectors `rows` and `cols` (whose dimension match the number of
non-zero in the Jacobian matrix).

"""
function jacobian_structure end

"""
    jacobian!(nlp::AbstractNLPEvaluator, jac::AbstractMatrix, u)

Evaluate the Jacobian of the constraints, at variable `u`.
Store the result inplace, in the `m x n` dense matrix `jac`.

"""
function jacobian! end

"""
    jacobian_coo!(nlp::AbstractNLPEvaluator, jac::AbstractVector, u)

Evaluate the (sparse) Jacobian of the constraints at variable `u`
in COO format.
Store the result inplace, in the `nnzj` vector `jac`.

"""
function jacobian_coo! end

"""
    jprod!(nlp::AbstractNLPEvaluator, jv, u, v)

Evaluate the Jacobian-vector product ``J v`` of the constraints.
The vector `jv` is modified inplace.

Let `(n, m) = n_variables(nlp), n_constraints(nlp)`.

* `u` is a vector with dimension `n`
* `v` is a vector with dimension `n`
* `jv` is a vector with dimension `m`

"""
function jprod! end

"""
    jtprod!(nlp::AbstractNLPEvaluator, jv, u, v)

Evaluate the transpose Jacobian-vector product ``J^{T} v`` of the constraints.
The vector `jv` is modified inplace.

Let `(n, m) = n_variables(nlp), n_constraints(nlp)`.

* `u` is a vector with dimension `n`
* `v` is a vector with dimension `m`
* `jv` is a vector with dimension `n`

"""
function jtprod! end

"""
    ojtprod!(nlp::AbstractNLPEvaluator, jv, u, σ, v)

Evaluate the transpose Jacobian-vector product `J' * [σ ; v]`,
with `J` the Jacobian of the vector `[f(x); h(x)]`.
`f(x)` is the current objective and `h(x)` constraints.
The vector `jv` is modified inplace.

Let `(n, m) = n_variables(nlp), n_constraints(nlp)`.

* `jv` is a vector with dimension `n`
* `u` is a vector with dimension `n`
* `σ` is a scalar
* `v` is a vector with dimension `m`

"""
function ojtprod! end

"""
    hessian!(nlp::AbstractNLPEvaluator, H, u)

Evaluate the Hessian `∇²f(u)` of the objective function `f(u)`.
Store the result inplace, in the `n x n` dense matrix `H`.

"""
function hessian! end

"""
    hessian_coo!(nlp::AbstractNLPEvaluator, hess::AbstractVector, u)

Evaluate the (sparse) Hessian of the constraints at variable `u`
in COO format.
Store the result inplace, in the `nnzh` vector `hess`.

"""
function hessian_coo! end

"""
    hessprod!(nlp::AbstractNLPEvaluator, hessvec, u, v)

Evaluate the Hessian-vector product `∇²f(u) * v` of the objective
evaluated at variable `u`.
Store the result inplace, in the vector `hessvec`.

## Note
The vector `hessprod` should have the same length as `u`.

"""
function hessprod! end

@doc raw"""
    hessian_lagrangian_prod!(nlp::AbstractNLPEvaluator, hessvec, u, y, σ, v)

Evaluate the Hessian-vector product of the Lagrangian
function ``L(u, y) = f(u) + \sum_i y_i c_i(u)`` with a vector `v`:
```math
∇²L(u, y) ⋅ v  = σ ∇²f(u) ⋅ v + \sum_i y_i ∇²c_i(u) ⋅ v
```

Store the result inplace, in the vector `hessvec`.

### Arguments

* `hessvec` is a `AbstractVector` with dimension `n`, which is modified inplace.
* `u` is a `AbstractVector` with dimension `n`, storing the current variable.
* `y` is a `AbstractVector` with dimension `n`, storing the current constraints' multipliers
* `σ` is a scalar, encoding the objective's scaling
* `v` is a vector with dimension `n`.
"""
function hessian_lagrangian_prod! end

@doc raw"""
    hessian_lagrangian_penalty_prod!(nlp::AbstractNLPEvaluator, hessvec, u, y, σ, d, v)

Evaluate the Hessian-vector product of the Augmented Lagrangian
function ``L(u, y) = f(u) + \sum_i y_i c_i(u) + \frac{1}{2} d_i c_i(u)^2`` with a vector `v`:
```math
∇²L(u, y) ⋅ v  = σ ∇²f(u) ⋅ v + \sum_i (y_i + d_i) ∇²c_i(u) ⋅ v + \sum_i d_i ∇c_i(u)^T ∇c_i(u)
```

Store the result inplace, in the vector `hessvec`.

### Arguments

* `hessvec` is a `AbstractVector` with dimension `n`, which is modified inplace.
* `u` is a `AbstractVector` with dimension `n`, storing the current variable.
* `y` is a `AbstractVector` with dimension `n`, storing the current constraints' multipliers
* `σ` is a scalar
* `v` is a vector with dimension `n`.
* `d` is a vector with dimension `m`.
"""
function hessian_lagrangian_penalty_prod! end

"Return `true` if the problem is constrained, `false` otherwise."
is_constrained(nlp::AbstractNLPEvaluator) = n_constraints(nlp) > 0

"""
    constraints_type(nlp::AbstractNLPEvaluator)

Return the type of the non-linear constraints of the evaluator `nlp`,
as a `Symbol`. Result could be `:inequality` if problem has only
inequality constraints, `:equality` if problem has only equality
constraints, or `:mixed` if problem has both types of constraints.
"""
function constraints_type end

"Check if Hessian of objective is implemented."
has_hessian(nlp::AbstractNLPEvaluator) = false
"Check if Hessian of Lagrangian is implemented."
has_hessian_lagrangian(nlp::AbstractNLPEvaluator) = false

"""
    reset!(nlp::AbstractNLPEvaluator)

Reset evaluator `nlp` to default configuration.

"""
function reset! end

include("common.jl")
include("scalers.jl")

# Basic evaluators
include("reduced_evaluator.jl")
include("full_space_evaluator.jl")
include("stoch_evaluator.jl")
include("slack_evaluator.jl")
include("feasibility_evaluator.jl")
# include("proxal_evaluators.jl")
include("bridge_evaluator.jl")

# Penalty evaluators
include("penalty.jl")
include("auglag.jl")
include("qp_model.jl")

