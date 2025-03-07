
@enum(ProxALTime,
    Origin,
    Final,
    Normal,
)

abstract type AbstractTimeStep end
struct Current <: AbstractTimeStep end
struct Previous <: AbstractTimeStep end
struct Next <: AbstractTimeStep end


"""
    ProxALEvaluator{T, VI, VT, MT} <: AbstractNLPEvaluator

Evaluator wrapping a `ReducedSpaceEvaluator` for use inside the
decomposition algorithm implemented in [ProxAL.jl](https://github.com/exanauts/ProxAL.jl).

"""
mutable struct ProxALEvaluator{T, VI, VT, MT, Pullback, Hess} <: AbstractNLPEvaluator
    inner::ReducedSpaceEvaluator{T, VI, VT, MT}
    obj_stack::Pullback
    hessian_obj::Hess
    s_min::VT
    s_max::VT
    nu::Int
    ng::Int
    # Augmented penalties parameters
    time::ProxALTime
    scale_objective::T
    τ::T
    λf::VT
    λt::VT
    ρf::T
    ρt::T
    pg_f::VT
    pg_ref::VT
    pg_t::VT
end
function ProxALEvaluator(
    nlp::ReducedSpaceEvaluator{T, VI, VT, MT},
    time::ProxALTime;
    τ=0.1, ρf=0.1, ρt=0.1, scale_obj=1.0, want_hessian=true,
) where {T, VI, VT, MT}
    nu = n_variables(nlp)
    ng = Base.get(nlp, PS.NumberOfGenerators())

    s_min = fill!(VT(undef, ng), zero(T))
    s_max = fill!(VT(undef, ng), one(T))
    λf = fill!(VT(undef, ng), zero(T))
    λt = fill!(VT(undef, ng), zero(T))

    pgf = fill!(VT(undef, ng), zero(T))
    pgc = fill!(VT(undef, ng), zero(T))
    pgt = fill!(VT(undef, ng), zero(T))

    pbm = ExaPF.pullback_ramping(nlp.model, nothing)

    hess = nothing
    if want_hessian
        hess = AutoDiff.Hessian(nlp.model, ExaPF.cost_penalty_ramping_constraints; tape=pbm)
    end
    return ProxALEvaluator(
        nlp, pbm, hess, s_min, s_max, nu, ng, time, scale_obj,
        τ, λf, λt, ρf, ρt,
        pgf, pgc, pgt,
    )
end
function ProxALEvaluator(
    pf::PS.PowerNetwork,
    time::ProxALTime;
    device=ExaPF.CPU(),
    options...
)
    # Build network polar formulation
    model = PolarForm(pf, device)
    # Build reduced space evaluator
    nlp = ReducedSpaceEvaluator(model; options...)
    return ProxALEvaluator(nlp, time)
end
function ProxALEvaluator(
    datafile::String;
    time::ProxALTime=Normal,
    device=ExaPF.CPU(),
    options...
)
    nlp = ReducedSpaceEvaluator(datafile; device=device, options...)
    return ProxALEvaluator(nlp, time; options...)
end

n_variables(nlp::ProxALEvaluator) = nlp.nu + nlp.ng
n_constraints(nlp::ProxALEvaluator) = n_constraints(nlp.inner)


constraints_type(::ProxALEvaluator) = :inequality
has_hessian(::ProxALEvaluator) = true
backend(nlp::ProxALEvaluator) = backend(nlp.inner)
inner_evaluator(nlp::ProxALEvaluator) = inner_evaluator(nlp.inner)

# Getters
Base.get(nlp::ProxALEvaluator, attr::AbstractNLPAttribute) = Base.get(nlp.inner, attr)
Base.get(nlp::ProxALEvaluator, attr::ExaPF.AbstractVariable) = Base.get(nlp.inner, attr)
Base.get(nlp::ProxALEvaluator, attr::PS.AbstractNetworkValues) = Base.get(nlp.inner, attr)
Base.get(nlp::ProxALEvaluator, attr::PS.AbstractNetworkAttribute) = Base.get(nlp.inner, attr)

# Setters
function setvalues!(nlp::ProxALEvaluator, attr::PS.AbstractNetworkValues, values)
    setvalues!(nlp.inner, attr, values)
end
function transfer!(nlp::ProxALEvaluator, vm, va, pg, qg)
    transfer!(nlp.inner, vm, va, pg, qg)
end

# Initial position
function initial(nlp::ProxALEvaluator)
    u0 = initial(nlp.inner)
    s0 = copy(nlp.s_min)
    return [u0; s0]
end

# Bounds
function bounds(nlp::ProxALEvaluator, ::Variables)
    u♭, u♯ = bounds(nlp.inner, Variables())
    return [u♭; nlp.s_min], [u♯; nlp.s_max]
end
bounds(nlp::ProxALEvaluator, ::Constraints) = bounds(nlp.inner, Constraints())

function update!(nlp::ProxALEvaluator, w)
    u = @view w[1:nlp.nu]
    s = @view w[nlp.nu+1:end]
    conv = update!(nlp.inner, u)
    pg = Base.get(nlp.inner, PS.ActivePower())
    return conv
end

## Update penalty terms
function update_multipliers!(nlp::ProxALEvaluator, ::Current, λt)
    copyto!(nlp.λf, λt)
end
function update_multipliers!(nlp::ProxALEvaluator, ::Next, λt)
    copyto!(nlp.λt, λt)
end

function update_primal!(nlp::ProxALEvaluator, ::Previous, pgk)
    copyto!(nlp.pg_f, pgk)
end
function update_primal!(nlp::ProxALEvaluator, ::Current, pgk)
    copyto!(nlp.pg_ref, pgk)
end
function update_primal!(nlp::ProxALEvaluator, ::Next, pgk)
    copyto!(nlp.pg_t, pgk)
end

## Objective
function objective(nlp::ProxALEvaluator, w)
    u = @view w[1:nlp.nu]
    s = @view w[nlp.nu+1:end]
    model = nlp.inner.model
    buffer = Base.get(nlp.inner, ExaPF.PhysicalState())
    return ExaPF.cost_penalty_ramping_constraints(
        model, buffer, s, Int(nlp.time),
        nlp.scale_objective, nlp.τ, nlp.λf, nlp.λt, nlp.ρf, nlp.ρt, nlp.pg_f, nlp.pg_ref, nlp.pg_t
    )
end

## Gradient
function _gradient_control!(nlp::ProxALEvaluator, grad, w)
    # Gradient wrt u
    gu = @view grad[1:nlp.nu]
    u = @view w[1:nlp.nu]
    s = @view w[nlp.nu+1:end]
    buffer = Base.get(nlp.inner, ExaPF.PhysicalState())
    ∂obj = nlp.obj_stack
    ju = ∂obj.stack.∇fᵤ ; jx = ∂obj.stack.∇fₓ
    # Evaluate adjoint of cost function and update inplace AdjointStackObjective
    ExaPF.adjoint_penalty_ramping_constraints!(
        nlp.inner.model, ∂obj, buffer, s, Int(nlp.time),
        nlp.scale_objective, nlp.τ, nlp.λf, nlp.λt, nlp.ρf, nlp.ρt, nlp.pg_f, nlp.pg_ref, nlp.pg_t
    )
    _adjoint_solve!(nlp.inner, gu, jx, ju, nlp.inner.λ, w)
    nlp.inner.is_adjoint_objective_updated = true
end

function _gradient_slack!(nlp::ProxALEvaluator, grad, w)
    s = @view w[nlp.nu+1:end]
    pg = Base.get(nlp.inner, PS.ActivePower())
    # Gradient wrt s
    gs = @view grad[nlp.nu+1:end]
    if nlp.time != Origin
        gs .= nlp.λf .+ nlp.ρf .* (nlp.pg_f .- pg .+ s)
    else
        gs .= 0.0
    end
end

function gradient!(nlp::ProxALEvaluator, g, w)
    _gradient_control!(nlp, g, w)
    _gradient_slack!(nlp, g, w)
    return nothing
end

function constraint!(nlp::ProxALEvaluator, cons, w)
    u = @view w[1:nlp.nu]
    constraint!(nlp.inner, cons, u)
end

function jacobian_structure(nlp::ProxALEvaluator)
    m, n = n_constraints(nlp), n_variables(nlp)
    nnzj = m * n
    rows = zeros(Int, nnzj)
    cols = zeros(Int, nnzj)
    jacobian_structure!(nlp, rows, cols)
    return rows, cols
end

## Jacobian
function jacobian_structure!(nlp::ProxALEvaluator, rows, cols)
    m, n = n_constraints(nlp), n_variables(nlp)
    idx = 1
    for i in 1:n # number of variables
        for c in 1:m #number of constraints
            rows[idx] = c ; cols[idx] = i
            idx += 1
        end
    end
end

function jacobian!(nlp::ProxALEvaluator, jac, w)
    m = n_constraints(nlp)
    u = @view w[1:nlp.nu]
    Jᵤ = @view jac[:, 1:nlp.nu]
    jacobian!(nlp.inner, Jᵤ, u)
end

function jprod!(nlp::ProxALEvaluator, jv, w, v)
    u = @view w[1:nlp.nu]
    vu = @view v[1:nlp.nu]
    jprod!(nlp.inner, jv, u, vu)
end

## Transpose Jacobian-vector product
## ProxAL does not add any constraint to the reduced model
function jtprod!(nlp::ProxALEvaluator, jv, w, v)
    u = @view w[1:nlp.nu]
    jvu = @view jv[1:nlp.nu]
    jtprod!(nlp.inner, jvu, u, v)
end

function ojtprod!(nlp::ProxALEvaluator, jv, w, σ, v)
    gu = @view jv[1:nlp.nu]
    s = @view w[nlp.nu+1:end]
    u = @view w[1:nlp.nu]

    # w.r.t. u
    buffer = Base.get(nlp.inner, ExaPF.PhysicalState())
    ∂obj = nlp.obj_stack
    jvx = ∂obj.stack.jvₓ ; fill!(jvx, 0)
    jvu = ∂obj.stack.jvᵤ ; fill!(jvu, 0)

    # Evaluate adjoint of cost function and update inplace AdjointStackObjective
    ExaPF.adjoint_penalty_ramping_constraints!(
        nlp.inner.model, ∂obj, buffer, s, Int(nlp.time),
        nlp.scale_objective, nlp.τ, nlp.λf, nlp.λt, nlp.ρf, nlp.ρt, nlp.pg_f, nlp.pg_ref, nlp.pg_t
    )
    copyto!(jvx, ∂obj.stack.∇fₓ)
    copyto!(jvu, ∂obj.stack.∇fᵤ)
    # compute gradient of objective
    jvx .*= σ
    jvu .*= σ
    # compute transpose Jacobian vector product of constraints
    full_jtprod!(nlp.inner, jvx, jvu, u, v)
    # Evaluate gradient in reduced space
    _adjoint_solve!(nlp.inner, gu, jvx, jvu, nlp.inner.μ, w)
    nlp.inner.is_adjoint_lagrangian_updated = true

    # w.r.t. s
    _gradient_slack!(nlp, jv, w)
end

#=
    For ProxAL, we have:
    H = [ H_xx  H_ux  J_x' ]
        [ H_xu  H_uu  J_u' ]
        [ J_x   J_u   ρ I  ]

    so, if `v = [v_x; v_u; v_s]`, we get

    H * v = [ H_xx v_x  +   H_ux v_u  +  J_x' v_s ]
            [ H_xu v_x  +   H_uu v_u  +  J_u' v_s ]
            [  J_x v_x  +    J_u v_u  +   ρ I     ]

=#
function hessprod!(nlp::ProxALEvaluator, hessvec, w, v)
    @assert nlp.inner.has_hessian
    @assert nlp.inner.has_jacobian

    model = nlp.inner.model
    nx = ExaPF.get(model, ExaPF.NumberOfState())
    nu = ExaPF.get(model, ExaPF.NumberOfControl())

    u = @view w[1:nlp.nu]
    vᵤ = @view v[1:nlp.nu]
    vₛ = @view v[1+nlp.nu:end]

    fill!(hessvec, 0.0)

    hvu = @view hessvec[1:nlp.nu]
    hvs = @view hessvec[1+nlp.nu:end]

    ## OBJECTIVE HESSIAN
    σ = 1.0
    hessprod!(nlp.inner, hvu, u, vᵤ)

    # Contribution of slack node
    if nlp.time != Origin
        hvs .+= nlp.ρf .* vₛ
        # TODO: implement block corresponding to Jacobian
        # and transpose-Jacobian
    end
    return
end

## Utils function
function reset!(nlp::ProxALEvaluator)
    reset!(nlp.inner)
    # Reset multipliers
    fill!(nlp.λf, 0)
    fill!(nlp.λt, 0)
    # Reset proximal centers
    fill!(nlp.pg_f, 0)
    fill!(nlp.pg_ref, 0)
    fill!(nlp.pg_t, 0)
end

function primal_infeasibility!(nlp::ProxALEvaluator, cons, w)
    @assert length(w) == nlp.nu + nlp.ng
    u = @view w[1:nlp.nu]
    return primal_infeasibility!(nlp.inner, cons, u)
end
function primal_infeasibility(nlp::ProxALEvaluator, w)
    @assert length(w) == nlp.nu + nlp.ng
    u = @view w[1:nlp.nu]
    return primal_infeasibility(nlp.inner, u)
end
