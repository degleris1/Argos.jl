
using LazyArtifacts
using DelimitedFiles

DEMANDS = joinpath(artifact"ExaData", "ExaData", "mp_demand")

function build_scopf_model(polar, buffer, solver; line_constraints=true)
    nbus = ExaPF.get(polar, PS.NumberOfBuses())
    ngen = ExaPF.get(polar, PS.NumberOfGenerators())
    nlines = ExaPF.get(polar, PS.NumberOfLines())

    network = polar.network
    pq, pv, ref = network.pq, network.pv, network.ref

    K = div(length(buffer.pload), nbus)

    pf = polar.network
    baseMVA = pf.baseMVA

    pg_min, pg_max = PS.bounds(pf, PS.Generators(), PS.ActivePower())
    qg_min, qg_max = PS.bounds(pf, PS.Generators(), PS.ReactivePower())
    vm_min, vm_max = PS.bounds(pf, PS.Buses(), PS.VoltageMagnitude())

    flow_min, flow_max = PS.bounds(pf, PS.Lines(), PS.ActivePower())

    flow_max = min.(1e5, flow_max)

    vm0 = buffer.vmag
    va0 = buffer.vang
    pg0 = buffer.pgen

    yff_re = real.(pf.lines.Yff)
    yff_im = imag.(pf.lines.Yff)
    yft_re = real.(pf.lines.Yft)
    yft_im = imag.(pf.lines.Yft)
    ytf_re = real.(pf.lines.Ytf)
    ytf_im = imag.(pf.lines.Ytf)
    ytt_re = real.(pf.lines.Ytt)
    ytt_im = imag.(pf.lines.Ytt)

    # Pd = buffer.pload
    # Qd = buffer.qload

    cost_coefs = PS.get_costs_coefficients(pf)

    bus2gen = PS.get_bus_generators(pf.buses, pf.generators, pf.bus_to_indexes)
    pvgen = Int[]
    for g in 1:ngen
        if network.gen2bus != ref[1]
            push!(pvgen, g)
        end
    end

    # Power flow data
    Ybus = pf.Ybus
    rows = Ybus.rowval
    yvals = Ybus.nzval
    g_ij = real.(yvals)
    b_ij = imag.(yvals)

    #=
        Build model
    =#

    opfmodel = Model(solver)

    # VARIABLES
    @variable(opfmodel, pg_min[i] <= Pg[i=1:ngen, k=1:K] <= pg_max[i], start=pg0[i])
    @variable(opfmodel, qg_min[i] <= Qg[i=1:ngen, k=1:K] <= qg_max[i])
    @variable(opfmodel, vm_min[i] <= Vm[i=1:nbus, k=1:K] <= vm_max[i], start=vm0[i])
    @variable(opfmodel, Va[i=1:nbus, k=1:K], start=va0[i])

    # @variable(opfmodel, Pd[i=1:nbus])
    # @variable(opfmodel, Qd[i=1:nbus])
    # JuMP.fix.(Pd, buffer.pload)
    # JuMP.fix.(Qd, buffer.qload)
    Pd = reshape(buffer.pload, nbus, K)
    Qd = reshape(buffer.qload, nbus, K)

    # Power-flow constraints
    ## active
    opfmodel.ext[:active_pf] = @NLconstraint(
        opfmodel, [b=1:nbus, k=1:K],
        Vm[b, k] * sum(
            Vm[rows[c], k] * (g_ij[c] * cos(Va[b, k] - Va[rows[c], k]) + b_ij[c] * sin(Va[b, k] - Va[rows[c], k]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)
        ) == (sum(Pg[g, k] for g in get(bus2gen, b, Int[])) - Pd[b, k])
    )
    ## reactive
    opfmodel.ext[:reactive_pf] = @NLconstraint(
        opfmodel, [b=1:nbus, k=1:K],
        Vm[b, k] * sum(
            Vm[rows[c], k] * (g_ij[c] * sin(Va[b, k] - Va[rows[c], k]) - b_ij[c] * cos(Va[b, k] - Va[rows[c], k]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)) == (sum(Qg[g, k] for g in get(bus2gen, b, Int[])) - Qd[b, k])
    )

    # Line constraints
    if line_constraints
        f = pf.lines.from_buses
        t = pf.lines.to_buses
        ## from lines
        yff_abs = yff_re.^2 .+ yff_im.^2
        yft_abs = yft_re.^2 .+ yft_im.^2
        yre_fr =   yff_re .* yft_re .+ yff_im .* yft_im
        yim_fr = - yff_re .* yft_im .+ yff_im .* yft_re

        opfmodel.ext[:line_fr] = @NLconstraint(
            opfmodel, [ℓ=1:nlines, k=1:K],
            Vm[f[ℓ], k]^2 * (
                yff_abs[ℓ] * Vm[f[ℓ], k]^2 + yft_abs[ℓ] * Vm[t[ℓ], k]^2 +
                2.0 * Vm[f[ℓ], k] * Vm[t[ℓ], k] * (yre_fr[ℓ] * cos(Va[f[ℓ], k]-Va[t[ℓ], k]) - yim_fr[ℓ] * sin(Va[f[ℓ], k]-Va[t[ℓ], k]))
            ) <= flow_max[ℓ]
        )

        ## to lines
        ytf_abs = ytf_re.^2 .+ ytf_im.^2
        ytt_abs = ytt_re.^2 .+ ytt_im.^2
        yre_to =   ytf_re .* ytt_re .+ ytf_im .* ytt_im
        yim_to = - ytf_re .* ytt_im .+ ytf_im .* ytt_re

        opfmodel.ext[:line_to] = @NLconstraint(
            opfmodel, [ℓ=1:nlines, k=1:K],
            Vm[t[ℓ], k]^2 * (
                ytf_abs[ℓ] * Vm[f[ℓ], k]^2 + ytt_abs[ℓ] * Vm[t[ℓ], k]^2 +
                2.0 * Vm[f[ℓ], k] * Vm[t[ℓ], k] * (yre_to[ℓ] * cos(Va[f[ℓ], k]-Va[t[ℓ], k]) - yim_to[ℓ] * sin(Va[f[ℓ], k]-Va[t[ℓ], k]))
            ) <= flow_max[ℓ]
        )
    end

    # Coupling constraints
    # for k in 2:K
    #     @constraint(opfmodel, Vm[pv, k] .== Vm[pv, 1])
    #     @constraint(opfmodel, Vm[ref, k] .== Vm[ref, 1])
    #     @constraint(opfmodel, Pg[pvgen, k] .== Pg[pvgen, 1])
    # end

    # Objective
    @objective(
        opfmodel,
        Min,
        1e-3 / K * sum(
            cost_coefs[g, 4] * Pg[g, k]^2 + cost_coefs[g, 3] * Pg[g, k] + cost_coefs[g, 2]
            for g in 1:ngen, k=1:K
        )
    )

    opfmodel.ext[:exapf] = pf

    return opfmodel
end

function scadnlp(datafile::String, demands::String, nscen)
    polar = ExaPF.PolarForm(datafile)

    pload = readdlm(joinpath(DEMANDS, "$(demands)_onehour_60.Pd"))[:, 1:nscen] ./ 100.0
    qload = readdlm(joinpath(DEMANDS, "$(demands)_onehour_60.Qd"))[:, 1:nscen] ./ 100.0

    stack = ExaPF.BlockNetworkStack(polar, pload, qload)

    optimizer = () -> MadNLP.Optimizer(
        linear_solver=MadNLPMa27,
        print_level=MadNLP.DEBUG,
        nlp_scaling=true,
        tol=1e-10,
        max_iter=10,
        dual_initialized=true,
    )
    m = build_scopf_model(polar, stack, optimizer; line_constraints=true)
    # attach_callback!(m)
    JuMP.optimize!(m; _differentiation_backend = SymbolicAD.DefaultBackend())
    # JuMP.optimize!(m)
    return m
end

function scopf(datafile::String, demands::String, nscen)
    polar = ExaPF.PolarForm(datafile)
    pload = readdlm(joinpath(DEMANDS, "$(demands)_onehour_60.Pd"))[:, 1:nscen] ./ 100.0
    qload = readdlm(joinpath(DEMANDS, "$(demands)_onehour_60.Qd"))[:, 1:nscen] ./ 100.0

    stack = ExaPF.BlockNetworkStack(polar, pload, qload)
    m = build_scopf_model(polar, stack, Ipopt.Optimizer; line_constraints=true)
    # attach_callback!(m)
    #
    JuMP.optimize!(m)
    return m
end
