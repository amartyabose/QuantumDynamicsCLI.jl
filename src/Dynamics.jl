module Dynamics

using LinearAlgebra
using QuantumDynamics
using ..QDSimUtilities, ..ParseInput

const PILD_reference = """
- A. Bose, “Incorporation of Empirical Gain and Loss Mechanisms in Open Quantum Systems through Path Integral Lindblad Dynamics,” J. Phys. Chem. Lett. 15(12), 3363–3368 (2024)."""

function dynamics(::QDSimUtilities.Method"TEMPO-TTM", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a TEMPO dynamics calculation with TTM. Please cite:"
        QDSimUtilities.print_citation(TEMPO.references)
        QDSimUtilities.print_citation(TTM.references)
    end
    rmax = sim_node["rmax"]
    rmax_group = Utilities.create_and_select_group(dt_group, "rmax=$(rmax)")
    kmax::Union{Nothing,Int} = get(sim_node, "kmax", nothing)
    cutoff = get(sim_node, "cutoff", 1e-10)
    maxdim = get(sim_node, "maxdim", 1000)
    algorithm = get(sim_node, "algorithm", "naive")
    @info "Running with $(BLAS.get_num_threads()) threads."

    if !isnothing(kmax)
        @assert kmax <= rmax "kmax = $(kmax) should be less than rmax = $(rmax)."
        rmax_group = Utilities.create_and_select_group(rmax_group, "kmax=$(kmax)")
    end
    maxdim_group = Utilities.create_and_select_group(rmax_group, "maxdim=$(maxdim)")
    cutoff_group = Utilities.create_and_select_group(maxdim_group, "cutoff=$(cutoff)")
    data = Utilities.create_and_select_group(cutoff_group, "algorithm=$(algorithm)")

    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:rmax*sim.dt/units.time_unit |> collect)
        flush(data)
        extraargs = TEMPO.TEMPOArgs(; cutoff, maxdim, algorithm)
        if haskey(sim_node, "lindblad")
            @info "Using the PILD method."
            QDSimUtilities.print_citation(PILD_reference)
            decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
            L = [ParseInput.parse_operator(sim_node["lindblad"][i], sys.Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
        else
            L = nothing
        end

        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, L)
        Utilities.check_or_insert_value(data, "fbU", fbU)
        flush(data)
        TEMPO.build_augmented_propagator(; fbU, Jw=bath.Jw, β=bath.β, dt=sim.dt, ntimes=rmax, kmax, extraargs, svec=bath.svecs, verbose=true, output=data)
        @info "After this run, please run a propagate-using-tmats calculation to obtain the time evolution of a particular density matrix."
    end
    data
end

function dynamics(::QDSimUtilities.Method"QuAPI-TTM", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a QuAPI calculation with TTM. Please cite:"
        QDSimUtilities.print_citation(QuAPI.references)
        QDSimUtilities.print_citation(TTM.references)
    end
    rmax = sim_node["rmax"]
    rmax_group = Utilities.create_and_select_group(dt_group, "rmax=$(rmax)")
    cutoff = get(sim_node, "cutoff", 1e-10)
    data = Utilities.create_and_select_group(rmax_group, "cutoff=$(cutoff)")
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:rmax*sim.dt/units.time_unit |> collect)
        flush(data)

        extraargs = QuAPI.QuAPIArgs(; cutoff)
        if haskey(sim_node, "lindblad")
            @info "Using the PILD method."
            QDSimUtilities.print_citation(PILD_reference)
            decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
            L = [ParseInput.parse_operator(sim_node["lindblad"][i], sys.Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
        else
            L = nothing
        end

        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, L)
        Utilities.check_or_insert_value(data, "fbU", fbU)
        flush(data)
        QuAPI.build_augmented_propagator(; fbU, Jw=bath.Jw, β=bath.β, dt=sim.dt, ntimes=rmax, extraargs, svec=bath.svecs, verbose=true, output=data, exec=QDSimUtilities.parse_exec(exec))
        @info "After this run, please run a propagate-using-tmats calculation to obtain the time evolution of a particular density matrix."
    end
    data
end

function dynamics(::QDSimUtilities.Method"adaptive-kinks-QuAPI-TTM", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running an adaptive kinks QuAPI calculation with TTM. Please cite:"
        QDSimUtilities.print_citation(QuAPI.references)
        QDSimUtilities.print_citation(TTM.references)
    end
    rmax = sim_node["rmax"]
    rmax_group = Utilities.create_and_select_group(dt_group, "rmax=$(rmax)")
    cutoff = get(sim_node, "cutoff", 1e-10)
    cutoff_group = Utilities.create_and_select_group(rmax_group, "cutoff=$(cutoff)")
    num_kinks = get(sim_node, "num_kinks", -1)
    kink_group = Utilities.create_and_select_group(cutoff_group, "num_kinks=$(num_kinks)")
    num_blips = get(sim_node, "num_blips", -1)
    blip_group = Utilities.create_and_select_group(kink_group, "num_blips=$(num_blips)")
    prop_cutoff = get(sim_node, "propagator_cutoff", 0.0)
    data = Utilities.create_and_select_group(blip_group, "prop_cutoff=$(prop_cutoff)")
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:rmax*sim.dt/units.time_unit |> collect)
        flush(data)

        extraargs = QuAPI.QuAPIArgs(; cutoff, prop_cutoff, num_kinks, num_blips)
        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, forward_backward=true)
        Utilities.check_or_insert_value(data, "fbU", fbU)
        flush(data)
        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, forward_backward=false)
        QuAPI.build_augmented_propagator_kink(; fbU, Jw=bath.Jw, β=bath.β, dt=sim.dt, ntimes=rmax, extraargs, svec=bath.svecs, verbose=true, output=data, exec=QDSimUtilities.parse_exec(exec))
        @info "After this run, please run a propagate-using-tmats calculation to obtain the time evolution of a particular density matrix."
    end
    data
end

function dynamics(::QDSimUtilities.Method"Blip-TTM", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a Blip calculation with TTM. Please cite:"
        QDSimUtilities.print_citation(Blip.references)
        QDSimUtilities.print_citation(TTM.references)
    end
    rmax = sim_node["rmax"]
    rmax_group = Utilities.create_and_select_group(dt_group, "rmax=$(rmax)")
    max_blips = get(sim_node, "max_blips", -1)
    num_changes = get(sim_node, "num_changes", -1)
    if max_blips == -1
        mb_data = Utilities.create_and_select_group(rmax_group, "max_blips=all")
    else
        mb_data = Utilities.create_and_select_group(rmax_group, "max_blips=$(max_blips)")
    end
    if num_changes == -1
        data = Utilities.create_and_select_group(mb_data, "num_changes=all")
    else
        data = Utilities.create_and_select_group(mb_data, "num_changes=$(num_changes)")
    end
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:rmax*sim.dt/units.time_unit |> collect)
        flush(data)

        extraargs = Blip.BlipArgs(; max_blips, num_changes)
        if haskey(sim_node, "lindblad")
            @info "Using the PILD method."
            QDSimUtilities.print_citation(PILD_reference)
            decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
            L = [ParseInput.parse_operator(sim_node["lindblad"][i], sys.Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
        else
            L = nothing
        end

        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, L)
        Utilities.check_or_insert_value(data, "fbU", fbU)
        flush(data)
        Blip.build_augmented_propagator(; fbU, Jw=bath.Jw, β=bath.β, dt=sim.dt, ntimes=rmax, extraargs, svec=bath.svecs, verbose=true, output=data, exec=QDSimUtilities.parse_exec(exec))
    end
    data
end

function dynamics(::QDSimUtilities.Method"adaptive-kinks-QuAPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running an adaptive kinks QuAPI calculation with TTM. Please cite:"
        QDSimUtilities.print_citation(QuAPI.references)
    end
    rmax = sim_node["rmax"]
    rmax_group = Utilities.create_and_select_group(dt_group, "rmax=$(rmax)")
    cutoff = get(sim_node, "cutoff", 1e-10)
    cutoff_group = Utilities.create_and_select_group(rmax_group, "cutoff=$(cutoff)")
    num_kinks = get(sim_node, "num_kinks", -1)
    kink_group = Utilities.create_and_select_group(cutoff_group, "num_kinks=$(num_kinks)")
    num_blips = get(sim_node, "num_blips", -1)
    blip_group = Utilities.create_and_select_group(kink_group, "num_blips=$(num_blips)")
    prop_cutoff = get(sim_node, "propagator_cutoff", 0.0)
    data = Utilities.create_and_select_group(blip_group, "prop_cutoff=$(prop_cutoff)")
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end
    outgroup = sim_node["outgroup"]
    if !dry
        data = Utilities.create_and_select_group(data, outgroup)
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:rmax*sim.dt/units.time_unit |> collect)
        flush(data)

        extraargs = QuAPI.QuAPIArgs(; cutoff, prop_cutoff, num_kinks, num_blips)
        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, forward_backward=true)
        Utilities.check_or_insert_value(data, "fbU", fbU)
        fbU = Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=rmax, forward_backward=false)
        flush(data)
        ρ0 = ParseInput.parse_operator(sim_node["rho0"], sys.Hamiltonian)
        QuAPI.propagate_kink(; fbU, Jw=bath.Jw, β=bath.β, ρ0, dt=sim.dt, ntimes=rmax, extraargs, svec=bath.svecs, verbose=true, output=data, exec=QDSimUtilities.parse_exec(exec))
    end
    data
end

function dynamics(::QDSimUtilities.Method"TEMPO", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a TEMPO dynamics calculation. Please cite:"
        QDSimUtilities.print_citation(TEMPO.references)
    end
    kmax = sim_node["kmax"]
    cutoff = get(sim_node, "cutoff", 1e-10)
    maxdim = get(sim_node, "maxdim", 1000)
    algorithm = get(sim_node, "algorithm", "naive")
    @info "Running with $(BLAS.get_num_threads()) threads."

    kmax_group = Utilities.create_and_select_group(dt_group, "kmax=$(kmax)")
    maxdim_group = Utilities.create_and_select_group(kmax_group, "maxdim=$(maxdim)")
    cutoff_group = Utilities.create_and_select_group(maxdim_group, "cutoff=$(cutoff)")
    data = Utilities.create_and_select_group(cutoff_group, "algorithm=$(algorithm)")

    outgroup = sim_node["outgroup"]

    if !dry
        outgrouphdf5 = Utilities.create_and_select_group(data, outgroup)
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect)
        Utilities.check_or_insert_value(outgrouphdf5, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(outgrouphdf5, "time", 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect)
        flush(data)
        extraargs = TEMPO.TEMPOArgs(; cutoff, maxdim, algorithm)
        if haskey(sim_node, "lindblad")
            @info "Using the PILD method."
            QDSimUtilities.print_citation(PILD_reference)
            decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
            L = [ParseInput.parse_operator(sim_node["lindblad"][i], sys.Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
        else
            L = nothing
        end

        fbU = if isnothing(sys.external_fields)
            Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=sim.nsteps, L)
        else
            ndivs = get(sim_node, "ndivs", 100)
            Propagators.calculate_bare_propagators(; Hamiltonian=sys.Hamiltonian, dt=sim.dt, ntimes=sim.nsteps, L, ndivs, external_fields=sys.external_fields)
        end
        Utilities.check_or_insert_value(data, "fbU", fbU)
        flush(data)
        ρ0 = ParseInput.parse_operator(sim_node["rho0"], sys.Hamiltonian)
        TEMPO.propagate(; fbU, Jw=bath.Jw, β=bath.β, ρ0, dt=sim.dt, ntimes=sim.nsteps, kmax, extraargs, svec=bath.svecs, verbose=true, output=data, outgroup=outgroup)
    end
    data
end

function dynamics(::QDSimUtilities.Method"HEOM", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a HEOM calculation."
        QDSimUtilities.print_citation(HEOM.references)
    end
    num_modes = sim_node["num_modes"]
    Lmax = sim_node["Lmax"]
    num_modes_group = Utilities.create_and_select_group(dt_group, "num_modes=$(num_modes)")
    Lmax_group = Utilities.create_and_select_group(num_modes_group, "Lmax=$(Lmax)")
    threshold = get(sim_node, "threshold", 0.0)
    threshold_group = Utilities.create_and_select_group(Lmax_group, "threshold=$(threshold)")
    reltol = get(sim_node, "reltol", 1e-6)
    abstol = get(sim_node, "abstol", 1e-6)
    data = Utilities.create_and_select_group(threshold_group, "reltol=$(reltol); abstol=$(abstol)")
    outgroup = sim_node["outgroup"]
    if !dry
        data = Utilities.create_and_select_group(data, outgroup)
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect)
        flush(data)

        Hamiltonian = sys.Hamiltonian .+ diagm(sum([SpectralDensities.reorganization_energy(j) * bath.svecs[nb, :] .^ 2 for (nb, j) in enumerate(bath.Jw)]))
        sys_ops = [diagm(complex(bath.svecs[nb, :])) for nb = 1:size(bath.svecs, 1)]
        ρ0 = ParseInput.parse_operator(sim_node["rho0"], sys.Hamiltonian)

        if haskey(sim_node, "lindblad")
            @info "Using the PILD method."
            QDSimUtilities.print_citation(PILD_reference)
            decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
            L = [ParseInput.parse_operator(sim_node["lindblad"][i], Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
        else
            L = nothing
        end

        @time _, ρs = if isnothing(sys.external_fields)
            HEOM.propagate(; Hamiltonian, ρ0, sys_ops, Jw=bath.Jw, β=bath.β, num_modes, Lmax, dt=sim.dt, ntimes=sim.nsteps, threshold, L, extraargs=Utilities.DiffEqArgs(; reltol, abstol))
        else
            HEOM.propagate(; Hamiltonian, ρ0, sys_ops, Jw=bath.Jw, β=bath.β, num_modes, Lmax, dt=sim.dt, ntimes=sim.nsteps, threshold, L, external_fields=sys.external_fields, extraargs=Utilities.DiffEqArgs(; reltol, abstol))
        end
        Utilities.check_or_insert_value(data, "rho", ρs)
        flush(data)
    end
    data
end

function dynamics(::QDSimUtilities.Method"BlochRedfield", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a Bloch-Redfield calculation."
    end
    reltol = get(sim_node, "reltol", 1e-6)
    abstol = get(sim_node, "abstol", 1e-6)
    data = Utilities.create_and_select_group(dt_group, "reltol=$(reltol); abstol=$(abstol)")
    outgroup = sim_node["outgroup"]
    if !dry
        data = Utilities.create_and_select_group(data, outgroup)
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect)
        flush(data)

        Hamiltonian = sys.Hamiltonian .+ diagm(sum([SpectralDensities.reorganization_energy(j) * bath.svecs[nb, :] .^ 2 for (nb, j) in enumerate(bath.Jw)]))
        sys_ops = [diagm(complex(bath.svecs[nb, :])) for nb = 1:size(bath.svecs, 1)]
        ρ0 = ParseInput.parse_operator(sim_node["rho0"], sys.Hamiltonian)
        @time _, ρs = BlochRedfield.propagate(; Hamiltonian, ρ0, sys_ops, Jw=bath.Jw, β=bath.β, dt=sim.dt, ntimes=sim.nsteps, extraargs=Utilities.DiffEqArgs(; reltol, abstol))
        Utilities.check_or_insert_value(data, "rho", ρs)
        flush(data)
    end
    data
end

function dynamics(::QDSimUtilities.Method"Forster", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group}, sim_node; dry=false)
    if !dry
        @info "Running a Forster calculation. Please cite:"
    end
    data = Utilities.create_and_select_group(dt_group, "Forster")
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect)
        flush(data)

        Jw = Vector{SpectralDensities.SpectralDensity}(undef, length(bath.Jw))
        for (ind, jw) in enumerate(bath.Jw)
            target_ind = argmax(bath.svecs[ind, :])
            Jw[target_ind] = deepcopy(jw)
        end
        k, U = Forster.build_incoherent_propagator(; H=sys.Hamiltonian, Jw, dt=sim.dt, β=bath.β, verbose=true)
        Utilities.check_or_insert_value(data, "k", k)
        Utilities.check_or_insert_value(data, "U", U)
    end
    data
end

function dynamics(::QDSimUtilities.Method"Spin-LSC", units::QDSimUtilities.Units,
                  sys::QDSimUtilities.System, bath::QDSimUtilities.Bath,
                  sim::QDSimUtilities.Simulation, dt_group::Union{Nothing,HDF5.Group},
                  sim_node; dry=false)
    if !dry
        @info "Running a Spin-LSC calculation. Please cite:"
        QDSimUtilities.print_citation(SpinLSC.references)
    end

    transforms = Dict{String,Type{<:SpinLSC.SWTransform}}(
        "QTransform" => SpinLSC.QTransform,
        "WTransform" => SpinLSC.WTransform,
        "PTransform" => SpinLSC.PTransform)
    solvers = Dict{String,Type{<:SpinLSC.SpinLSCSolver}}(
        "RK4" => SpinLSC.RK4,
        "Verlet" => SpinLSC.Verlet,
    )

    transform = get(sim_node, "SW_transform", "QTransform")
    solver = get(sim_node, "solver", "RK4")
    transform_group = Utilities.create_and_select_group(dt_group, "SW_transform=$transform")
    solver_group = Utilities.create_and_select_group(transform_group, "solver=$solver")

    ρ0 = ParseInput.parse_operator(sim_node["rho0"], sys.Hamiltonian)

    nbins = get(sim_node, "num_bins", 1)
    nmc = sim_node["num_mc"]

    for n in 1:nbins
        Utilities.create_and_select_group(solver_group, "bin #$n")
    end

    outgroup = sim_node["outgroup"]

    if !dry
        @info "Running with $(Threads.nthreads()) threads."
        ρs = Vector{AbstractArray{<:Complex,3}}(undef, nbins)
        U0es = Vector{AbstractArray{<:Complex,3}}(undef, nbins)
        T0es = Vector{AbstractArray{<:Complex,3}}(undef, nbins)

        time = 0:sim.dt/units.time_unit:sim.nsteps*sim.dt/units.time_unit |> collect
        for n in 1:nbins
            data = Utilities.create_and_select_group(solver_group, "bin #$n")
            Utilities.check_or_insert_value(data, "num_mc", nmc)

            outgrouphdf5 = Utilities.create_and_select_group(data, outgroup)
            Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
            Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
            Utilities.check_or_insert_value(data, "time", time)
            Utilities.check_or_insert_value(outgrouphdf5, "time", time)
            Utilities.check_or_insert_value(outgrouphdf5, "time_unit", units.time_unit)
            flush(data)

            @info "Calculating bin $n of $nbins"
            U0e, ρ = SpinLSC.propagate(; Hamiltonian=sys.Hamiltonian, Jw=bath.Jw,
                                       β=bath.β, num_bath_modes=bath.num_osc,
                                       ρ0=ρ0, dt=sim.dt, ntimes=sim.nsteps,
                                       svec=bath.svecs,
                                       transform=transforms[transform],
                                       nmc=nmc, solver=solvers[solver],
                                       verbose=true, outgroup=outgroup,
                                       output=data)
            ρs[n] = ρ
            U0es[n] = U0e
            T0es[n] = read(data["T0e"])
        end

        U0e_mean, U0e_std = QDSimUtilities.matrix_avg_std(U0es)
        T0e_mean, T0e_std = QDSimUtilities.matrix_avg_std(T0es)
        ρ_mean, ρ_std = QDSimUtilities.matrix_avg_std(ρs)

        outgrouphdf5 = Utilities.create_and_select_group(solver_group, outgroup)
        Utilities.check_or_insert_value(solver_group, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(solver_group, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(solver_group, "time", time)
        Utilities.check_or_insert_value(outgrouphdf5, "time", time)
        Utilities.check_or_insert_value(outgrouphdf5, "time_unit", units.time_unit)

        Utilities.check_or_insert_value(solver_group, "U0e", U0e_mean)
        Utilities.check_or_insert_value(solver_group, "U0e_std", U0e_std)

        Utilities.check_or_insert_value(solver_group, "T0e", T0e_mean)
        Utilities.check_or_insert_value(solver_group, "T0e_std", T0e_std)

        Utilities.check_or_insert_value(outgrouphdf5, "rho", ρ_mean)
        Utilities.check_or_insert_value(outgrouphdf5, "rho_std", ρ_std)
        flush(outgrouphdf5)

        flush(solver_group)
    end
    solver_group
end

end
