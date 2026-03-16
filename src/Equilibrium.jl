module Equilibrium

using QuantumDynamics
using LinearAlgebra
using DelimitedFiles
using ..QDSimUtilities, ..ParseInput

function rho(::QDSimUtilities.Method"TNPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dat_group::HDF5.Group, sim_node; dry=false)
    if !dry
        @info "Running an imaginary time TNPI simulation. Please cite:"
        QDSimUtilities.print_citation(ComplexTNPI.references)
    end
    cutoff = get(sim_node, "cutoff", 1e-10)
    maxdim = get(sim_node, "maxdim", 1000)
    algorithm = get(sim_node, "algorithm", "naive")
    @info "Running with $(BLAS.get_num_threads()) threads."

    maxdim_group = Utilities.create_and_select_group(dat_group, "maxdim=$(maxdim)")
    cutoff_group = Utilities.create_and_select_group(maxdim_group, "cutoff=$(cutoff)")
    data = Utilities.create_and_select_group(cutoff_group, "algorithm=$(algorithm)")
    if !dry
        flush(data)
        extraargs = Utilities.TensorNetworkArgs(; cutoff, maxdim, algorithm)

        nsites = size(sys.Hamiltonian, 1)

        idmat = Matrix(1.0I, nsites, nsites)
        At, _ = ComplexTNPI.A_of_t(; Hamiltonian=sys.Hamiltonian, β=bath.β, t=0.0, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A=idmat, extraargs)
        Z = tr(At)
        Utilities.check_or_insert_value(data, "eqm_rho", real.(At / Z))
    end
    data
end

function complex_time_correlation_function(::QDSimUtilities.Method"TNPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dat_group::HDF5.Group, sim_node; dry=false)
    if !dry
        @info "Running a complex time TNPI simulation. Please cite:"
        QDSimUtilities.print_citation(ComplexTNPI.references)
    end
    cutoff = get(sim_node, "cutoff", 1e-10)
    maxdim = get(sim_node, "maxdim", 1000)
    algorithm = get(sim_node, "algorithm", "naive")
    tfinal = sim_node["tfinal"] * units.time_unit
    @info "Running with $(BLAS.get_num_threads()) threads."

    type_corr = get(sim_node, "corr_type", "symm")

    maxdim_group = Utilities.create_and_select_group(dat_group, "maxdim=$(maxdim)")
    cutoff_group = Utilities.create_and_select_group(maxdim_group, "cutoff=$(cutoff)")
    data = Utilities.create_and_select_group(cutoff_group, "algorithm=$(algorithm)")
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:tfinal/units.time_unit |> collect)
        flush(data)
        extraargs = Utilities.TensorNetworkArgs(; cutoff, maxdim, algorithm)

        nsites = size(sys.Hamiltonian, 1)

        idmat = Matrix(1.0I, nsites, nsites)
        @info "Calculating the partition function for normalization."
        At, _ = ComplexTNPI.A_of_t(; Hamiltonian=sys.Hamiltonian, β=bath.β, t=0.0, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A=idmat, extraargs)
        norm_op = get(sim_node, "partition_function", "id")
        Z = real(tr(At * ParseInput.parse_operator(norm_op, sys.Hamiltonian)))
        @info "Partition function = $(Z)."
        @info "Saving the equilibrium density matrix."
        Utilities.check_or_insert_value(data, "eqm_rho", real.(At / tr(At)))
        A = ParseInput.parse_operator(sim_node["A"], sys.Hamiltonian)
        B = ParseInput.parse_operator(sim_node["B"], sys.Hamiltonian)
        ts, corr, _ = ComplexTNPI.complex_correlation_function(; Hamiltonian=sys.Hamiltonian, β=bath.β, tfinal, dt=sim.dt, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A, B=[B], Z, verbose=true, extraargs, output=data, type_corr)
        
        ft = get(sim_node, "fourier_transform", false)
        if ft
            sig = vcat(reverse(conj.(corr[2:end])), corr)
            ts = vcat(reverse(-ts[2:end]), ts)
            conjugated = get(sim_node, "conjugate", false)
            ωs, spectrum = conjugated ? Utilities.fourier_transform(ts, sig; neg=true) : Utilities.fourier_transform(ts, sig; neg=false)
            Utilities.check_or_insert_value(data, "frequency", ωs ./ units.energy_unit)
            Utilities.check_or_insert_value(data, "spectrum", real.(spectrum))
        end
    end
    data
end

function rho(::QDSimUtilities.Method"QuAPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dat_group::HDF5.Group, sim_node; dry=false)
    if !dry
        @info "Running an imaginary time TNPI simulation. Please cite:"
        QDSimUtilities.print_citation(ComplexTNPI.references)
    end
    cutoff = get(sim_node, "cutoff", 1e-10)
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end

    data = Utilities.create_and_select_group(dat_group, "cutoff=$(cutoff)")
    if !dry
        flush(data)
        extraargs = QuAPI.QuAPIArgs(; cutoff)

        nsites = size(sys.Hamiltonian, 1)

        idmat = Matrix(1.0I, nsites, nsites)
        At, _ = ComplexQuAPI.A_of_t(; Hamiltonian=sys.Hamiltonian, β=bath.β, t=0.0, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A=idmat, extraargs, exec=QDSimUtilities.parse_exec(exec))
        Z = tr(At)
        Utilities.check_or_insert_value(data, "eqm_rho", real.(At / Z))
    end
    data
end

function complex_time_correlation_function(::QDSimUtilities.Method"QuAPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dat_group::HDF5.Group, sim_node; dry=false)
    if !dry
        @info "Running a complex time QuAPI simulation. Please cite:"
        QDSimUtilities.print_citation(ComplexTNPI.references)
    end
    cutoff = get(sim_node, "cutoff", 1e-10)
    tfinal = sim_node["tfinal"] * units.time_unit
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end

    type_corr = get(sim_node, "corr_type", "symm")

    data = Utilities.create_and_select_group(dat_group, "cutoff=$(cutoff)")
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:tfinal/units.time_unit |> collect)
        flush(data)
        extraargs = QuAPI.QuAPIArgs(; cutoff)

        nsites = size(sys.Hamiltonian, 1)

        idmat = Matrix(1.0I, nsites, nsites)
        @info "Calculating the partition function for normalization."
        At, _ = ComplexQuAPI.A_of_t(; Hamiltonian=sys.Hamiltonian, β=bath.β, t=0.0, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A=idmat, extraargs)
        norm_op = get(sim_node, "partition_function", "id")
        Z = real(tr(At * ParseInput.parse_operator(norm_op, sys.Hamiltonian)))
        @info "Partition function = $(Z)."
        @info "Saving the equilibrium density matrix."
        Utilities.check_or_insert_value(data, "eqm_rho", real.(At / tr(At)))
        A = ParseInput.parse_operator(sim_node["A"], sys.Hamiltonian)
        B = ParseInput.parse_operator(sim_node["B"], sys.Hamiltonian)
        ts, corr, _ = ComplexQuAPI.complex_correlation_function(; Hamiltonian=sys.Hamiltonian, β=bath.β, tfinal, dt=sim.dt, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A, B=[B], Z, verbose=true, extraargs, output=data, type_corr)
        ft = get(sim_node, "fourier_transform", false)
        if ft
            sig = vcat(reverse(conj.(corr[2:end])), corr)
            ts = vcat(reverse(-ts[2:end]), ts)
            conjugated = get(sim_node, "conjugate", false)
            ωs, spectrum = conjugated ? Utilities.fourier_transform(ts, sig; neg=true) : Utilities.fourier_transform(ts, sig; neg=false)
            Utilities.check_or_insert_value(data, "frequency", ωs ./ units.energy_unit)
            Utilities.check_or_insert_value(data, "spectrum", real.(spectrum))
        end
    end
    data
end

function complex_time_correlation_function(::QDSimUtilities.Method"adaptive-kinks-QuAPI", units::QDSimUtilities.Units, sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, dat_group::HDF5.Group, sim_node; dry=false)
    if !dry
        @info "Running a complex time QuAPI simulation. Please cite:"
        QDSimUtilities.print_citation(ComplexTNPI.references)
    end
    cutoff = get(sim_node, "cutoff", 1e-10)
    tfinal = sim_node["tfinal"] * units.time_unit
    exec = get(sim_node, "exec", "ThreadedEx")
    if exec != "SequentialEx"
        @info "Running with $(Threads.nthreads()) threads."
    end

    type_corr = get(sim_node, "corr_type", "symm")

    cutoff_group = Utilities.create_and_select_group(dat_group, "cutoff=$(cutoff)")
    prop_cutoff = get(sim_node, "propagator_cutoff", 0.0)
    data = Utilities.create_and_select_group(cutoff_group, "prop_cutoff=$(prop_cutoff)")
    if !dry
        Utilities.check_or_insert_value(data, "dt", sim.dt / units.time_unit)
        Utilities.check_or_insert_value(data, "time_unit", units.time_unit)
        Utilities.check_or_insert_value(data, "time", 0:sim.dt/units.time_unit:tfinal/units.time_unit |> collect)
        flush(data)
        extraargs = QuAPI.QuAPIArgs(; cutoff, prop_cutoff)

        nsites = size(sys.Hamiltonian, 1)

        idmat = Matrix(1.0I, nsites, nsites)
        @info "Calculating the partition function for normalization."
        At, _ = ComplexQuAPI.adaptive_kink_A_of_t(; Hamiltonian=sys.Hamiltonian, β=bath.β, t=0.0, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A=idmat, extraargs)
        norm_op = get(sim_node, "partition_function", "id")
        Z = real(tr(At * ParseInput.parse_operator(norm_op, sys.Hamiltonian)))
        @info "Partition function = $(Z)."
        @info "Saving the equilibrium density matrix."
        Utilities.check_or_insert_value(data, "eqm_rho", real.(At / tr(At)))
        A = ParseInput.parse_operator(sim_node["A"], sys.Hamiltonian)
        B = ParseInput.parse_operator(sim_node["B"], sys.Hamiltonian)
        ts, corr, _ = ComplexQuAPI.adaptive_kink_complex_correlation_function(; Hamiltonian=sys.Hamiltonian, β=bath.β, tfinal, dt=sim.dt, N=sim.nsteps, Jw=bath.Jw, svec=bath.svecs, A, B=[B], Z, verbose=true, extraargs, output=data, type_corr)
        ft = get(sim_node, "fourier_transform", false)
        if ft
            sig = vcat(reverse(conj.(corr[2:end])), corr)
            ts = vcat(reverse(-ts[2:end]), ts)
            conjugated = get(sim_node, "conjugate", false)
            ωs, spectrum = conjugated ? Utilities.fourier_transform(ts, sig; neg=true) : Utilities.fourier_transform(ts, sig; neg=false)
            Utilities.check_or_insert_value(data, "frequency", ωs ./ units.energy_unit)
            Utilities.check_or_insert_value(data, "spectrum", real.(spectrum))
        end
    end
    data
end

end
