"""
Various routines for post-processing and analysing the simulation results.
"""
@cast module Post

using Comonicon
using QuantumDynamics
using DelimitedFiles
using TOML
using ..ParseInput, ..Simulate, ..QDSimUtilities

"""
Combine the source files `sources` into `output`. If `output` does not exist, it is created.

# Args
- `sources`: source output files
- `output`: destination output file
"""
@cast function merge_into(sources::String...; output::String)
    for input in sources
        Utilities.merge_into(input, output)
    end
end

function calculate_print_observable(::QDSimUtilities.Calculation"dynamics", sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, units::QDSimUtilities.Units, sim_node)
    out = h5open(sim.output, "r")
    method_group = out["$(sim.name)/$(sim.calculation)/$(sim.method)"]
    if sim.method == "Forster"
        data_node = Simulate.calc(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node, method_group; dry=true)
        ts = read_dataset(data_node, "time")
        U = read_dataset(data_node, "U")
        ρ = zeros(Float64, size(U, 1))
        for j = 1:size(U, 1)
            ρ .= 0
            ρ[j] = 1.0
            open("populations_init$(j).dat", "w") do io
                write(io, "# (1)t ")
                for k = 1:size(U, 1)
                    write(io, "($(k+1))P_$(k) ")
                end
                write(io, "\n")
                for t in ts
                    write(io, "$(round(t; sigdigits=10)) ")
                    for k = 1:size(U, 1)
                        write(io, "$(round(ρ[k]; sigdigits=10)) ")
                    end
                    write(io, "\n")
                    ρ = U * ρ
                end
            end
        end
    else
        outputdir = sim_node["outgroup"]
        data_node = Simulate.calc(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node, method_group; dry=true)[outputdir]
        ts = read_dataset(data_node, "time")
        dt = (ts[2] - ts[1]) * units.time_unit
        ωlim = π/dt
        dω = π/(ts[end] * units.time_unit)
        ω = -ωlim:dω:ωlim
        ρs = read_dataset(data_node, "rho")
        num_obs = length(sim_node["observable"])
        names = String[]
        ft = get(sim_node, "fourier_transform", false)
        full = true
        if ft
            full = get(sim_node, "full_transform", true)
        end
        vals = ft ? zeros(ComplexF64, length(ω), num_obs) : zeros(ComplexF64, length(ts), num_obs)
        values = ft ? zeros(ComplexF64, length(ω)) : zeros(ComplexF64, length(ts))
        for (os, obs) in enumerate(sim_node["observable"])
            push!(names, obs["observable"])
            if obs["observable"] == "trace"
                values .= [tr(ρs[j, :, :]) for j in axes(ρs, 1)]
            elseif obs["observable"] == "purity"
                values .= [tr(ρs[j, :, :] * ρs[j, :, :]) for j in axes(ρs, 1)]
            elseif obs["observable"] == "vonNeumann_entropy"
                values = [-tr(ρs[j, :, :] * log(ρs[j, :, :])) for j in axes(ρs, 1)]
            else
                obs = ParseInput.parse_operator(obs["observable"], sys.Hamiltonian;
                                                mat_type=get(obs, "type", "real"))
                values = Utilities.expect(ρs, obs)
            end
            _, valft = ft ? Utilities.fourier_transform(ts, values; full=full) : (ts, values)
            vals[:, os] .= ft ? valft : valft
        end

        obs_file = sim_node["observable_output"]
        fname, ext = splitext(obs_file)

        open("$(fname)_real$(ext)", "w") do io
            if ft
                write(io, "# (1)w ")
            else
                write(io, "# (1)t ")
            end
            for (j, n) in enumerate(names)
                write(io, "($(j+1))$(n) ")
            end
            write(io, "\n")
            if ft
                writedlm(io, [round.(ω ./ units.energy_unit; sigdigits=10) round.(real.(vals); sigdigits=10)])
            else
                writedlm(io, [round.(ts; sigdigits=10) round.(real.(vals); sigdigits=10)])
            end
        end

        open("$(fname)_imag$(ext)", "w") do io
            if ft
                write(io, "# (1)w ")
            else
                write(io, "# (1)t ")
            end
            for (j, n) in enumerate(names)
                write(io, "($(j+1))$(n) ")
            end
            write(io, "\n")
            if ft
                writedlm(io, [round.(ω ./ units.energy_unit; sigdigits=10) round.(imag.(vals); sigdigits=10)])
            else
                writedlm(io, [round.(ts; sigdigits=10) round.(imag.(vals); sigdigits=10)])
            end
        end
    end
end

function calculate_print_observable(::QDSimUtilities.Calculation"complex_corr", sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, units::QDSimUtilities.Units, sim_node)
    out = h5open(sim.output, "r")
    method_group = out["$(sim.name)/$(sim.calculation)/$(sim.method)"]
    data_node = Simulate.calc(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node, method_group; dry=true)
    obs_file = sim_node["observable_output"]
    fname, ext = splitext(obs_file)
    ts = read_dataset(data_node, "time")
    corr = read_dataset(data_node, "corr")

    open("$(fname)$(ext)", "w") do io
        writedlm(io, [round.(ts; sigdigits=10) round.(real.(corr); sigdigits=10) round.(imag.(corr); sigdigits=10)])
    end

    ft = get(sim_node, "fourier_transform", false)
    if ft
        freq = read_dataset(data_node, "frequency")
        spect = read_dataset(data_node, "spectrum")
        open("$(fname)_spectrum$(ext)", "w") do io
            writedlm(io, [round.(freq; sigdigits=10) round.(real.(spect); sigdigits=10) round.(imag.(spect); sigdigits=10)])
        end
    end
end

@cast function get_observable(system_input, simulate_input)
    QDSimUtilities.print_banner()
    units, sys, bath = ParseInput.parse_system_bath(system_input)
    sim_file = TOML.parsefile(simulate_input)
    for (ns, sim_node) in enumerate(sim_file["simulation"])
        @info "Getting observables for simulation number $(ns)."
        sim = ParseInput.parse_sim(sim_node, units)
        @assert isfile(sim.output) "File not present."
        calculate_print_observable(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node)
    end
end

end
