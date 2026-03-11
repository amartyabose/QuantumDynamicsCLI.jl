"""
Various routines for post-processing and analysing the simulation results.
"""
@cast module Post

using Comonicon
using QuantumDynamics
using DelimitedFiles
using TOML
using LinearAlgebra
using Statistics: mean, stdm
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

function calculate_observable(sys::QDSimUtilities.System, ρs::AbstractArray{<:Complex,3}, obs::String;
                              mat_type::String="real")
    if obs == "trace"
        [tr(ρs[j, :, :]) for j in axes(ρs, 1)]
    elseif obs == "purity"
        [tr(ρs[j, :, :] * ρs[j, :, :]) for j in axes(ρs, 1)]
    elseif obs == "vonNeumann_entropy"
        [-tr(ρs[j, :, :] * log(ρs[j, :, :])) for j in axes(ρs, 1)]
    else
        op = ParseInput.parse_operator(obs, sys.Hamiltonian; mat_type)
        Utilities.expect(ρs, op)
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
        root_node = Simulate.calc(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node, method_group; dry=true)
        # This is a MC method.
        if !haskey(root_node, outputdir)
            nbins = read_dataset(root_node, "num_bins")
            ρs = [ read(root_node["bin #$b"][outputdir]["rho"]) for b in 1:nbins ]
            ts = read(root_node["bin #1"][outputdir]["time"])
        else
            nbins = 1
            ρs = [ read(root_node[outputdir]["rho"]) ]
            ts = read(root_node[outputdir]["time"])
        end
        dt = (ts[2] - ts[1]) * units.time_unit
        ωlim = π/dt
        dω = π/(ts[end] * units.time_unit)
        ω = -ωlim:dω:ωlim
        num_obs = length(sim_node["observable"])
        names = String[]
        ft = get(sim_node, "fourier_transform", false)
        full = true
        if ft
            full = get(sim_node, "full_transform", true)
        end
        vals = zeros(ComplexF64, length(ft ? ωs : ts), num_obs)
        vals_std = similar(vals)
        values = zeros(ComplexF64, length(ft ? ωs : ts), nbins)
        for (os, obs) in enumerate(sim_node["observable"])
            push!(names, obs["observable"])
            for b in 1:nbins
                values[:,b] .= calculate_observable(sys, ρs[b], obs["observable"];
                                                    mat_type=get(obs, "type", "real"))
                if ft
                    _, values[:,b] = Utilities.fourier_transform(ts, values[:, b]; full)
                end
            end
            if nbins > 1
                vals[:, os] .= mean(real.(values); dims=2)[:,1] .+
                         im .* mean(imag.(values); dims=2)[:,1]
                vals_std[:, os] .= (stdm(real.(values), real(vals[:, os]); dims=2)[:,1] .+
                              im .* stdm(imag.(values), imag(vals[:, os]); dims=2)[:,1]) / sqrt(nbins)
            else
                vals[:, os] .= values[:,1]
            end
        end

        obs_file = sim_node["observable_output"]
        fname, ext = splitext(obs_file)

        open("$(fname)_real$(ext)", "w") do io
            if ft
                write(io, "# (1)w")
            else
                write(io, "# (1)t")
            end
            for (j, n) in enumerate(names)
                write(io, "\t($(j+1))$(n)")
            end
            nbins > 1 && for (j,n) in enumerate(names)
                write(io, "\t($(j+1+num_obs))$n std")
            end
            write(io, "\n")
            if nbins > 1
                writedlm(io, hcat(round.(ft ? ω ./ units.energy_unit : ts; sigdigits=10),
                                  round.(real.(vals); sigdigits=10),
                                  round.(real.(vals_std); sigdigits=10)))
            else
                writedlm(io, hcat(round.(ft ? ω ./ units.energy_unit : ts; sigdigits=10),
                                  round.(real.(vals); sigdigits=10)))
            end
        end

        open("$(fname)_imag$(ext)", "w") do io
            if ft
                write(io, "# (1)w")
            else
                write(io, "# (1)t")
            end
            for (j, n) in enumerate(names)
                write(io, "\t($(j+1))$(n)")
            end
            nbins > 1 && for (j,n) in enumerate(names)
                write(io, "\t($(j+1+num_obs))$n std")
            end
            write(io, "\n")
            if nbins > 1
                writedlm(io, hcat(round.(ft ? ω ./ units.energy_unit : ts; sigdigits=10),
                                  round.(imag.(vals); sigdigits=10),
                                  round.(imag.(vals_std); sigdigits=10)))
            else
                writedlm(io, hcat(round.(ft ? ω ./ units.energy_unit : ts; sigdigits=10),
                                  round.(imag.(vals); sigdigits=10)))
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

function calculate_print_statetostate(::QDSimUtilities.Calculation"dynamics", sys::QDSimUtilities.System, bath::QDSimUtilities.Bath, sim::QDSimUtilities.Simulation, units::QDSimUtilities.Units, sim_node)
    out = h5open(sim.output, "r")
    method_group = out["$(sim.name)/$(sim.calculation)/$(sim.method)"]
    data_node = Simulate.calc(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node, method_group; dry=true)
    outputdir = sim_node["outgroup"]
    obs_file = sim_node["observable_output"]
    fname, ext = splitext(obs_file)

    ts = read(data_node[outputdir]["time"])
    ρs = read(data_node[outputdir]["rho"])
    if haskey(sim_node, "lindblad")
        decayconstant = [sim_node["decay_constant"][i] for i in 1:length(sim_node["decay_constant"])]
        L = [ParseInput.parse_operator(sim_node["lindblad"][i], sys.Hamiltonian) / sqrt(decayconstant[i] * units.time_unit) for i in 1:length(sim_node["decay_constant"])]
    else
        L = nothing
    end
    derivative = get(sim_node,"derivative", false)
    display(ts)
    display(sys.Hamiltonian / units.energy_unit)
    ddt_flows, flows = Utilities.statetostate(;t=(ts * units.time_unit), ρs=ρs, H0=sys.Hamiltonian, L=L)

    for i in axes(flows,1)
        header = []
        push!(header,"#t")
        for j in axes(flows,1)
            push!(header,"#$(j)")
        end

        open("flows_in_state_$(i)_$(fname)_real$(ext)", "w") do io
            writedlm(io, vcat(reshape(header,1,:), hcat(ts,real.(flows[i,:,:]))))
        end
        open("flows_in_state_$(i)_$(fname)_imag$(ext)", "w") do io
            writedlm(io, vcat(reshape(header,1,:), hcat(ts,imag.(flows[i,:,:]))))
        end

        if derivative
            open("derivative_flows_in_state_$(i)_$(fname)_real$(ext)", "w") do io
                writedlm(io, vcat(reshape(header,1,:), hcat(ts,real.(ddt_flows[i,:,:]))))
            end
            open("derivative_flows_in_state_$(i)_$(fname)_imag$(ext)", "w") do io
                writedlm(io, vcat(reshape(header,1,:), hcat(ts,imag.(ddt_flows[i,:,:]))))
            end
        end
    end
end

@cast function state_to_state(system_input, simulate_input)
    QDSimUtilities.print_banner()
    units, sys, bath = ParseInput.parse_system_bath(system_input)
    # sys_file = TOML.parsefile(system_input)
    # is_QuAPI = get(sys_file["system"], "is_QuAPI", true)
    # if !is_QuAPI
    #     sys.Hamiltonian .+= diagm(sum([SpectralDensities.reorganization_energy(j) * bath.svecs[nb, :] .^ 2 for (nb, j) in enumerate(bath.Jw)]))
    # end
    sim_file = TOML.parsefile(simulate_input)
    for (ns, sim_node) in enumerate(sim_file["simulation"])
        @info "Getting the state-to-state flows for simulation number $(ns). Please cite:"
        QDSimUtilities.print_citation(Utilities.statetostate_references)
        sim = ParseInput.parse_sim(sim_node, units)
        @assert isfile(sim.output) "File not present."
        calculate_print_statetostate(QDSimUtilities.Calculation(sim.calculation)(), sys, bath, sim, units, sim_node)
    end
end

end
