module ParseInput

using DelimitedFiles
using LinearAlgebra
using TOML
using Unitful, UnitfulAtomic

using QuantumDynamics
using ..QDSimUtilities

function read_matrix(fname::String, mat_type::String="real")
    if mat_type == "real"
        Matrix{ComplexF64}(readdlm(fname))
    elseif mat_type == "complex"
        readdlm(fname, ComplexF64)
    end
end

"""
    parse_unit(input_dict)
Parses the `[units]` section of the `system` TOML file. It takes two variables:

Parameters obtained:
- `energy_unit_name` [Default: `ha`]: Specifies the energy units. Typically `eV`, `meV`, `cm^-1`, or `ha`.
- `time_unit_name` [Default: `au`]: Specifies the time units. Typically `fs`, or `au`.
"""
function parse_unit(input_dict)
    energy_unit = 1.0
    energy_unit_name = "ha"
    time_unit = 1.0
    time_unit_name = "au"

    if haskey(input_dict, "units")
        energy_unit_name = get(input_dict["units"], "energy_unit_name", "ha")
        time_unit_name = get(input_dict["units"], "time_unit_name", "au")

        energy_unit_val = uparse(energy_unit_name)
        if energy_unit_name == "cm^-1"
            energy_unit_val *= Unitful.c * Unitful.h
        end
        time_unit_val = time_unit_name != "au" ? uparse(time_unit_name) : Unitful.hbar * 1u"ha^-1"

        energy_unit = austrip(1 * energy_unit_val)
        time_unit = austrip(1 * time_unit_val)
    end

    QDSimUtilities.Units(energy_unit, energy_unit_name, time_unit, time_unit_name)
end

"""
    parse_system(sys_inp, unit)
Parses the `[system]` portion of the system file for the Hamiltonian and converts it to atomic units for internal use.

Parameters obtained:
- `type` [Default: `real`, Optionally: `complex`]: the type of elements in the Hamiltonian.
- `Htype` [Default: `file`]: How to parse the Hamiltonian file. Typically `file` to read a file, `nearest_neighbor` to specify the Hamiltonian as a nearest-neighbor tight-binding model, or `nearest_neighbor_cavity` for a Hamiltonian consisting of a nearest-neighbor tight-binding part and a cavity with interacts with all the sites.

Then the Hamiltonian needs to be specified in the energy units that are being used. This can be done in several ways depending on the value of `Htype`:
- if `Htype` = `file`, put in the variable, `Hamiltonian` with the name of a file containing the elements of the Hamiltonian matrix.
- if `Htype = "nearest_neighbor"`, specify the following variables:
    - `site_energy`: the energy of each site
    - `coupling`: the intersite coupling element
    - `num_sites`: the number of sites
- if `Htype = "nearest_neighbor_cavity"`, in addition to the same variables as the `Htype = "nearest_neighbor"` case, also specify:
    - `cavity_energy`: the energy of the cavity mode
    - `cavity_coupling`: coupling of the cavity to the monomers

The parameter `is_QuAPI` [Default: `true`] specifies if the specified system Hamiltonian should be interpreted to be a part of the QuAPI system-bath Hamiltonian form or not.

If external time-dependent fields have to be incorporated, the subsection `[external_field]` should be used. By default no external fields are involved. If used, the subsection would have two variables:
- `coupling_op`: the operator through which the external time-dependent field interacts with the system
- `V_t`: this is a single line expression which evaluates to the value of the field as a function of time, t. This variable needs to be called `t`.
"""
function parse_system(sys_inp, unit)
    Htype = get(sys_inp, "Htype", "file")
    H0 = if Htype == "file"
        read_matrix(sys_inp["Hamiltonian"], get(sys_inp, "type", "real")) * unit.energy_unit
    elseif Htype == "nearest_neighbor"
        site_energy = sys_inp["site_energy"] * unit.energy_unit
        coupling = sys_inp["coupling"] * unit.energy_unit
        num_sites = sys_inp["num_sites"]
        Utilities.create_nn_hamiltonian(; site_energies=repeat([site_energy], num_sites), couplings=repeat([coupling], num_sites-1), periodic=false)
    elseif Htype == "nearest_neighbor_cavity"
        site_energy = sys_inp["site_energy"] * unit.energy_unit
        coupling = sys_inp["coupling"] * unit.energy_unit
        num_sites = sys_inp["num_sites"]
        cavity_energy = sys_inp["cavity_energy"] * unit.energy_unit
        cavity_coupling = sys_inp["cavity_coupling"] * unit.energy_unit
        H = Matrix{ComplexF64}(undef, num_sites+1, num_sites+1)
        H[1:num_sites, 1:num_sites] .= Utilities.create_nn_hamiltonian(; site_energies=repeat([site_energy], num_sites), couplings=repeat([coupling], num_sites-1), periodic=false)
        H[end, :] .= cavity_coupling / sqrt(num_sites)
        H[:, end] .= cavity_coupling / sqrt(num_sites)
        H[end, end] = cavity_energy
        H
    end

    external_fields = if haskey(sys_inp, "external_field")
        EFs = Vector{Utilities.ExternalField}()
        for EF in sys_inp["external_field"]
            coupling_op = read_matrix(EF["coupling_op"], "real") .+ 0.0im
            Vstr = "t -> " * EF["V_t"]
            V(t::Float64) = Base.invokelatest(eval(Meta.parse(Vstr)), t)::Float64
            push!(EFs, Utilities.ExternalField(V, coupling_op))
        end
        EFs
    else
        nothing
    end

    QDSimUtilities.System(Htype, H0, external_fields)
end

"""
    get_bath(b, unit)
Parse individual baths

Parameters:
- `type` [Default: `ohmic`]: Type of harmonic bath. Can be `ohmic`, `drude_lorentz`, `tabular`, `tabular_jw_over_w`, `huang_rhys`

The parameters for the bath are specified in different ways for each different bath:
- if `type="ohmic"`, the bath spectral density ``J(\\omega)=\\frac{2\\pi}{\\Delta s^2}\\hbar\\xi\\omega_c\\exp\\left(-\\frac{\\omega}{\\omega_c}\\right)``
    - `xi`: dimensionless Kondo parameter (``\\xi``)
    - `omegac`: cutoff frequency (``\\omega_c``) in energy units
    - `Ds` [Default: 2]: ``\\Delta s`` value. Typically set to 1 for exciton transfer problems.
    - `npoints` [Default: 100000]: Number of points used in integration for the influence functional coefficients.
    - `omega_max` [Default: `30.0 * omegac`]: Upper limit of influence functional coefficient integrations
- if `type="drude_lorentz"`, the bath spectral density ``J(\\omega) = \\frac{2\\lambda}{\\Delta s^2}\\frac{\\omega \\gamma}{\\omega^2+\\gamma^2}``
    - `lambda`: reorganization energy, ``\\lambda``, in energy units
    - `gamma`: cutoff frequency, ``\\gamma``, in energy units
    - `Ds` [Default: 2]: ``\\Delta s`` value. Typically set to 1 for exciton transfer problems.
    - `npoints` [Default: 100000]: Number of points used in integration for the influence functional coefficients.
    - `omega_max` [Default: `1000.0 * gamma`]: Upper limit of influence functional coefficient integrations
- if `type="tabular"`, the bath spectral density is provided as a table in the file specified in `jw_file`
- if `type="huang_rhys"`, the bath spectral density is provided in the form of a table of frequency-dependent Huang-Rhys factors in the file specified in `huang_rhys_file`
"""
function get_bath(b, unit)
    sd_type = get(b, "type", "ohmic")
    J = nothing
    if sd_type == "ohmic"
        ξ = b["xi"]
        ωc = b["omegac"]
        n = get(b, "n", 1.0)
        Δs = get(b, "Ds", 2.0)
        npoints = get(b, "npoints", 100000)
        ωmax = get(b, "omega_max", 30.0 * ωc)
        ωc *= unit.energy_unit
        ωmax *= unit.energy_unit
        classical = get(b, "classical", false)
        J = SpectralDensities.ExponentialCutoff(; ξ, ωc, n, Δs, ωmax, npoints, classical)
    elseif sd_type == "drude_lorentz"
        λ = b["lambda"] * unit.energy_unit
        γ = b["gamma"]
        ωmax = get(b, "omega_max", 1000.0 * γ)
        γ *= unit.energy_unit
        ωmax *= unit.energy_unit
        npoints = get(b, "npoints", 100000)
        Δs = get(b, "Ds", 2.0)
        classical = get(b, "classical", false)
        J = SpectralDensities.DrudeLorentz(; λ, γ, Δs, ωmax, npoints, classical)
    elseif sd_type == "tabular"
        inpfile = b["jw_file"]
        skipstart = get(b, "skipstart", 1)
        classical = get(b, "classical", false)
        J = SpectralDensities.read_jw(inpfile; skipstart, classical)
    elseif sd_type == "tabular_jw_over_w"
        inpfile = b["jw_over_w_file"]
        skipstart = get(b, "skipstart", 1)
        classical = get(b, "classical", false)
        J = SpectralDensities.read_jw(inpfile; skipstart, classical)
    elseif sd_type == "huang_rhys"
        inpfile = b["huang_rhys_file"]
        skipstart = get(b, "skipstart", 1)
        classical = get(b, "classical", false)
        J = SpectralDensities.read_huang_rhys(inpfile; skipstart, classical)
    else
        throw(ArgumentError("Spectral density of type $(sd_type) not supported."))
    end
    J
end

"""
    parse_bath(baths, sys, unit)
This function parses the bath(s) interacting with the system. Currently only harmonic baths are supported.

Parameters:
- The temperature for the simulation needs to be specfied. This can be done in one of two ways:
    - `beta` in units of 1/ha for the inverse temperature ``\\beta = \\frac{1}{k_BT}``.
    - or `temperature` in the units of Kelvin.
- A list of baths that interact with the system specified under `[[baths.bath]]` heading in the TOML file. Each bath is parsed by [`QuantumDynamicsCLI.ParseInput.get_bath`](@ref).
"""
function parse_bath(baths, sys, unit)
    β = 0.0
    if haskey(baths, "beta")
        β = baths["beta"]
    elseif haskey(baths, "temperature")
        β = 1.0 / austrip(baths["temperature"] * uparse("K") * Unitful.k)
    end
    Jw = Vector{SpectralDensities.SpectralDensity}()
    svecs = zeros(size(sys.Hamiltonian, 1), size(sys.Hamiltonian, 1))
    btype = get(baths, "baths_type", "list")
    num_osc = Vector{Integer}()
    if btype == "list"
        svecs = zeros(length(baths["bath"]), size(sys.Hamiltonian, 1))
        for (nb, b) in enumerate(baths["bath"])
            push!(Jw, get_bath(b, unit))
            svecs[nb, :] .= b["svec"]
            haskey(b, "num_osc") && push!(num_osc, b["num_osc"])
        end
    elseif btype == "site_based"
        bath = get_bath(baths["bath"][1], unit)
        nsites = sys.Htype == "nearest_neighbor_cavity" ? size(sys.Hamiltonian, 1) - 1 : size(sys.Hamiltonian, 1)
        svecs = zeros(nsites, size(sys.Hamiltonian, 1))
        for nb = 1:nsites
            push!(Jw, bath)
            svecs[nb, nb] = 1.0
        end
        haskey(b, "num_osc") && (num_osc = fill(b["num_osc"], nsites))
    end
    QDSimUtilities.Bath(β, Jw,svecs,(isempty(num_osc) ? nothing : num_osc))
end

function parse_system_bath(input_file)
    input_dict = TOML.parsefile(input_file)
    unit = parse_unit(input_dict)
    sys = parse_system(input_dict["system"], unit)
    bath = parse_bath(input_dict["baths"], sys, unit)
    is_QuAPI = get(input_dict["system"], "is_QuAPI", true)
    if !is_QuAPI
        sys.Hamiltonian .-= diagm(sum([SpectralDensities.reorganization_energy(j) * bath.svecs[nb, :] .^ 2 for (nb, j) in enumerate(bath.Jw)]))
    end
    unit, sys, bath
end

function parse_sim(sim, unit)
    name = sim["name"]
    calculation = get(sim, "calculation", "dynamics")
    method = sim["method"]
    output = sim["output"]
    nsteps = sim["nsteps"]
    dt = get(sim, "dt", 0.0)
    QDSimUtilities.Simulation(name, calculation, method, output, dt, nsteps)
end

function parse_operator(op, Hamiltonian)
    obs = zero(Hamiltonian)
    if startswith(op, "P_")
        state = parse(Int64, split(op, "_")[2])
        obs[state, state] = 1
        obs
    elseif startswith(op, "F_")
        state = parse(Int64, split(op, "_")[2])
        obs[state, state] = 1
        1im * Utilities.commutator(Hamiltonian, obs)
    elseif startswith(op, "O_")
        states = parse.(Int64, split(split(op, "_")[2], ","))
        obs[states[1], states[2]] = 1
        obs
    elseif op == "id"
        for j = 1:size(Hamiltonian, 1)
            obs[j, j] = 1.0
        end
        obs
    else
        ParseInput.read_matrix(op)
    end
end

end
