# Inputs to qdsim

A simulation depends upon the specification of the system. This is done in the `system` TOML file, which specifies three different things:
- the `units` in use
- the Hamiltonian for the problem being simulated:
    - the description of the `system`
    - the description of the `bath` or environment

Along with this `system` TOML file, a `simulation` TOML file also needs to be prepared that provides the details of the simulation.

## System File
Each of the three sections in the `system` file has a dedicated function for parsing the data. Before giving a full example of a `system` input file, we discuss the parameters accepted by the various sections.

### Specifying the Units
```@docs
QuantumDynamicsCLI.ParseInput.parse_unit
```

### Specifying the Hamiltonian
The basic problem under study has a general form, $$\hat{H} = \hat{H}_0 + \hat{H}_\text{env}$$, where $\hat{H}_0$ forms the system Hamiltonian and the $\hat{H}_\text{env}$ is the environment or bath Hamiltonian.

#### System Hamiltonian
```@docs
QuantumDynamicsCLI.ParseInput.parse_system
```

#### Bath Hamiltonian
```@docs
QuantumDynamicsCLI.ParseInput.parse_bath
```

```@docs
QuantumDynamicsCLI.ParseInput.get_bath
```

### Specifying Operators

The specification of operators need to be specified in multiple places. Convenient shorthands for the specification of most commonly used operators. These can be used to specify the initial reduced density matrix as well.

```@docs
QuantumDynamicsCLI.ParseInput.parse_operator
```

## Simulation File
The simulation file has only a single TOML section `[simulation]`. Every simulation file should have a `name` to identify the simulation and an `output` file that specifies an HDF5 file for storing all the data.

There are broadly three types of simulations that are currently supported:
- `dynamics` simulations that simulate the non-equilibrium dynamics of the given problem
- `equilibrium_rho` simulations that simulate the equilibrium density at the given temperature
- `complex_corr` simulations for calculating equilibrium correlation functions of various flavors
These are specified in the `calculation` field which, if unspecified, is taken to be `dynamics` by default.

The rest of the parameters required for a simulation file are specific to the type of simulation being run.
