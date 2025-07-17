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

## Simulation File
The simulation file has only a single TOML section `[simulation]`. Every simulation file should have a `name` to identify the simulation and an `output` file that specifies an HDF5 file for storing all the data.

There are broadly three types of simulations that are currently supported:
- `dynamics` simulations that simulate the non-equilibrium dynamics of the given problem
- `equilibrium_rho` simulations that simulate the equilibrium density at the given temperature
- `complex_corr` simulations for calculating equilibrium correlation functions of various flavors
These are specified in the `calculation` field which, if unspecified, is taken to be `dynamics` by default.

The rest of the parameters required for a simulation file are specific to the type of simulation being run.

### Dynamics Simulations
Various methods of simulation are supported:
- Path Integral Methods using Feynman-Vernon Influence Functional[feynmanTheoryGeneralQuantum1963](@cite):
    - Quasi-adiabatic Propagator Path Integrals (QuAPI) [makriTensorPropagatorIterativeI1995, makriTensorPropagatorIterativeII1995](@cite)
    - Blip QuAPI [makriBlipDecompositionPath2014](@cite)
    - Time-Evolved Matrix Product Operators (TEMPO) [strathearnEfficientNonMarkovianQuantum2018](@cite)
    - Pairwise-Connected Tensor Network Path Integral (PC-TNPI) [bosePairwiseConnectedTensor2022](@cite)
- Hierarchical Equations of Motion (HEOM) [tanimuraNumericallyExactApproach2020](@cite)
- Generalized Quantum Master Equation
- Multichromophore Incoherest Forster Theory
- Bloch-Redfield Master Equation
- Transfer Tensor Method (TTM) [cerrilloNonMarkovianDynamicalMaps2014](@cite) coupled with any of the path integral methods

All of these dynamics methods require some core common parameters and then more specfic method-dependent parameters. The core parameters of all the dynamics methods are:
- `dt`: for the time-step in the units specified in the system file
- `nsteps`: for the number of steps of simulation of the dynamics

#### Feynman-Vernon Influence Functional Simulations
There are two ways of incorporating the effect of non-Markovian memory in path integral simulations: iterative propagation beyond memory or using the transfer tensor method. To use TTM, one can choose one of the following:
- `method = "QuAPI-TTM"`: for using QuAPI within memory
- `method = "Blip-TTM"`: for using Blip within memory
- `method = "TEMPO-TTM"`: for using TEMPO within memory

Finally, if one does not intend on using TTM, we suggest using TEMPO for accessing long memory lengths efficiently. This is chosen by setting `method = "TEMPO"`.

For all the TTM-based methods, the memory length in terms of number of time-steps is set through the parameter, `rmax`. For `method = "TEMPO"`, the memory length is set through the parameter, `kmax`.