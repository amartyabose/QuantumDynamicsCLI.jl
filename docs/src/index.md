# QuantumDynamicsCLI

| **Documentation** |
|:-----------------:|
|[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amartyabose.github.io/QuantumDynamicsCLI.jl/dev/)|

## What is QuantumDynamicsCLI?
Simulating the dynamics of quantum systems is a challenging task with a
multitude of complicated computational methods. The
[QuantumDynamics.jl](https://github.com/amartyabose/QuantumDynamics.jl) package
provides modular open-source implementations of an increasingly growing number
of these methods, while remaining a flexible platform for further development.
However, owing primarily to its exceptionally flexible nature, the usage of
QuantumDynamics.jl happens through short Julia scripts. This means that for the
most common simulation jobs, one needs to effectively rewrite the same code
multiple times increasing the chances of errors. As a means to making some of
the common types of simulations more facile, we now offer the
QuantumDynamicsCLI.jl package which installs the qdsim application as a sister
code of the QuantumDynamics.jl package. As the framework grows, so will this
application grow to accommodate the new methods and their most common use cases.

## Installation
QuantumDynamicsCLI.jl is a registered package. Installation is a simple procedure. It can be done either through the Pkg REPL:
```bash
~ julia
```

```
julia> ]
pkg> add QuantumDynamicsCLI
```

or by using the `Pkg` package manager in a script as follows:
```julia
julia> using Pkg
julia> Pkg.add("QuantumDynamicsCLI")
```

After the package gets built, an executable called `qdsim` will be placed in `$HOME/.julia/bin` along with the code completions for the shell in `$HOME/.julia/completions`. Please add `$HOME/.julia/bin` to your path and source the correct completions file.

While QuantumDynamicsCLI.jl builds on top of the QuantumDynamics.jl package, separate installation of that package is unnecessary. Just installing QuantumDynamicsCLI.jl would install QuantumDynamics.jl as a dependency.

## Basic Usage
`qdsim` comes as a single program with multiple sub-components. These components can call each other, but are mostly meant for the end-user, and are used for running simulations and post-processing the data. The general syntax for running any particular component is as follows:
```bash
> qdsim <component_name> <command_name> <arguments>
```

Currently, two `command`s are supported:
- `simulate`: gives access to various techniques for simulating the dynamics
- `post`: provides post-processing tools for the output

The most important sub-command of `simulate` is `run`, and for `post` is `get-observable`.

## Types of Simulations
The primary focus of the qdsim application provided by the QuantumDynamicsCLI.jl package and the underlying QuantumDynamics.jl package is the simulation of dynamics and spectra of open quantum systems. New methods are consistently added and the support for old methods improved. Currently the following methods are supported:
- Path Integral Methods using Feynman-Vernon Influence Functional[feynmanTheoryGeneralQuantum1963](@cite):
    - Quasi-adiabatic Propagator Path Integrals (QuAPI) [makriTensorPropagatorIterativeI1995, makriTensorPropagatorIterativeII1995](@cite)
    - Blip QuAPI [makriBlipDecompositionPath2014](@cite)
    - Time-Evolved Matrix Product Operators (TEMPO) [strathearnEfficientNonMarkovianQuantum2018](@cite)
    - Pairwise-Connected Tensor Network Path Integral (PC-TNPI) [bosePairwiseConnectedTensor2022](@cite)
- Hierarchical Equations of Motion (HEOM) [tanimuraNumericallyExactApproach2020](@cite)
- Generalized Quantum Master Equation (GQME) [nakajimaQuantumTheoryTransport1958, zwanzigIdentityThreeGeneralized1964](@cite)
- Multichromophore Incoherest Forster Resonance Energy Transfer [forsterZwischenmolekulareEnergiewanderungUnd1948, jangMultichromophoricForsterResonance2004](@cite)
- Bloch-Redfield Master Equation
- Transfer Tensor Method [cerrilloNonMarkovianDynamicalMaps2014](@cite) coupled with any of the path integral methods