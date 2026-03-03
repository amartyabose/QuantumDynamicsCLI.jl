# Dynamics Simulations

Various methods of simulation are supported:
- Path Integral Methods using Feynman-Vernon Influence Functional[feynmanTheoryGeneralQuantum1963](@cite):
    - Quasi-adiabatic Propagator Path Integrals (QuAPI) [makriTensorPropagatorIterativeI1995, makriTensorPropagatorIterativeII1995](@cite)
    - Blip QuAPI [makriBlipDecompositionPath2014](@cite)
    - Time-Evolved Matrix Product Operators (TEMPO) [strathearnEfficientNonMarkovianQuantum2018](@cite)
    - Pairwise-Connected Tensor Network Path Integral (PC-TNPI) [bosePairwiseConnectedTensor2022](@cite)
- Semiclassical Dynamics using Meyer-Miller-Stock-Thoss Hamiltonian
    - Partially Linearized Semiclassical Dynamics (PLDM) [huoCommunicationPartialLinearized2011, huoConsistentSchemesNonadiabatic2012](@cite)
    - Linearized Semiclassical Initial Value Representation (LSC-IVR)
    - Spin-Mapping Approaches to PLDM [mannouchPartiallyLinearizedSpinmapping2020a, mannouchPartiallyLinearizedSpinmapping2020](@cite) and LSC-IVR
- Hierarchical Equations of Motion (HEOM) [tanimuraNumericallyExactApproach2020](@cite)
- Generalized Quantum Master Equation (GQME) [nakajimaQuantumTheoryTransport1958, zwanzigIdentityThreeGeneralized1964](@cite)
- Multichromophore Incoherest Forster Resonance Energy Transfer [forsterZwischenmolekulareEnergiewanderungUnd1948, jangMultichromophoricForsterResonance2004](@cite)
- Bloch-Redfield Master Equation
- Transfer Tensor Method (TTM) [cerrilloNonMarkovianDynamicalMaps2014](@cite) coupled with any of the path integral methods

All of these dynamics methods require some core common parameters and then more specfic method-dependent parameters. The core parameters of all the dynamics methods are:
- `dt`: for the time-step in the units specified in the system file
- `nsteps`: for the number of steps of simulation of the dynamics
- `rho0`: the initial reduced density matrix
- `outgroup`: where to store the computed density matrix in the HDF5 file (i.e., the group name)

## Specifying the initial reduced density matrix

The simplest way to specify the initial density matrix is to set the `rho0` parameter to a file which will be parsed as a matrix. However, convenient shortcuts exist to specify most commonly used values of `rho0`, these are explained in the documentation of [QuantumDynamicsCLI.ParseInput.parse_operator](@ref).

## Feynman-Vernon Influence Functional Simulations

There are two ways of incorporating the effect of non-Markovian memory
in path integral simulations: iterative propagation beyond memory or
using the transfer tensor method.

For the traditional iterative propagation technique, one can use:
- `method = "QuAPI"`: for using QuAPI with standard iteration
- `method = "TEMPO"`: for using TEMPO with standard iteration

To use TTM, one can choose one of the following:
- `method = "QuAPI-TTM"`: for using QuAPI within memory
- `method = "Blip-TTM"`: for using Blip within memory
- `method = "TEMPO-TTM"`: for using TEMPO within memory

Finally, if one does not intend on using TTM, we suggest using TEMPO
for accessing long memory lengths efficiently, though QuAPI is also
available. This is chosen by setting `method = "TEMPO"`.

In this family of methods, the non-Markovian memory is incorporated
explicitly by specifying the number of time-steps it spans. For all
the TTM-based methods, this memory length is set through the
parameter, `rmax`. For the methods with traditional iteration, it is
set through the parameter, `kmax`.

Every method has a separate set of parameters that handle the balance
between accuracy and efficiency. In case of QuAPI-based methods, that
parameter is called `cutoff` and it is by default set to
$10^{-10}$. All paths with absolute value of amplitude below this
cutoff are ignored.

In the case of "Blip-TTM", the accuracy is set via the maximum number
of blips allowed, `max_blips`. If this is set to $-1$, that means all
blips are included.

For "TEMPO-TTM", the accuracy is set through the tensor network
evaluation parameters. The `cutoff` of the singular value
decomposition, and the maximum bond dimension, `maxdim`, used in
obtaining the matrix product state representation of the path
amplitude tensor are the two parameters that are used for controlling
the accuracy. The default values of these two parameters are
$10^{-10}$ and $1000$ respectively.

## Semiclassical Simulations

We provide two mapping Hamiltonian based class of semiclassical methods: Meyer-Miller-Stock-Thoss mapping, and Spin-mapping. Both the PLDM and the LSC-IVR family of semiclassical methods are available. To use these methods, say `method = "$METHOD"` where `$METHOD` is one of:
- `LSC` or `PLDM` for the fully or partially linearized semiclassical dynamics using MMST mapping for the system degrees of freedom
- `Spin-LSC` or `Spin-PLDM` for the corresponding spin-mapped variants

As these methods perform a Monte-Carlo average to calculate the reduced density matrix, the following parameters must be set for all these methods:
- `num_bins`: the number of independent bins to calculate the average and standard deviation of the observables with
- `num_mc`: the number of trajectories for each bin

Apart from this, the spin-mapped semiclassical methods offer the choice of choosing the corresponding Stratonovich--Weyl kernel to use for transformation of the Hamiltonian of the system and the initial (reduced) density matrix. This is set by the `SW_transform` keyword, and the supported values are `QTransform`, `PTransform` and `WTransform`. By default, Spin-LSC uses `QTransform` and Spin-PLDM uses `WTransform`. **NOTE:** Spin-PLDM performs poorly with P and Q Stratonovich--Weyl kernels.

Focused initial sampling is currently only supported by Spin-LSC at the moment. Moreover, only `rho0` of the form $\ket{n}\bra{n}$ are supported currently. Focused sampling may be enabled for such initial reduced density matrices in Spin-LSC by setting the `focused_sampling` parameter to `true`.

These methods _must_ specify the number of discrete oscillators for each `bath` mode via the `num_osc` parameter as specified in the [Bath Hamiltonian](@ref) section.
