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
	- Spin-Mapping Approaches to PLDM [mannouchPartiallyLinearizedSpinMapping2020a, mannouchPartiallyLinearizedSpinMapping2020a](@cite) and LSC-IVR
- Hierarchical Equations of Motion (HEOM) [tanimuraNumericallyExactApproach2020](@cite)
- Generalized Quantum Master Equation (GQME) [nakajimaQuantumTheoryTransport1958, zwanzigIdentityThreeGeneralized1964](@cite)
- Multichromophore Incoherest Forster Resonance Energy Transfer [forsterZwischenmolekulareEnergiewanderungUnd1948, jangMultichromophoricForsterResonance2004](@cite)
- Bloch-Redfield Master Equation
- Transfer Tensor Method (TTM) [cerrilloNonMarkovianDynamicalMaps2014](@cite) coupled with any of the path integral methods

All of these dynamics methods require some core common parameters and then more specfic method-dependent parameters. The core parameters of all the dynamics methods are:
- `dt`: for the time-step in the units specified in the system file
- `nsteps`: for the number of steps of simulation of the dynamics

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

We provide both the PLDM family of methods and the LSC-IVR family of
methods.

### Partially Linearized Density Matrix

Both PLDM and Spin-PLDM are implemented --- chosen by `method =
"PLDM"` and `method = "Spin-PLDM"` respectively.

### Linearized Semiclassical Initial Value Representation

A choice between standard LSC and spin LSC is provided: `method =
"LSC"` and `method = "Spin-LSC"` respectively.
