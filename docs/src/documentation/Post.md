# qdsim post

Once the simulation is done using the `simulate` subcommand of `qdsim`, various post-processing can be done using its `post` subcommand.

## qdsim post get-observable

### Post-processing of `dynamics` simulation
This subcommand calculates the expectation of the specified observables using the reduced density matrix previously calculated with the `dynamics` simulation using [the `simulate` module](@ref "qdsim simulate"). As with the `simulate` subcommand, this command also requires the `system` TOML file. Along with it, one needs the `observables` TOML file which specifies both the simulation details used to compute the density matrix beforehand, and the observables themselves. A typical `observables` TOML file looks like the following:

```toml
[[simulation]]
name = "model_a"
method = "Spin-LSC"
output = "model_a.hdf5"

dt = 0.1
nsteps = 200

SW_transform = "WTransform"
focused_sampling = true

num_bins = 100
num_mc = 1000

rho0 = "P_1"
outgroup = "init1"

observable_output = "observables.txt"

[[simulation.observable]]
observable = "P_1"
[[simulation.observable]]
observable = "./sx"
[[simulation.observable]]
observable = "./sy"
type = "complex"
```

The TOML file must contain the `simulation` section, within with each observable is specified. Apart from the simulation details in the `simulation` section (such as the specification of the Stratonovich--Weyl kernel, time step, etc.), the specified parameter `observable_output` tells the output filename for the observables.

The subsections `simulation.observable` specify the observable to calculate the expectation value for. These can contain the following parameters:
- `observable`: this denotes the observable to compute and is like the `rho0` keyword in the `simulation` section (see also [QuantumDynamicsCLI.ParseInput.parse_operators](@ref)), but other quantities may be computed by specifying:
  - `observable = "trace"`: computes the trace of the reduced density matrix
  - `observable = "purity"`: computes the purity of the reduced density matrix
  - `observable = "vonNeumann_entropy"`: computes the von Neumann entropy for the calculated density matrix
- `type`: this parameter may be either "real" or "complex" to denote whether the `observable` is a real or complex matrix. By default, it is assumed that it specifies a real matrix.
<!-- TODO: - `fourier_transform`: whether to take the Fourier transform of the requested observables -->

Once the observables are calculated, the real and imaginary parts are separated stored in separate files e.g., `observables_real.txt` and `observables_imag.txt` in the given example. For the methods that do a Monte-Carlo average to compute the reduced density matrix, such as the semiclassical methods, the standard error of the expectation values are inserted as separate columns after their mean value for both the real and the imaginary files.

In the example above, the first observable is the population of state 1 specified via the shorthand. Whereas the other two observables are specified in their matrix form in the given file.

<!-- TODO: ### Post-processing of `complex_corr` simulation -->

## qdsim post merge-into

```@docs
QuantumDynamicsCLI.Post.merge_into
```
