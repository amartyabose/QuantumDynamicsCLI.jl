# Spin-Boson Example
Let us say we are trying to simulate a typical Spin-Boson problem, where all parameters are specified in atomic units.

```toml
[system]
Hamiltonian = "Hamiltonian"

[baths]
beta = 5.0
[[baths.bath]]
type = "ohmic"
xi = 0.1
omegac = 7.5
svec = [1.0, -1.0]
```

Notice that the `[units]` section is completely skipped over because the default values specify atomic units.