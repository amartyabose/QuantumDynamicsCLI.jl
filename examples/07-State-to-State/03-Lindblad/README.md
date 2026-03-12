# Lindblad State-to-State Transport Analysis for a Pumped-Drained Excitonic Dimer

Lindblad jump operators can be used to describe excitation pumps and drains in excitonic aggregates. Here we use the TEMPO-TTM method to first get the dynamics for an excitonic dimer coupled to site-based baths where the first monomer has a pump of timescale $300\,\,\text{fs}$ while the second monomer has a drain of timescale of $300\,\,\text{fs}$. With the dynamics in hand, we generate the state-to-state transport flows, $P_{j\leftarrow k}(t)$, for the dimer. 

1. Run the TEMPO-TTM simulation of the dynamics:
    ```bash
    qdsim simulate run system.toml simulate.toml
    qdsim simulate propagate-using-tmats system.toml propagate.toml
    ``` 
2. Finally, use the dynamics obtained to get the state-to-state flows:
    ```bash
    qdsim post state-to-state system.toml observables.toml
    ```

The Linblad jump operators are specified using the label `lindblad` with the corresponding timescale decribed by `decay_constant` in the TOML files. It is important to note that the valid Lindblad jump operators are restricted to form:

$$
L = \sum_n c_n |f_n\rangle\langle i_n|,
\;
n \neq m \Rightarrow |f_n\rangle \neq |f_m\rangle
$$

Also, note that the Hamiltonian here describes the full Hilbert space of the dimer.