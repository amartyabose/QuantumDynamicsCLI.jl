# Non-Hermitian State-to-State Transport Analysis for a Leaky Excitonic Dimer

Non-hermitian Hamiltonians can be used to describe loss of excitation from an excitonic aggregate. Here we use the TEMPO-TTM method to first get the dynamics for an excitonic dimer coupled to site-based baths. The second monomer loses excitation with a site-specific timescale of $600\ \text{fs}$. With the dynamics in hand, we generate the state-to-state transport flows, $P_{j\leftarrow k}(t)$, for the dimer. 

1. Run the TEMPO-TTM simulation of the dynamics:
    ```bash
    qdsim simulate run system.toml simulate.toml
    qdsim simulate propagate-using-tmats system.toml propagate.toml
    ``` 
2. Finally, use the dynamics obtained to get the state-to-state flows:
    ```bash
    qdsim post state-to-state system.toml observables.toml
    ```

Note that $P_{j\leftarrow j}(t)$ terms represent the site-specific loss.