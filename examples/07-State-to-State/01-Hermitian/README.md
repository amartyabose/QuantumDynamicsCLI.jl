# State-to-State Transport Analysis for Fenna-Matthews-Olson Complex

Here we use the HEOM method to first get the dynamics for the FMO model studied by Ishizaki and Fleming in 2009. With the dynamics in hand, we generate the state-to-state transport flows, $P_{j\leftarrow k}(t)$, for FMO. 

1. Run the HEOM simulation of the dynamics:
    ```bash
    qdsim simulate run system.toml simulate.toml
    ```
2. Finally, use the dynamics obtained to get the state-to-state flows:
    ```bash
    qdsim post state-to-state system.toml observables.toml
    ```

Notice that setting `derivative = true` in the `observables.toml` file also generates the time derivative of state-to-state flows, $\dot P_{j\leftarrow k}(t)$.