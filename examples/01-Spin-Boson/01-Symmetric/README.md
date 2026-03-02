# Asymmetric Spin-Boson example

In this folder, we have all the files necessary for running an asymmetric spin-boson calculation using various methods

1. Run the simulation with different methods:
    ```bash
    qdsim simulate run system.toml <simulation-file>
    ```
   The following are the different simulation files:
   1. TEMPO: simulate_TEMPO.toml
   2. PLDM: simulate_PLDM.toml
   3. Spin-PLDM: simulate_spinPLDM.toml
   4. LSC: simulate_LSC.toml
   5. Spin-LSC: simulate_spinLSC.toml
2. Finally, use the dynamics obtained to generate observables:
    ```bash
    qdsim post get-observable system.toml observables.toml
    ```
