using Documenter
using QuantumDynamicsCLI

makedocs(
    modules=[QuantumDynamicsCLI],
    sitename="QuantumDynamicsCLI.jl",
    pages=[
        "Introduction" => "index.md",
        "Documentation" => [
            "Simulate Module" => "./documentation/Simulate.md",
            "Post Module" => "./documentation/Post.md",
            "Comonicon" => "./documentation/ComonIcon.md",
        ]
        # "Tutorial" => [
        #     "Empirical Approaches" => "./tutorial/EmpiricalApproaches.md"
        #     "Path Integrals" => "./tutorial/BasicPI.md"
        #     "Bloch-Redfield Master Equation" => "./tutorial/Bloch-Redfield.md"
        #     "Quantum-Classical Path Integral" => "./tutorial/QCPI.md"
        #     "Dynamics under External Fields" => "./tutorial/ExternalFieldDynamics.md"
        #     "Quantum Correlation Functions" => "./tutorial/QuantumCorrelationFunction.md"
        #     "Hierarchical Equations of Motion" => "./tutorial/HEOM.md"
        # ],
        # "Documentation" => [
        #     "Bare System Propagation" => "./documentation/Bare.md",
        #     "Spectral Densities" => "./documentation/SpectralDensities.md",
        #     "Bloch-Redfield Master Equation" => "./documentation/BlochRedfield.md",
        #     "Incoherent Forster Theory" => "./documentation/Forster.md",
        #     "Influence Functional Coefficients" => "./documentation/EtaCoefficients.md",
        #     "Propagators" => "./documentation/Propagators.md",
        #     "Quasi-Adiabatic Propagator Path Integral" => "./documentation/QuAPI.md",
        #     "Blip Decomposition" => "./documentation/Blip.md",
        #     "Quantum-Classical Path Integral" => "./documentation/QCPI.md",
        #     "Time-Evolving Matrix Product Operators" => "./documentation/TEMPO.md",
        #     "Pairwise-Connected Tensor Network Path Integral" => "./documentation/PCTNPI.md",
        #     "Transfer Tensor Method" => "./documentation/TTM.md",
        #     "Generalized Quantum Master Equation" => "./documentation/GQME.md",
        #     "Complex-time Path Integrals" => "./documentation/CorrelationFunction.md",
        #     "Hierarchical Equations of Motion" => "./documentation/HEOM.md",
        #     "Utilities" => "./documentation/Utilities.md",
        #     "Compiling QuantumDynamics.jl" => "./documentation/Compile.md",
        # ],
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamicsCLI.jl.git"
)