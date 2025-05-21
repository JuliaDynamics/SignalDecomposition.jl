[![docsdev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/signaldecomposition/dev/)
[![docsstable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/signaldecomposition/stable/)
[![CI](https://github.com/JuliaDynamics/SignalDecomposition.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/SignalDecomposition.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/SignalDecomposition.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/SignalDecomposition.jl)

SignalDecomposition.jl is a Julia package providing an interface and dozens of algorithm implementations for signal decomposition.
Given a signal (or timeseries), the function `decompose` splits it into two components.
These may be:

- structure and noise (de-noising or smoothing)
- seasonal and residual (climatologies or anomalies)
- trend and residual (de-trending)

It can be used as a standalone package, or as part of
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/).

To install it, run `import Pkg; Pkg.add("SignalDecomposition")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/signaldecomposition/stable/) or build locally by running the `docs/make.jl` file.
