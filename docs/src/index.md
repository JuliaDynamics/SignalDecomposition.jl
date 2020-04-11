# SignalDecomposition.jl

`SignalDecomposition` is a Julia package with one purpose.

The function [`estimate_period`](https://juliadynamics.github.io/DynamicalSystems.jl/dev/chaos/periodicity/#ChaosTools.estimate_period) of `ChaosTools` is useful for some methods.

```@contents
Depth = 3
```

## Overview
SignalDecomposition.jl offers a single function:
```@docs
decompose
```

## Methods
Some of the following methods are linear and you could of course get more control of the process by directly using [DSP.jl](https://github.com/JuliaDSP/DSP.jl/) (but then you won't be using the quite convenient `decompose` interface).
```@autodocs
Modules = [SignalDecomposition]
Order   = [:type]
```

## Utilities
Simple utility functions to check how well the decomposition is for your data:
```@docs
rmse
nrmse
```
