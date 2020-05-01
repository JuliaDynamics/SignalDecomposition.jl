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
Feel free to add more methods to `decmpose` by opening a PR.
SignalDecomposition.jl is very open in the methods it offers and there is no plan to limit its dependencies.

In the `test` folder, there are examples of all these methods. You can comment out the plotting code to see the visual result of applying each method.

## Linear methods
Here "linear" means linear in frequency space. For most methods you could of course get more control of the process by directly using [DSP.jl](https://github.com/JuliaDSP/DSP.jl/) (but then you won't be using the quite convenient `decompose` interface).
```@docs
Fourier
FrequencySplit
SinusoidalFit
```

## Nonlinear methods
```@docs
ExtremelySimpleNL
```

## Product methods
```@docs
ProductInversion
```

## Utilities
Simple utility functions to check how well the decomposition is for your data:
```@docs
rmse
nrmse
σrmse
```
