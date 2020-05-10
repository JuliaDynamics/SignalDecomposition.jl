# SignalDecomposition.jl

`SignalDecomposition` is a simple Julia package that offers a single function:
```@docs
decompose
```
Subtypes of `Decomposition` are listed in the rest of this page.
Feel free to add more methods to `decompose` by opening a PR,
SignalDecomposition.jl is very open in the methods it offers.

## Linear methods
Here "linear" means linear in frequency space. For some methods you could of course get more control of the process by directly using [DSP.jl](https://github.com/JuliaDSP/DSP.jl/).    
```@docs
Fourier
FrequencySplit
Sinusoidal
```

## Nonlinear methods
```@docs
ExtremelySimpleNL
ManifoldProjection
```

## Product methods
```@docs
ProductInversion
```

## Miscellaneous methods
```@docs
TimeAnomaly
```

## Utilities
Simple utility functions to check how good the decomposition is for your data:
```@docs
rmse
nrmse
σrmse
```
