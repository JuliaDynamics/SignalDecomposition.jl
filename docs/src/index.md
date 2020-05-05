# SignalDecomposition.jl

`SignalDecomposition` is a simple Julia package that offers a single function:
```@docs
decompose
```
Feel free to add more methods to `decompose` by opening a PR.
SignalDecomposition.jl is very open in the methods it offers and there is no plan to limit its dependencies.

In the `test` folder, there are examples of all these methods. You can comment out the plotting code to see the visual result of applying each method.

## Linear methods
Here "linear" means linear in frequency space. For most methods you could of course get more control of the process by directly using [DSP.jl](https://github.com/JuliaDSP/DSP.jl/) (but then you won't be using the quite convenient `decompose` interface).
```@docs
Fourier
FrequencySplit
Sinusoidal
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
Simple utility functions to check how good the decomposition is for your data:
```@docs
rmse
nrmse
σrmse
```
