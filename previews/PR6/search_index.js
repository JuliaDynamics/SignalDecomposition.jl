var documenterSearchIndex = {"docs":
[{"location":"#SignalDecomposition.jl-1","page":"Documentation","title":"SignalDecomposition.jl","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"SignalDecomposition is a simple Julia package that offers a single function:","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"decompose","category":"page"},{"location":"#SignalDecomposition.decompose","page":"Documentation","title":"SignalDecomposition.decompose","text":"decompose([t, ] s, method::Decomposition) → x, r\n\nDecompose an 1D input signal or timeseries s(t) into components, x, r, using the given method. t (the axis of s) defaults to 1:length(s).\n\nWhat are x and r really depend on your point of view and your application. They can be structure x and noise r (i.e. noise reduction). They can be seasonal/periodic x and residual component r. They can even be multiplier x and input r.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Documentation","title":"Documentation","text":"Feel free to add more methods to decompose by opening a PR. SignalDecomposition.jl is very open in the methods it offers and there is no plan to limit its dependencies.","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"In the test folder, there are examples of all these methods. You can comment out the plotting code to see the visual result of applying each method.","category":"page"},{"location":"#Linear-methods-1","page":"Documentation","title":"Linear methods","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"Here \"linear\" means linear in frequency space. For most methods you could of course get more control of the process by directly using DSP.jl (but then you won't be using the quite convenient decompose interface).","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"Fourier\r\nFrequencySplit\r\nSinusoidalFit","category":"page"},{"location":"#SignalDecomposition.Fourier","page":"Documentation","title":"SignalDecomposition.Fourier","text":"Fourier([s, ] frequencies) <: Decomposition\n\nDecompose a timeseries s into a sum x + r, by identifying specific frequencies at the Fourier space and removing them from the signal. x is the removed periodic component while r is the residual. If a given frequency is not exactly matching the Fourier frequencies, the closest one is removed.\n\nIf you provide s the method plans the forward and inverse fourier transforms (so that it is efficient to re-use it for s of same type and length).\n\nThis method works well when a periodic signal P is superimposed on fluctuations S, and you have a good educated guess of what frequencies compose P. This method works if the given signal has length multiple of the largest period given.\n\nImportant: periods/frequencies are defined with respect to the time axis length, the actual time axis t is not used in this method. So, frequency 1/12 (a period of 12) means 12 time points (whose actual time value depends on t).\n\n\n\n\n\n","category":"type"},{"location":"#SignalDecomposition.FrequencySplit","page":"Documentation","title":"SignalDecomposition.FrequencySplit","text":"FrequencySplit([s, ] f::Real) <: Decomposition\n\nSimilar to the Fourier method, but now the \"residual\" signal is the part of s with frequencies higher than f, while the \"seasonal\" part has frequencies ≤ f.\n\n\n\n\n\n","category":"type"},{"location":"#SignalDecomposition.SinusoidalFit","page":"Documentation","title":"SignalDecomposition.SinusoidalFit","text":"SinusoidalFit(s, fs) <: Decomposition\n\nDecompose a timeseries s into a sum x + r (with x the periodic component), by fitting sinuisoidals with given frequencies fs to the signal s using the package LsqFit. Specifically, fit\n\ns - bars approx sum_i A_i sin(2pi f_i t + phi_i)\n\nwith bars the mean. The fit happens on the amplitudes and phases A_i phi_i. After the decomposition you can find these in the struct's fields amps, phases. The fit is done on s versus t, so be sure that you have transformed t appropriately (e.g. if t is \"days\" but your frequencies are multiples of years, then you should give t/365.26).\n\nSinusoidalFit(fs, φ0s, A0s [, ub, lb])\n\nThe quality of the fit depends dramatically on the initial guesses for the phases and amplitudes, φ0s, A0s. This second constructor gives full control over initial phases and amplitudes, as well as upper and lower bounds on the amplitudes ub, lb (also vectors of length(fs)). In the first constructor φ0 = 0, A0 = abs(-(extrema(s)...))/2, ub = Inf, lb = -Inf for all frequencies. (the bounds for the phases are always ± π)\n\nNotice: LsqFit performs poorly for fitting sinusoidals and Fourier should be preferred over this method if the signal given is in multiples of the expected periods (and of course t is equally spaced).\n\n\n\n\n\n","category":"type"},{"location":"#Nonlinear-methods-1","page":"Documentation","title":"Nonlinear methods","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"ExtremelySimpleNL","category":"page"},{"location":"#SignalDecomposition.ExtremelySimpleNL","page":"Documentation","title":"SignalDecomposition.ExtremelySimpleNL","text":"ExtremelySimpleNL(k::Int, ℓ::Int, w::Int, τ::Int, ε::Real) <: Decomposition\n\nQuite literally the \"extremely simple nonlinear noise-reduction method\"[Schreiber1993]. It decomposes s into the sum x + r with x being the \"noiseless\" timeseries. This is the average position of neighbors in the delay embedded space.\n\nThis method works well if your timeseries is composed by the addition of a structured component (which follows deterministic and stationary dynamics which the embedding should approximate) and some noise.\n\nThe cited paper has some info on choosing optimal ε.\n\nArguments:\n\nk amount of past delay\nℓ amount of forward delay\nw Theiler window\nτ delay time\nε radius of the neighborhood in the embedded space\n\n[Schreiber1993]: Schreiber, (1993) Extremely simple nonlinear noise-reduction method. Physical Review E, 47(4)\n\n\n\n\n\n","category":"type"},{"location":"#Product-methods-1","page":"Documentation","title":"Product methods","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"ProductInversion","category":"page"},{"location":"#SignalDecomposition.ProductInversion","page":"Documentation","title":"SignalDecomposition.ProductInversion","text":"ProductInversion(r, μ; verbose=false) <: Decomposition\n\nDecompose a timeseries s into a product x * r, given that you have a good estimate of the second factor r (the \"input\") and you need x (the \"multiplier\") but you can't do simply x =  s ./ r because r contains zeros.\n\nThis method works well when the characteristic timescales of x are comparable, or larger than those of r but not much smaller.\n\nThe second argument μ is a regularization parameter. In short, we estimate r by minimizing a cost with two components: that x is close to s/r and that x is smooth. μ is the multiplier of the smoothness cost.\n\nYou can give a vector as μ. The process will be repeated for all μ and the rmse between s and the estimated x * r will be computed. The x that gives the least error will be returned finally. If verbose = true, the method also prints the pairs (μ, err) for each μ.\n\nUse the low level SignalDecomposition.matrix_invert(s, r, μ::Real) → x, err to get the error values.\n\n\n\n\n\n","category":"type"},{"location":"#Utilities-1","page":"Documentation","title":"Utilities","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"Simple utility functions to check how good the decomposition is for your data:","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"rmse\r\nnrmse\r\nσrmse","category":"page"},{"location":"#SignalDecomposition.rmse","page":"Documentation","title":"SignalDecomposition.rmse","text":"rmse(x, y) → e\n\nReturn the root mean square error e of the \"fit\" y into data x.\n\n\n\n\n\n","category":"function"},{"location":"#SignalDecomposition.nrmse","page":"Documentation","title":"SignalDecomposition.nrmse","text":"nrmse(x, y) → e\n\nReturn the normalized root mean square error of the \"fit\" y into data x. This number is the relative error of y to x versus mean(x) to x, i.e. if e < 1 the fit y is better than using mean(x) as a fit.\n\n\n\n\n\n","category":"function"},{"location":"#SignalDecomposition.σrmse","page":"Documentation","title":"SignalDecomposition.σrmse","text":"σrmse(x, y) = rmse(x, y)/std(x)\n\nRelative error of the fit y to data x with respect to std(x). \n\n\n\n\n\n","category":"function"}]
}
