var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Only a few examples are shown here. Every method has an example (and plotting code) in the test folder!","category":"page"},{"location":"examples/#Nonlinear-1","page":"Examples","title":"Nonlinear","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"using SignalDecomposition, DynamicalSystems, Random, Plots, Statistics\n\nhe = Systems.henon()\ntr = trajectory(he, 10000; Ttr = 100)\nRandom.seed!(151521)\nz = tr[:, 1]\ns = z .+ randn(10001)*0.1*std(z)\nm = 5\nw = 1 # theiler window\nmetric = Euclidean()\nk = 30\nQ = [2, 2, 2, 3, 3, 3, 3]\nx, r = decompose(s, ManifoldProjection(m, Q, k))","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"p1 = plot(s, label = \"input\")\nplot!(p1, z, color = :black, ls = :dash, label = \"real\")\nplot!(p1, x, alpha = 0.5, label = \"output\")\nxlims!(p1, 0, 50)\np1","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Alright, this doesn't seem much of a difference to be honest. One sees a big difference once going into the state space and looking at the attractor:","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"p2 = scatter(s[1:end-1], s[2:end], ms = 1, label = \"input\", msw = 0)\nscatter!(p2, z[1:end-1], z[2:end], ms = 1, label = \"real\", color = :black, msw = 0)\nscatter!(p2, x[1:end-1], x[2:end], ms = 1, label = \"output\", alpha = 0.5, msw = 0)\np2","category":"page"},{"location":"examples/#Time-and-Sinusoidal-1","page":"Examples","title":"Time and Sinusoidal","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"using SignalDecomposition, Dates, Random, Plots\nRandom.seed!(41516)\ny = Date(2001):Day(1):Date(2025)\ndy = dayofyear.(y)\ncy =  @. 4 + 7.2cos(2π*dy/365.26) + 5.6cos(4π*dy/365.26 + 3π/5)\nr0 = randn(length(dy))/2\nsy = cy .+ r0\n\nx, r = decompose(y, sy, TimeAnomaly())\n\nt = collect(1:length(y)) ./ 365.26 # true time in years\n\nx2, r2 = decompose(t, sy, Sinusoidal([1.0, 2.0]))\n\np3 = plot(t, sy, label = \"input\")\nplot!(p3, t, cy, label = \"true periodic\", color = :black, ls = :dash)\nplot!(p3, t, x, label = \"TimeAnomaly\", alpha = 1.0, color = :red)\nplot!(p3, t, x2, label = \"Sinusoidal\", alpha = 0.5, color = :green)\nxlabel!(p3, \"years\")\nxlims!(p3, 0, 1) # zoom in","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Allthough not immediatelly obvious from the figure, Sinusoidal performs better:","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"rmse(cy, x), rmse(cy, x2)","category":"page"},{"location":"#SignalDecomposition.jl-1","page":"Introduction","title":"SignalDecomposition.jl","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"(Image: SignalDecomposition.jl)","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"SignalDecomposition is a simple Julia package that offers a single function:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"decompose","category":"page"},{"location":"#SignalDecomposition.decompose","page":"Introduction","title":"SignalDecomposition.decompose","text":"decompose([t, ] s, method::Decomposition) → x, r\n\nDecompose an 1D input signal or timeseries s(t) into components, x, r, using the given method. t defaults to 1:length(s).\n\nWhat are x and r really depends on your point of view and your application. They can be structure x and noise r (i.e. noise reduction). They can be seasonal/periodic x and residual component r. They can even be multiplier x and input r.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Introduction","title":"Introduction","text":"Subtypes of Decomposition are listed in the rest of this page. Feel free to add more methods to decompose by opening a PR, SignalDecomposition.jl is very open in the methods it offers.","category":"page"},{"location":"#Linear-methods-1","page":"Introduction","title":"Linear methods","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"Here \"linear\" means linear in frequency space. For some methods you could of course get more control of the process by directly using DSP.jl.    ","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Fourier\r\nFrequencySplit\r\nSinusoidal","category":"page"},{"location":"#SignalDecomposition.Fourier","page":"Introduction","title":"SignalDecomposition.Fourier","text":"Fourier([s, ] frequencies, x=true) <: Decomposition\n\nDecompose a timeseries s into a sum x + r, by identifying specific frequencies at the Fourier space and removing them from the signal. x is the removed periodic component while r is the residual. If a given frequency is not exactly matching the Fourier frequencies, the closest one is removed.\n\nImportant: periods/frequencies are defined with respect to the t axis length, the actual t values are not used in this method. So, frequency 1/12 (a period of 12) means 12 data points (whose actual value depends on t).\n\nIf you provide s the method plans the forward and inverse Fourier transforms (so that it is efficient to re-use it for s of same type and length).\n\nThis method works well when a periodic signal P is superimposed on fluctuations S, and you have a good educated guess of what frequencies compose P. This method works well if the given signal has length multiple of the periods given.\n\nThere is arbitrarity of which part of the signal x, r gets the mean value of s, because it is deducted for a better fit. The argument x=true attributes it to x, use false for r.\n\n\n\n\n\n","category":"type"},{"location":"#SignalDecomposition.FrequencySplit","page":"Introduction","title":"SignalDecomposition.FrequencySplit","text":"FrequencySplit([s, ] f::Real) <: Decomposition\n\nSimilar to the Fourier method, but now the \"residual\" signal is the part of s with frequencies higher than f, while the \"seasonal\" part has frequencies ≤ f.\n\n\n\n\n\n","category":"type"},{"location":"#SignalDecomposition.Sinusoidal","page":"Introduction","title":"SignalDecomposition.Sinusoidal","text":"Sinusoidal(fs)\n\nDecompose a timeseries s into a sum x + r, where x are sinusoidal components with the given frequencies fs that minimize coefficients A phi of the expression\n\ns approx A_0 + sum_i A_i cos(2pi f_i t + phi_i)\n\nwith bars the mean.\n\nThis method uses a new least-squares algorithm in frequency domain using the package LPVSpectral.jl, see[Bagge2017]. It works for non-equispaced t axis (and also normal), is generally very accurate (if choosen frequencies are not too close), but has performance scaling of O(N^2.4) instead of O(n log(n)) of Fourier.\n\nBecause it can work with arbitrary signal length the method always estimates the zero-frequency Fourier component, and attributes it to x. The fitted coefficients A phi are available as fields .A and .φ of the struct (first entry is zero-frequency component, i.e. the mean with respect to the sinusoidals).\n\n[Bagge2017]: F. Bagge Carlson et al., Linear Parameter-Varying Spectral Decomposition.\n\n\n\n\n\n","category":"type"},{"location":"#Nonlinear-methods-1","page":"Introduction","title":"Nonlinear methods","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"ExtremelySimpleNL\r\nManifoldProjection","category":"page"},{"location":"#SignalDecomposition.ExtremelySimpleNL","page":"Introduction","title":"SignalDecomposition.ExtremelySimpleNL","text":"ExtremelySimpleNL(k::Int, ℓ::Int, w::Int, τ::Int, ε::Real) <: Decomposition\n\nQuite literally the \"extremely simple nonlinear noise-reduction method\"[Schreiber1993]. It decomposes s into the sum x + r with x being the \"noiseless\" timeseries. This is the average position of neighbors in the delay embedded space.\n\nThis method works well if your timeseries is composed by the addition of a structured component (which follows deterministic and stationary dynamics which the embedding should approximate) and some noise.\n\nThe cited paper has some info on choosing optimal ε.\n\nArguments:\n\nk amount of past delay\nℓ amount of forward delay\nw Theiler window\nτ delay time\nε radius of the neighborhood in the embedded space\n\n[Schreiber1993]: Schreiber, (1993) Extremely simple nonlinear noise-reduction method. Physical Review E, 47(4)\n\n\n\n\n\n","category":"type"},{"location":"#SignalDecomposition.ManifoldProjection","page":"Introduction","title":"SignalDecomposition.ManifoldProjection","text":"ManifoldProjection(m, Q, k, τ=1, w=0, r=1000.0) <: Decomposition\n\nA nonlinear noise reduction method, also known as \"locally linear projections\", which works by bringing a noisy signal closer to a multi-dimensional manifold that represents the deterministic dynamics of the signal. The method is method \"IV\" of [Grassberger1993].\n\nm::Int is the same as in [Grassberger1993], the embedding dimension - 1. Q is related with the expected dimension d of the manifold of the deterministic dynamics, with d = m-Q+1. If given a Vector{Int} as Q the algorithm will iteratively do noise reduction to the resulting outputs (thus a vector is strongly recommended). Duplicate entries can exist in Q.\n\nk can be either Int or a SearchType from Neighborhood.jl. If Int, the k nearest neighbors are choosen as the neighborhood 𝓤 of each point. The paper contains an involved process for determining optimal k::Int, see eq.(5.4). w is just the Theiler window, while r is the value of the edge entries of vector R, and probably has not much impact (the rest of the entries are 1).\n\nIn the paper too big correction vectors were rescaled to the average magnitude of corrections, using as a criterion the distribution of their size. This is not implemented here (as it is not clear exactly what it means computationally, what is \"too big\"?) Contributing it is welcomed if you know how...\n\nSee also [Schreiber1996] for an application of the same algorithm in real ECG data.\n\n[Schreiber1996]: Schreiber & Kaplan (1996). Nonlinear noise reduction for electrocardiograms. Chaos, 6(1), 87–92\n\n[Grassberger1993]: Grassberger et al., (1993). On noise reduction methods for chaotic data. Chaos 3(2), 127–141\n\n\n\n\n\n","category":"type"},{"location":"#Product-methods-1","page":"Introduction","title":"Product methods","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"ProductInversion","category":"page"},{"location":"#SignalDecomposition.ProductInversion","page":"Introduction","title":"SignalDecomposition.ProductInversion","text":"ProductInversion(r, μ; verbose=false) <: Decomposition\n\nDecompose a timeseries s into a product x * r, given that you have a good estimate of the second factor r (the \"input\") and you need x (the \"multiplier\") but you can't do simply x =  s ./ r because r contains zeros.\n\nThis method works well when the characteristic timescales of x are comparable, or larger than those of r but not much smaller.\n\nThe second argument μ is a regularization parameter. In short, we estimate r by minimizing a cost with two components: that x is close to s/r and that x is smooth. μ is the multiplier of the smoothness cost.\n\nYou can give a vector as μ. The process will be repeated for all μ and the rmse between s and the estimated x * r will be computed. The x that gives the least error will be returned finally. If verbose = true, the method also prints the pairs (μ, err) for each μ.\n\nUse the low level SignalDecomposition.matrix_invert(s, r, μ::Real) → x, err to get the error values.\n\n\n\n\n\n","category":"type"},{"location":"#Miscellaneous-methods-1","page":"Introduction","title":"Miscellaneous methods","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"TimeAnomaly","category":"page"},{"location":"#SignalDecomposition.TimeAnomaly","page":"Introduction","title":"SignalDecomposition.TimeAnomaly","text":"TimeAnomaly()\n\nDecompose timeseries s into its temporal average x and the anomalies r so that s = x + r. Each unique day+month combination in t is identified, and the values of s for each year that has this day+month combination are averaged. As a result, the time vector t must be <:AbstractVector{<:TimeType}.\n\nThis method is very common in climate science, referred to as simply \"anomalies\".\n\n\n\n\n\n","category":"type"},{"location":"#Utilities-1","page":"Introduction","title":"Utilities","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"Simple utility functions to check how good the decomposition is for your data:","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"rmse\r\nnrmse\r\nσrmse","category":"page"},{"location":"#SignalDecomposition.rmse","page":"Introduction","title":"SignalDecomposition.rmse","text":"rmse(x, y) → e\n\nReturn the root mean square error e of the \"fit\" y into data x.\n\n\n\n\n\n","category":"function"},{"location":"#SignalDecomposition.nrmse","page":"Introduction","title":"SignalDecomposition.nrmse","text":"nrmse(x, y) → e\n\nReturn the normalized root mean square error of the \"fit\" y into data x. This number is the relative error of y to x versus mean(x) to x, i.e. if e < 1 the fit y is better than using mean(x) as a fit.\n\n\n\n\n\n","category":"function"},{"location":"#SignalDecomposition.σrmse","page":"Introduction","title":"SignalDecomposition.σrmse","text":"σrmse(x, y) = rmse(x, y)/std(x)\n\nRelative error of the fit y to data x with respect to std(x). \n\n\n\n\n\n","category":"function"}]
}
