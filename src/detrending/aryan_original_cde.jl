using Polynomials
using Statistics
using LinearAlgebra
using EmpiricalModeDecomposition
using HPFilter
using Loess
using MarketData
using ComplexityMeasures

function closing_stock_timeseries(name::String, start::DateTime, finish::DateTime, freq::String)
    opt = YahooOpt(period1 = start, period2 = finish, interval = freq)
    stockdata = yahoo(name, opt)
    stock_close = stockdata["Close"]
    actual_numbers = values(stock_close)
    return float.(actual_numbers)
end

function Complexity_Measures(y; r = std(y)*0.2, a=6)
    permutation = entropy_normalized(OrdinalPatterns{4}(), y)
    spectral = entropy_normalized(PowerSpectrum(), y)
    dispersion = entropy_normalized(Dispersion(c = a), y)
    approximate = complexity(ApproximateEntropy(; r), y)
    sample = complexity(SampleEntropy(; r), y)
    return permutation, spectral, dispersion, approximate, sample
end

function detrend_ols(y)
    n = length(y)
    X = hcat(ones(n), collect(1:n))
    β = X \ y  # OLS estimation
    trend = X * β  # Fitted trend
    y_detrended = y - trend  # de-trend data
    return y_detrended
end

function polynomial_detrend(y, deg=2)
    x = 1:length(y)
    p = fit(x, y, deg)
    trend = p.(x)
    mse = sum((y - trend).^2) / length(y)
    y_detrended = y .-trend
    return y_detrended
end

function moving_average(y, window = 30)
    n = length(y)
    trend = zeros(n)

    for i in 1:n
        half = (window - 1) ÷ 2
        start_idx = max(1, i - half)
        end_idx = min(n, i + half)
        trend[i] = mean(y[start_idx:end_idx])
    end
    mse = sum((y - trend).^2) / n
    detrend = y .- trend
    return detrend
end

function detrend_emd(y)
    n = length(y)
    t = 1:n
    imf = emd(y, t)
    trend = imf[7] + imf[6]
    residual = sum(imf[i] for i in 1:5)
    detrend = y .- trend
    return detrend
end

function detrend_hp(y, lambda = 1600)
    n = length(y)
    trend = HP(y, lambda)
    detrend = y .- trend
    mse = sum((y .- trend).^2) / n
    return detrend
end

function detrend_loess(y, s = 0.2)
    n = length(y)
    x = 1:n
    model = loess(x, y, span = s)
    trend = predict(model, x)
    detrend = y .- trend
    mse = sum((y .- trend).^2) / n
    return detrend
end

function mutual_information(x, y, n=1000)
    x = copy(x)
    y = copy(y)
    Hx = entropy(ValueHistogram(10), Dataset(x))
    Hy = entropy(ValueHistogram(10), Dataset(y))
    Hxy = entropy(ValueHistogram(10), Dataset(x, y))
    mi = Hx + Hy - Hxy

    null = zeros(n)
    for i in 1:n
        shuffle!(x), shuffle!(y)
        Hxy = entropy(ValueHistogram(10), Dataset(x, y))
        null[i] = Hx + Hy - Hxy
    end

    mu = mean(null)
    sd = std(null)
    return mi, mu, sd
end


function main_analysis(timeseries, detrending, complexity_measure)
    detrended = detrending(timeseries)[1]
    measure = complexity_measure(detrended)
    return measure
end