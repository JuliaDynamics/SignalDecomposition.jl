using SignalDecomposition, Test

@testset "detrending" begin
    detrending_methods = [
        PolynomialDetrending(1), PolynomialDetrending(2), NoDecomposition(),
        MovingAverageSmoothing(), LoessSmoothing(),
    ]

    N = 10000
    s = float.(1:N)
    s .+= 0.1rand(N)

    @testset "$(nameof(typeof(m)))" for m in detrending_methods
        trend, residual = decompose(s, m)
        err = rmse(s, trend)
        @test abs(err) < 1e-1 # noise is this magnitude.
    end
end
