using SignalDecomposition, Test
using Dates, Random

@testset "TimeAnomaly" begin
    Random.seed!(41516)
    y = Date(2001):Day(1):Date(2025)
    dy = dayofyear.(y)
    cy =  @. 4 + 7.2cos(2π*dy/365.26) + 5.6cos(4π*dy/365.26 + 3π/5)
    r0 = randn(length(dy))/2
    sy = cy .+ r0

    x, r = decompose(y, sy, TimeAnomaly())

    errper = nrmse(cy, x)
    @test errper < 0.1
    errres = nrmse(r0, r)
    @test nrmse(r0, r) < 0.5

    # figure()
    # title("TimeAnomaly")
    # ax1 = subplot(211)
    # plot(y, cy; alpha = 0.75, label = "original")
    # plot(y, x; alpha = 0.75, ls = "dashed", label = "err=$errper")
    # ylabel("x"); legend()
    # subplot(212; sharex = ax1)
    # plot(y, r0; alpha = 0.75, label = "original")
    # plot(y, r; alpha = 0.75, ls = "dashed", label = "err=$errres")
    # ylabel("r"); legend()
end
