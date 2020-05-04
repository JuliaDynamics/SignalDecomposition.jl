peri_with_zero = copy(periodic)
peri_with_zero[periodic .≤ 0] .= 0

ro_with_zero = copy(10roesslerz)
ro_with_zero[10roesslerz .≤ 0.1] .= 0

for j in 1000:2000:9000
    ro_with_zero[j:1:j+200] .= 0
end

cases = [
    (lorenzx_slow, peri_with_zero)
    (lorenzx_slow, ro_with_zero)
]

@testset "ProductInversion" begin
    for (i, case) in enumerate(cases)

        s = case[1] .* case[2]
        m = ProductInversion(case[2], 0.1:0.1:10)
        x, r = decompose(s, m)
        errper = nrmse(case[1], x)
        @test errper < 0.5
        errres = nrmse(case[2], r) # by definition 0
        fullerr = nrmse(s, x .* r)
        @test fullerr < 0.1
        # println("  case $i errper=$errper, errres=$errres ")
        #
        # figure()
        # ax1 = subplot(311)
        # title(string(nameof(typeof(m)))*" case $i")
        # plot(s, label = "original")
        # plot( x .* r, label = "err=$(fullerr)")
        # ylabel("x*r"); legend()
        # subplot(312; sharex = ax1)
        # plot(case[1]; alpha = 0.75, label = "original")
        # plot(x; alpha = 0.75, ls = "dashed", label = "err=$errper")
        # ylabel("x (multiplier)"); legend()
        # subplot(313; sharex = ax1)
        # plot(case[2]; alpha = 0.75, label = "original")
        # plot(r; alpha = 0.75, ls = "dashed", label = "err=$errres")
        # ylabel("r (input)"); legend()
    end
end
