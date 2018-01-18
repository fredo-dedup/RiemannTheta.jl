≈(a,b) = isapprox(a, b, rtol=1e-3, atol=1e-8)

@testset "values of riemanntheta" begin

    # using the Wolfram example
    # https://reference.wolfram.com/language/ref/SiegelTheta.html

    Ω = [ im -0.5 ; -0.5 im ]
    zs = [Complex128[0.5, 0.]]

    res = riemanntheta(zs, Ω)

    @test real(res[1]) ≈ 1.00748
    @test imag(res[1]) ≈ 0.

    zs = [ Complex128[x, 2x] for x in -1:0.01:1 ]
    res = riemanntheta(zs, Ω, eps=1e-3)

    @test maximum(real, res) ≈ 1.165
    @test minimum(real, res) ≈ 0.901
    @test maximum(imag, res) ≈ 0.
    @test minimum(imag, res) ≈ 0.

    # using example in abelfunctions doc
    # https://github.com/abelfunctions/abelfunctions/blob/master/doc/GettingStarted.md

    Ω = [ -1.309017+0.951057im -0.809017+0.587785im ;
          -0.809017+0.587785im -1.000000+1.175571im ]
    z = [0.5, 0.5im]
    res = riemanntheta([z], Ω, eps=1e-3)

    @test real(res[1]) ≈ 1.11415
    @test imag(res[1]) ≈ 0.8824
end

################################################################################

g = 3
δ = 1e-8
tmp = rand(g,g) ; Ω = Complex.(rand(g, g), tmp*tmp')
z₀ = rand(Complex128, g) - Complex(0.5, 0.5)

@testset "oscillatory_part derivs are correct" begin

    circvec = [Complex(1.,0.); zeros(Complex128, g-1)]

    # calculate function at slightly shifted z₀
    z = [[z₀] ; [z₀ + circshift(δ * circvec, i) for i in 0:g-1]]
    res = RiemannTheta.oscillatory_part(z, Ω)
    dres = [ (res[i] - res[1]) / δ for i in 2:g+1 ]

    # and compare to calculated derivates
    for i in 0:g-1
        derivs = [ circshift(circvec, i) ]
        res2 = RiemannTheta.oscillatory_part([z₀], Ω, derivs=derivs)
        # println(res2[1], dres[i+1])
        # res2[1]/δ , dres[i+1]/δ
        @test res2[1] ≈ dres[i+1]
    end

end


#### does not work, FIXME

@testset "riemanntheta derivs are correct" begin

    circvec = [Complex(1.,0.); zeros(Complex128, g-1)]

    # calculate function at slightly shifted z₀
    z = [[z₀] ; [z₀ + circshift(δ * circvec, i) for i in 0:g-1]]
    res = riemanntheta(z, Ω)
    dres = [ (res[i] - res[1]) / δ for i in 2:g+1 ]

    # and compare to calculated derivates
    for i in 0:g-1
        derivs = [ circshift(circvec, i) ]
        res2 = riemanntheta([z₀], Ω, derivs=derivs)
        # println(res2[1], dres[i+1])
        # res2[1]/δ , dres[i+1]/δ
        @test res2[1] ≈ dres[i+1]
    end
end
