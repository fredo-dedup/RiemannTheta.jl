using Base.Test

g, num_vectors = 5, 10
z = [ rand(Complex128, g) for i in 1:num_vectors]
mode, ϵ = 0, 1e-6
accuracy_radius = 5.
tttt = rand(5,5)
Y = tttt * tttt'
Ω = Complex.(rand(g, g), Y)

derivs = Vector{Complex128}[]
@btime oscillatory_part($z, $Ω, 0, $ϵ, $derivs, $accuracy_radius) # 380-390ms / 430 Mb

derivs = [ rand(Complex128, 5) for i in 1:5 ]
@btime oscillatory_part($z, $Ω, 0, $ϵ, $derivs, $accuracy_radius) # 6.9s / 5.6Gb


@testset "oscillatory_part derivs are correct" begin
    z₀ = rand(Complex128, g) - Complex(0.5, 0.5)
    δ = [Complex(1e-8,1e-8); zeros(Complex128, g-1)]
    derivs = Vector{Complex128}[]
    ϵ = 1e-8
    # calculate function at slightly shifted z₀
    z = [[z₀] ; [z₀ + circshift(δ, i) for i in 0:g-1]]
    res = oscillatory_part(z, Ω, 0, ϵ, derivs, accuracy_radius)
    dres = [ res[i] - res[1] for i in 2:g+1 ]
    dres

    # and compare to calculated derivates
    for i in 0:g-1
        derivs = [ circshift(δ, i) ]
        res2 = oscillatory_part([z₀], Ω, 0, ϵ, derivs, accuracy_radius)
        @test isapprox(res2[1], dres[i+1], rtol=1e-4)
    end
end
