@testset "Pair distance histograms" begin
    @testset "Simple cubic 2D" begin
        L = 3
        δr = 0.01
        box = 2.0 * (L + δr)
        simplecubic2D = vec([Particle(r = SVector{Float64}(i,j,0))
                             for i in -L:L, j in -L:L])
        idx = findfirst(a -> a.r == [0, 0, 0], simplecubic2D)
        reference_atom = simplecubic2D[idx:idx]
        other_atoms = [simplecubic2D[i] for i in 1:(2*L+1)^2 if i ≠ idx]

        rpoints, hist2D = histradial(reference_atom, other_atoms, box, δr)
        @test hist2D[hist2D .≠ 0.0] ==
            [4.0, 4.0, 4.0, 8.0, 4.0, 4.0]
    end # testset

    @testset "Simple cubic 3D" begin
        L = 2
        δr = 0.01
        box = 10.0 * (L + δr)
        simplecubic3D = vec([Particle(r = SVector(i,j,k))
                             for i in -L:L, j in -L:L, k in -L:L])
        idx = findfirst(a -> a.r == [0, 0, 0], simplecubic3D)
        reference_atom = simplecubic3D[idx:idx]
        other_atoms = [simplecubic3D[i] for i in 1:(2*L+1)^3 if i ≠ idx]

        rpoints, hist3D = histradial(reference_atom, other_atoms, box, δr)
        @test hist3D[hist3D .≠ 0.0] ==
            [6.0, 12.0, 8.0, 6.0, 24.0, 24.0, 12.0, 24.0, 8.0]
    end # testset
end # testset
