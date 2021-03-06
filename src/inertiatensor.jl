export
    inertiatensor,
    asphericity,
    elongation,
    prolateness,
    gyrationradius,
    asphericity2,
    prolateness2

δ(i, j) = i == j ? 1 : 0

function inertiatensor(particles::T, box) where T <:AbstractVector{Particle{3,Float64}}
    com = centerofmass(particles, box)
    coords = [distancevector(com, p.r, box) for p in particles]
    Iijterm(r,i,j) = δ(i,j) * norm(r)^2 - r[i] * r[j]
    Iij(i, j) = sum([Iijterm(r, i, j) for r in coords])
    I = SMatrix{3,3}([Iij(i,j) for i in 1:3, j in 1:3])
    return I
end # function

function inertiatensor(atoms::T, box) where T<:AbstractVector{Atom}
    natoms = size(atoms, 1)
    com = centerofmass(atoms, box)
    coords = [distancevector(com, a.r, box) for a in atoms]
    Iijterm(r,i,j) = δ(i,j) * norm(r)^2 - r[i] * r[j]
    Iij(i, j) = sum([
        atoms[n].m * Iijterm(coords[i], i, j) for n in 1:natoms
    ])
    I = SMatrix{3,3}([Iij(i,j) for i in 1:3, j in 1:3])
    return I
end # function

function asphericity(I::SMatrix{3,3,Float64})
    λ = sqrt.(3 .* eigvals(I))
    a = (λ[1] - λ[2])^2
    b = (λ[2] - λ[3])^2
    c = (λ[3] - λ[1])^2
    α = (a + b + c) / (2 * sum(λ)^2)
    return α
end # function

δ(k::Integer,j::Integer) = k == j ? 1 : 0

function asphericity2(I::SMatrix{3,3,Float64})
    Δ = [δ(k,j) for j in 1:3, k in 1:3]
    Q = I .- Δ * tr(I) ./ 3
    α = 3/2 * tr(Q^2) / tr(I)^2
    return α
end # function

function prolateness2(I::SMatrix{3,3,Float64})
    Δ = [δ(k,j) for j in 1:3, k in 1:3]
    Q = I .- Δ * tr(I) ./ 3
    p = 4 * det(Q) / (2/3 * tr(Q^2))^(3/2)
    return p
end # function

function elongation(I::SMatrix{3,3,Float64})
    λ = sqrt.(3 .* eigvals(I))
    ϵ = (λ[2] - λ[1]) / λ[3]
    return ϵ
end # function

function prolateness(I::SMatrix{3,3,Float64})
    λ = sqrt.(3 .* eigvals(I))
    a = 2*λ[1]-λ[2]-λ[3]
    b = 2*λ[2]-λ[3]-λ[1]
    c = 2*λ[3]-λ[1]-λ[2]
    d = 2*(λ[1]^2 + λ[2]^2 + λ[3]^2 - λ[1]*λ[2] - λ[2]*λ[3] - λ[3]*λ[1])
    p = a * b * c / d^(3/2)
    return p
end # function

function gyrationradius(I::SMatrix{3,3,Float64})
    return sum(eigvals(I))
end # function
