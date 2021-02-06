export
    inertiatensor,
    asphericity,
    elongation,
    prolateness

δ(i, j) = i == j ? 1 : 0

function inertiatensor(config::T) where T <:AbstractVector{Particle{3,Float64}}
    Iij(i, j) = sum([norm(p.r)^2*δ(i,j) - p.r[i]*p.r[j] for p in config])
    I = SMatrix{3,3}([Iij(i,j) for i in 1:3, j in 1:3])
    return I
end # function

function inertiatensor(config::T) where T<:AbstractVector{Atom}
    Iij(i, j) = sum([
        p.m*(norm(p.r)^2*δ(i,j) - p.r[i]*p.r[j]) for p in config
    ])
    I = SMatrix{3,3}([Iij(i,j) for i in 1:3, j in 1:3])
    return I
end # function

function asphericity(I::SMatrix{3,3,Float64})
    λ = eigvals(I)
    α = sum(abs2.(diff(vcat(λ, [λ[1]])))) / (2 * sum(λ)^2)
    return α
end # function

function elongation(I::SMatrix{3,3,Float64})
    λ = eigvals(I)
    ϵ = (λ[2] - λ[1]) / λ[3]
    return ϵ
end # function

function prolateness(I::SMatrix{3,3,Float64})
    λ = eigvals(I)
    a = 2*λ[1]-λ[2]-λ[3]
    b = 2*λ[2]-λ[3]-λ[1]
    c = 2*λ[3]-λ[1]-λ[2]
    d = 2*(λ[1]^2 + λ[2]^2 + λ[3]^2 - λ[1]*λ[2] - λ[2]*λ[3] - λ[3]*λ[1])
    p = a * b * c / d^(3/2)
    return p
end # function
