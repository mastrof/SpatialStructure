export
    centerofmass

# bai, breen 2008 journal of graphical tools
function centerofmass(particles::P, box::AbstractVector) where P<:AbstractVector{AP} where AP<:AbstractParticle{D,T} where {D,T}
    rcom = zeros(D)
    M = size(particles, 1)
    for d in 1:D
        θ = 2π / box[d] .* [p.r[d] for p in particles]
        ξ = cos.(θ)
        ζ = sin.(θ)
        ξcom = sum(ξ) / M
        ζcom = sum(ζ) / M
        θcom = atan(-ζcom, -ξcom) + π
        rcom[d] = box[d] * θcom / (2π)
    end # for
    return rcom
end # function

function centerofmass(atoms::T, box::AbstractVector) where T<:AbstractVector{Atom}
    D = 3
    rcom = zeros(D)
    mass = [a.m for a in atoms]
    M = sum(mass)
    for d in 1:D
        θ = 2π / box[d] .* [a.r[d] for a in atoms]
        ξ = cos.(θ)
        ζ = sin.(θ)
        ξcom = sum(ξ .* mass) / M
        ζcom = sum(ζ .* mass) / M
        θcom = atan(-ζcom, -ξcom) + π
        rcom[d] = box[d] * θcom / (2π)
    end # for
    return rcom
end # function
