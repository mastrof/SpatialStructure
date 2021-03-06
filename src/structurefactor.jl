export structurefactor

function density_ft(particles::T, q::AbstractVector{Float64}) where {S<:AbstractParticle,T<:AbstractVector{S}}
    sum([cos(dot(q, p.r)) + im * sin(dot(q, p.r)) for p in particles])
end # function

function density_ft(particles::T, q::AbstractVector{Float64}) where {S<:AbstractVector,T<:AbstractVector{S}}
    sum([cos(dot(q, r)) + im * sin(dot(q, r)) for r in particles])
end # function

function structurefactor(particles, q::AbstractVector{Float64})
    nparticles = length(particles)
    ρft = density_ft(particles, q)
    Sq = abs2(ρft) / nparticles
    return Sq
end # function

function structurefactor(particles, qs::T) where {S<:AbstractVector{Float64}, T<:AbstractVector{S}}
    Sq = 0.0
    for q in qs
        Sq += structurefactor(particles, q)
    end # for
    Sq /= length(qs)
    return Sq
end # function

function structurefactor(particles, qs::AbstractMatrix{Float64})
    Sq = 0.0
    for q in eachrow(qs)
        Sq += structurefactor(particles, q)
    end # for
    Sq /= size(qs, 1)
    return Sq
end # function

function structurefactor(particles, allqs::AbstractDict)
    nqs = length(keys(allqs))
    Sq = zeros(nqs)
    sortedqs = sort(collect(keys(allqs)))
    for (i, qabs) in enumerate(sortedqs)
        qs = allqs[qabs]
        Sq[i] = structurefactor(particles, qs)
    end # for
    return [sortedqs Sq]
end # function
