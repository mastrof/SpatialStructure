export
    histradial,
    histangular,
    histspatial

AbstractCoordinates{T} = Union{AbstractVector{T}, AbstractParticle{D,T} where D}

getdimension(r::SVector{D,T}) where {D,T} = D

@doc raw"""
    histradial(particles::T, box, npoints::Int) where {S<:Real,T<:AbstractVector{S}}
    histradial(particles::T, box, npoints::Int) where {S<:AbstractParticle,T<:AbstractVector{S}}
    histradial(distances::T, box, npoints::Int) where {S<:Real,T<:AbstractMatrix{S}}

Evaluates the distribution of pair distances of `particles` in a box of edges `box` up to a maximum distance equal to half the shortest edge in order to avoid artefacts from periodic boundary conditions(assumed cubic).
`npoints` is the sampling density. Works with systems of arbitrary dimensionality.

---
    histradial(particles_A::Vector{T}, particles_B::Vector{T}, box, npoints::Int) where {S<:Real,T<:AbstractVector{S}}
    histradial(particles_A::Vector{T}, particles_B::Vector{T}, box, npoints::Int) where {S<:AbstractParticle,T<:AbstractVector{S}}

Evaluates the distribution of pair distances between `particles_A` and `particles_B` in a box of edges `box` up to a maximum distance equal to half the shortest edge in order to avoid artefacts from periodic boundary conditions(assumed cubic).
`npoints` is the sampling density.
Works with systems of arbitrary dimensionality.

"""
function histradial(particles::T, box, npoints::Int) where {S<:AbstractCoordinates,T<:AbstractVector{S}}
    nparticles = length(particles)
    bin_edges = collect(range(0, minimum(box/2); length = npoints + 1))
    δr = bin_edges[2] - bin_edges[1]
    rpoints = bin_edges[1:end-1] .+ δr / 2
    hist = zeros(npoints)

    @inbounds for i in 1:nparticles
        r_i = particles[i]
        for j in i+1:nparticles
            r_j = particles[j]
            r_ij = distance(r_i, r_j, box)
            k = floor(Int, r_ij / δr) + 1
            if k ≤ npoints
                hist[k] += 2.0
            end # if
        end # for
    end # for

    return rpoints, hist
end # function

function histradial(distances::T, box, npoints) where {S<:Real,T<:AbstractMatrix}
    nparticles = size(distances, 1)
    bin_edges = collect(range(0, minimum(box/2); length = npoints + 1))
    δr = bin_edges[2] - bin_edges[1]
    rpoints = bin_edges[1:end-1] .+ δr / 2
    hist = zeros(npoints)

    @inbounds for j in 1:nparticles
        for i in j+1:nparticles
            r_ij = distances[i,j]
            k = floor(Int, r_ij / δr) + 1
            if 1 ≤ k ≤ npoints
                hist[k] += 2.0
            end # if
        end # for
    end # for

    return rpoints, hist
end # function

function histradial(particles, box, δr::Float64)
    npoints = floor(Int, (minimum(box) / 2) ÷ δr) + 1
    histradial(particles, box, npoints)
end # function

function histradial(particles_A::T, particles_B::T, box, npoints::Int) where {S<:AbstractCoordinates,T<:AbstractVector{S}}
    nparticles_A = length(particles_A)
    nparticles_B = length(particles_B)
    bin_edges = collect(range(0, minimum(box/2); length = npoints + 1))
    δr = bin_edges[2] - bin_edges[1]
    rpoints = bin_edges[1:end-1] .+ δr / 2
    hist = zeros(npoints)

    @inbounds for i in 1:nparticles_A
        r_i = particles_A[i]
        for j in 1:nparticles_B
            r_j = particles_B[j]
            r_ij = distance(r_i, r_j, box)
            k = floor(Int, r_ij / δr) + 1
            if k ≤ npoints
                hist[k] += 1
            end # if
        end # for
    end # for

    return rpoints, hist
end # function

function histradial(particles_A, particles_B, box, δr::Float64)
    npoints = floor(Int, (minimum(box) / 2) ÷ δr) + 1
    histradial(particles_A, particles_B, box, npoints)
end # function






#======= NEEDS REVIEW =======#

@doc raw"""
    histangular(particles::Vector{T}, box::AbstractVector, anglerule::Function, binedges::AbstractVector) where {T<:AbstractParticle}

Evaluates the distribution of angles between pairs of `particles` in a box of edges `box`.
`anglerule` defines how the angle should be evaluated; it should take two AbstractParticle objects and `box` as inputs.
The binning is defined by `binedges`.
"""
function histangular(particles::Vector{T},
                     box::AbstractVector,
                     anglerule::Function,
                     binedges::AbstractVector) where {T<:AbstractParticle}
    nparticles = size(particles, 1)
    θpoints = binedges[1:end-1] .+ δθ / 2
    θs = zeros(nparticles * (nparticles - 1) ÷ 2)
    k = 1

    @inbounds for i in 1:nparticles
        r_i = particles[i]
        for j in i+1:nparticles
            r_j = particles[j]
            θs[k] = anglerule(r_i, r_j, box)
            k += 1
        end # for
    end # for

    hist = fit(Histogram, θs, binedges).weights
    return θpoints, hist
end # function




function histspatial(particles::Vector{T},
                     box::AbstractVector,
                     redges::AbstractVector,
                     θedges::AbstractVector,
                     ϕedges::AbstractVector,
                     anglerule) where {T<:AbstractParticle}
    nparticles = size(particles, 1)
    space = zeros(3, nparticles * (nparticles - 1) ÷ 2)
    k = 1
    
    @inbounds for i in 1:nparticles
        r_i = particles[i]
        for j in i+1:nparticles
            r_j = particles[j]
            rij = distance(r_i, r_j, box)
            θij, ϕij = anglerule(r_i, r_j, box)
            space[:, k] .= (rij, θij, ϕij)
            k += 1
        end # for
    end # for

    rs = space[1, :]
    θs = space[2, :]
    ϕs = space[3, :]
    hist = fit(Histogram, (rs, θs, ϕs), (redges, θedges, ϕedges))
    rpoints, θpoints, ϕpoints = midpoints.(hist.edges)
    
    return rpoints, θpoints, ϕpoints, hist.weights
end # function
