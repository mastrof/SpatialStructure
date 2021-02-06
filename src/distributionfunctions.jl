export
    radialdistribution,
    angulardistribution,
    spatialdistribution

@doc raw"""
    radialdistribution(particles, box, npoints::Int)
    radialdistribution(particles, box, δr::Float64)

Evaluate the radial distribution function of `particles` in a periodic
box `box`.
`npoints` specifies the number of equally spaced points to sample from `box`, otherwise `δr` specifies the spacing of the sampled points.

---
    radialdistribution(particles_A, particles_B, box, npoints::Int)
    radialdistribution(particles_A, particles_B, box, δr::Float64)

Evaluate the radial distribution function between species `particles_A` and `particles_B` in a periodic box `box`.
The total radial distribution function of the system is obtained by summing over all the partial radial distributions of each pair of species, weighted by their partial concentrations (x_A = N_A/N):
    RDF(tot) = ∑ x_A * x_B * RDF(A,B)
"""
function radialdistribution(particles, box, npoints::Int)
    nparticles = size(particles, 1)
    ρ = (nparticles - 1) / prod(box)

    rpoints, hist = histradial(particles, box, npoints)
    n = hist ./ nparticles
    δr = rpoints[2] - rpoints[1]
    n_idealgas = [4π * ρ * rpoints[k]^2 * δr for k in 1:npoints]

    rdf = n ./ n_idealgas
    return rpoints, rdf
end # function

function radialdistribution(particles, box, δr::Float64)
    npoints = floor(Int, (minimum(box) / 2) ÷ δr) + 1
    radialdistribution(particles, box, npoints)
end # function

function radialdistribution(particles_A, particles_B, box, npoints::Int)
    nparticles_A = size(particles_A, 1)
    nparticles_B = size(particles_B, 1)
    ρA = nparticles_A / prod(box)
    ρB = nparticles_B / prod(box)
    ρ = sqrt(ρA * ρB) # partial density

    rpoints, hist = histradial(particles_A, particles_B, box, npoints)
    n = hist ./ sqrt(nparticles_A * nparticles_B)
    δr = rpoints[2] - rpoints[1]
    n_idealgas = [4π * ρ * rpoints[k]^2 * δr for k in 1:npoints]

    rdf = n ./ n_idealgas
    return rpoints, rdf
end # function

function radialdistribution(particles_A, particles_B, box, δr::Float64)
    npoints = floor(Int, (minimum(box) / 2) ÷ δr) + 1
    radialdistribution(particles_A, particles_B, box, npoints)
end # function

