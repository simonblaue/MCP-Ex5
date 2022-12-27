using LinearAlgebra
using Statistics
using Plots

const M = 300
const β = 1/2
const κ = 2
const N = 30000
const n = 1000

mutable struct Walker
    posE1::Vector{Float64}
    posE2::Vector{Float64}
end

function localEnergy(walker::Walker,α)
    r₁ = walker.posE1
    r₂ = walker.posE2
    r₁₂ = norm(r₁-r₂)
    u = 1 + α*r₁₂
    return (κ-2)/norm(r₁) + (κ-2)/norm(r₂) + 1/r₁₂ * (1 - 2β/u^2) + 2*β*α/u^3 -κ^2 - β^2/u^4 + κ*β/u^2 * norm((r₁/norm(r₁) - r₂/norm(r₂)) .* (r₁-r₂)/r₁₂)
end

function initWalkers(M)
    walkers = Array{Walker}(undef, M)
    for i in 1:M
        pos1 = rand(Float64,3)
        pos2 = rand(Float64,3)
        walkers[i] = Walker(pos1,pos2)
    end
    return walkers
end


function randomStep(w::Walker, s)
    randomStep = rand(Float64,3).*s.-s/2
    if rand()<0.5
        w.posE1 += randomStep
    else
        w.posE2 += randomStep
    end
end

function vmc(walkers, s, α, N, n)
    
    averagingsteps = floor(Int,N/n)

    n_energies = Array{Float64}(undef,(M,n))
    av_e = Array{Float64}(undef,(M,averagingsteps))

    for i in 1:averagingsteps
        for j in 1:n
            n_energies[:,j] = localEnergy.(walkers, α)
            randomStep.(walkers, s)
        end

        av_e[:,i] = mean(n_energies, dims=2)
    end

    walkerMean = mean(av_e, dims=1)'
    walkerStd = std(av_e, dims=1)'
    
    return walkerMean, walkerStd
end
