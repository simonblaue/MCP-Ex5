using LinearAlgebra
using Statistics

const M = 300
const β = 1/2
const κ = 2
const N = 30000
const n = 1000
const n_equi = 2000

mutable struct Walker
    posE1::Vector{Float64}
    posE2::Vector{Float64}
end

function localEnergy(walker::Walker,α, β, κ)
    r₁ = norm(walker.posE1)
    r₂ = norm(walker.posE2)
    r₁₂ = norm(walker.posE1-walker.posE2)
    u = 1 + α*r₁₂
    return (κ-2)/r₁ + (κ-2)/r₂ + 1/r₁₂ * (1 - 2*β/u^2) + 2*β*α/u^3 - κ^2 - β^2/u^4 + κ*β/u^2 * sum((walker.posE1/r₁ - walker.posE2/r₂) .* (walker.posE1-walker.posE2)/r₁₂)
end

function wavefunction(walker::Walker, α, β, κ)
    r₁ = norm(walker.posE1)
    r₂ = norm(walker.posE2)
    r₁₂ = norm(walker.posE1-walker.posE2)
    return ℯ^(-κ*r₁) * ℯ^(-κ*r₂) * ℯ^(β*r₁₂/(1+α*r₁₂))
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


function randomStep(w::Walker, s, α)
    # propose random step
    randomStep = s.*(rand(Float64,3).-0.5)
    if rand()<0.5
        wpropose = Walker(w.posE1, w.posE2)
        wpropose.posE1 += randomStep
    else
        wpropose = Walker(w.posE1, w.posE2)
        wpropose.posE2 += randomStep
    end
    # acceptance
    wfold = wavefunction(w, α, β, κ)
    wfpropose = wavefunction(wpropose, α, β, κ)
    if rand() < wfpropose^2/wfold^2
        w.posE1 = wpropose.posE1
        w.posE2 = wpropose.posE2
    end
end

function vmc(walkers, s, α,β,κ, N, n)
    # mesurement as average over every n steps
    averagingsteps = floor(Int,N/n)

    n_energies = Array{Float64}(undef,(M,n))
    av_e = Array{Float64}(undef,(M,averagingsteps))

    for i in 1:averagingsteps
        for j in 1:n
            n_energies[:,j] = localEnergy.(walkers, α, β, κ)
            randomStep.(walkers, s, α)
        end

        av_e[:,i] = mean(n_energies, dims=2)
    end

    walkerMean = mean(av_e, dims=1)'
    walkerStd = std(av_e, dims=1)'
    
    return walkerMean, walkerStd
end

function vmc2(walkers, s, α, β, κ, N, n_equi)
    # measurment as mean after equilibration time

    av_e = Array{Float64}(undef,N-n_equi)
    std_e = Array{Float64}(undef,N-n_equi)

    # Equilibration
    for i in 1:n_equi
        randomStep.(walkers, s, α)
    end

    # Messurements
    for i in 1:N-n_equi
            m_energies = localEnergy.(walkers, α, β, κ)
            # mean over all walkers
            av_e[i] = mean(m_energies)
            std_e[i] = std(m_energies)
            randomStep.(walkers, s, α)
        end

    # mean over all already equilibrated steps
    # TODO: Is the std correctly calculated?
    return mean(av_e), mean(std_e)
end