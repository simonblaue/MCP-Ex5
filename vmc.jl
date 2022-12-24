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


function vmc(walkers, s, N, n)
    α = 0.15
    energies = Array{Float64}(undef,(M,N))
    energies_n = []
    mean_energies_n = []
    std_en_n = []
    for i in 1:N
        energies[:,i] = localEnergy.(walkers, α)
        # if i%n == 0
        #     push!(energies_n,mean(energies, dims=1))
        # end

        #random step for every walker
        for w in walkers
            randomStep = rand(Float64,3).*s.-s/2
            if rand()<0.5
                w.posE1 += randomStep
            else
                w.posE2 += randomStep
            end
        end
    end

    # walker_av_energy = [mean(e) for e in energies_n]
    # wlaker_std_energy = [std(e) for e in energies_n]


    return energies
end





walkers = initWalkers(300)
sList = [0.1, 1.0, 10.0]
s = sList[2]

@time all_energies = vmc(walkers,s,N,n)



# av_e = [mean(e) for e in energies_n]
# display(plot(av_e))
# display(plot(std_e))