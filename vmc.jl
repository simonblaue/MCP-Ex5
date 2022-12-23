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



@time walkers = initWalkers(300)

walkers
@time localEnergy.(walkers, 0.14)