
using LinearAlgebra
using ProgressMeter

function trialWaveFunction(R::Vector{Vector{Float64}}; α=0.16, β=0.5, κ=2.0)
    r₁ = norm(R[1])
    r₂ = norm(R[2])
    r₁₂ = norm(R[1]-R[2])
    ψ1 = exp(-κ*r₁)
    ψ2 = exp(-κ*r₂)
    ψ = ψ1 * ψ2 * exp(β*r₁₂/(1+α*r₁₂))
    return ψ
end


function localEnergy(R::Vector{Vector{Float64}}; α=0.16, β=0.5, κ=2.0)
    r₁ = norm(R[1])
    r₂ = norm(R[2])
    r₁₂ = norm(R[1]-R[2])
    u = 1+α*r₁₂
    E1 = (κ-2)/r₁
    E2 = (κ-2)/r₂
    Eia = 1/r₁₂ * (1 - 2*β/u^2) + 2*β*α/u^3 - κ^2 - β^2/u^4 + κ*β/u^2* sum((R[1]/r₁ - R[2]/r₂).*((R[1]-R[2])/r₁₂))
    E = E1 + E2 + Eia
    return E
end


function quantumForce(R::Vector{Vector{Float64}}; α=0.16, β=0.5, κ=2.0)
    r₁ = norm(R[1])
    r₂ = norm(R[2])
    r₁₂ = norm(R[1]-R[2])
    u = 1+α*r₁₂
    F1 = 2 * (-κ/r₁ * R[1] + β*( 1/(r₁₂*u)  - α/u^2) * (R[1]-R[2]))
    F2 = 2 * (-κ/r₂ * R[2] - β*( 1/(r₁₂*u)  - α/u^2) * (R[1]-R[2]))

    return [F1,F2]
end


function GreensFunction(R::Vector{Vector{Float64}}, Rnew::Vector{Vector{Float64}}; α=0.16, β=0.5, κ=2.0, Δτ=0.01)
    F = quantumForce(R)
    exponent = -norm(Rnew-R-(Δτ/2)*F)^2 / (2Δτ)
    pref = (2*π*Δτ)^-3 
    G = pref * exp(exponent)
    return G
end


function FPStep(R::Vector{Vector{Float64}};α=0.16, β=0.5, κ=2.0, Δτ=0.01)
    Rnew = copy(R)
    randomStep = sqrt(Δτ)*[randn(Float64,3),randn(Float64,3)]
    F = quantumForce(Rnew; α, β, κ)
    Rnew += randomStep + F*Δτ/2

    ψold = trialWaveFunction(R; α, β, κ)
    ψpropose = trialWaveFunction(Rnew; α, β, κ)
    acceptance = ψpropose^2/ψold^2
    wRnewR = GreensFunction(Rnew, R; α, β, κ, Δτ)
    wRRnew = GreensFunction(R, Rnew; α, β, κ, Δτ)
    acceptance *= wRnewR/wRRnew
    if rand() < acceptance
        return Rnew
    end
    return R
end

function branching(walkers::Vector{Vector{Vector{Float64}}}, ET::Float64, acc_idx::BitVector)
    # q <= 1 -> walker survives with prob q, death with 1-q
    # q >  1 -> birth to m new walkers m= floor(Int, q+r) with r∈[0,1] random

    walkers_iq = walkers[acc_idx]

    EL = localEnergy.(walkers_iq)
    q = exp.(-Δτ*(EL .- ET))
    r = rand(length(q))

    survivers = walkers_iq[ q .> r ]
    breeders = walkers_iq[q .> 1]

    breeders_q = q[q .> 1]

    m = floor.(breeders_q + rand(length(breeders)))

    for (i,b) in enumerate(breeders)
        for _ in 1:m[i]
            push!(survivers, copy(b))
        end
    end
    return append!(survivers, walkers[.!acc_idx])
end


function initWalkers(;M=300)
    walkers = [[rand(Float64,3).-0.5,rand(Float64,3).-0.5] for _ in 1:M]
    return walkers
end

function runSimulation(;M=300,N=10000,n=1000,n_eq=0,α=0.16,β=0.5,κ=2.0, Δτ=0.5, E₀=-2.99, p::Progress, walkersReturn=false)
    # Initalize M heliumAtoms
    walkers = initWalkers(;M)

    # Initalize Energy
    ET = E₀
    # run equilibration if any
    if n_eq>0
        for _ in 1:n_eq
            # Do FP step
            newwalkers = FPStep.(walkers; α, β, κ, Δτ)
            acc_idx = walkers .!= newwalkers
            # Do Brancing
            walkers = branching(newwalkers, ET, acc_idx)
            ET = E₀ + a/Δτ*log(M/length(walkers))            
            next!(p)
        end
    end
    
    #Initalize return arrays
    returnEnergies = zeros(floor(Int,(N-n_eq)/n))
    returnNumWalkers = zeros(floor(Int,(N-n_eq)/n))
    
    # run main steps
    for i in 1:(N-n_eq)
        # Do FP step
        newwalkers = FPStep.(walkers; α, β, κ, Δτ)
        acc_idx = walkers .!= newwalkers
        # Do Brancing
        walkers = branching(newwalkers, ET, acc_idx)
        ET = E₀ + a/Δτ*log(M/length(walkers))            
        next!(p)

        # push into return arrays 
        if i%n == 0
            # calc index
            j = floor(Int,i/n)
            returnEnergies[j] = ET
            returnNumWalkers[j] = length(walkers)
        end
        next!(p)
    end

    if walkersReturn
        return walkers
    end

    return returnEnergies, returnNumWalkers

end

