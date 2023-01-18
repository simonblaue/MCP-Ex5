
using LinearAlgebra
using ProgressMeter

function trialWaveFunction(R::Vector{Float64}; α=0.16, β=0.5, κ=2.0)
    Re1 = R[1:3]
    Re2 = R[4:6]
    r₁ = norm(Re1)
    r₂ = norm(Re2)
    r₁₂ = norm(Re1-Re2)

    ψ1 = exp(-κ*r₁)
    ψ2 = exp(-κ*r₂)

    ψ = ψ1 * ψ2 * exp(β*r₁₂/(1+α*r₁₂))
    return ψ
end


function localEnergy(R::Vector{Float64}; α=0.16, β=0.5, κ=2.0)
    Re1 = R[1:3]
    Re2 = R[4:6]
    r₁ = norm(Re1)
    r₂ = norm(Re2)
    r₁₂ = norm(Re1-Re2)

    u = 1+α*r₁₂

    E1 = (κ-2)/r₁
    E2 = (κ-2)/r₂

    Eia = 1/r₁₂ * (1 - 2*β/u^2) + 2*β*α/u^3 - κ^2 - β^2/u^4 + κ*β/u^2* sum((Re1/r₁ - Re2/r₂).*((Re1-Re2)/r₁₂))

    E = E1 + E2 + Eia
    return E
end


# Define Quantum Force 
function quantumForce(R::Vector{Float64}; α=0.16, β=0.5, κ=2.0)
    Re1 = R[1:3]
    Re2 = R[4:6]
    r₁ = norm(Re1)
    r₂ = norm(Re2)
    r₁₂ = norm(Re1-Re2)

    u = 1+α*r₁₂

    F1 = 2 * (-κ/r₁ * Re1 + β*( 1/(r₁₂*u)  - α/u^2) * (Re1-Re2))
    F2 = 2 * (-κ/r₂ * Re2 - β*( 1/(r₁₂*u)  - α/u^2) * (Re1-Re2))

    return [F1...,F2...]
end


function GreensFunction(R1::Vector{Float64}, R2::Vector{Float64}; α=0.16, β=0.5, κ=2.0, Δτ=0.5)
    F = quantumForce(R1)
    exponent = -norm(R1-R2-(Δτ/2)*F)^2 / (2Δτ)

    pref = (2*π*Δτ)^-3 
    G = pref * exp(exponent)

    return G
end


# Define random step
function FPStep(R::Vector{Float64}; s=1.0, α=0.16, β=0.5, κ=2.0, Δτ=0.5)

    # propose a step ( copy old one first)
    proposeR = copy(R)
 
    randomStep = sqrt(Δτ/2)*(rand(Float64,6))
    F = quantumForce(R)
    proposeR += randomStep + F*Δτ/2

    # Calculate acceptance ratio VMC
    ψold = trialWaveFunction(R)
    ψpropose = trialWaveFunction(proposeR)
    acceptance = ψpropose^2/ψold^2

    # Correct acceptance ratio
    wR2R1 = GreensFunction(proposeR, R)
    wR1R2 = GreensFunction(R, proposeR)
    acceptance *= wR2R1/wR1R2
    acceptance = min(1, acceptance)

    # Take or leave the propsed new position
    if rand() < acceptance
         R = proposeR
    end
end

function branching(walkers::Vector{Vector{Float64}}, ET::Float64)
    # q <= 1 -> walker survives with prob q, death with 1-q
    # q >  1 -> birth to m new walkers m= floor(Int, q+r) with r∈[0,1] random

    maxM = 600
    minM = 150

    EL = localEnergy.(walkers)
    q = exp.(-Δτ*(EL .- ET))
    r = rand(length(q))

    survivers = walkers[ q .> r ]
    breeders = walkers[q .> 1]

    breeders_q = q[q .> 1]

    m = floor.(Int,breeders_q + rand(length(breeders)))

    for (i,b) in enumerate(breeders)
        for _ in 1:m[i]
            push!(survivers, copy(b))
        end
    end

    # M_new = length(survivers)

    # if M_new>maxM
    #     diff = M_new-maxM
    #     indizes = rand(1:length(survivers), M_new-diff)
    #     survivers = survivers[indizes]
    # else
    #     for _ in minM-M_new
    #         push!(survivers, rand(Float64, 6).-0.5)
    #     end
    # end

    return survivers
end

function runSimulation(;M=300,N=10000,n=1000,n_eq=0,s=1.0,α=0.16,β=0.5,κ=2.0, Δτ=0.5, E₀=-2.99, p::Progress, walkersReturn=false)
    # Initalize M heliumAtoms
    walkers = Vector{Vector{Float64}}(undef, M)
    for i in 1:M
        walkers[i] = rand(Float64, 6).-0.5
    end

    # Initalize Energy and variance for each walker:

    ET = E₀
      # run equilibration if any
    if n_eq>0
        for _ in 1:n_eq
            # Do FP step
            FPStep.(walkers)
            # Do Brancing
            walkers = branching(walkers, ET)
            ET = E₀ + a/Δτ*log(M,length(walkers))            
            next!(p)
        end
    end
    
    #Initalize return arrays
    returnEnergies = zeros(floor(Int,(N-n_eq)/n))
    # returnStd = zeros(floor(Int,(N-n_eq)/n))
    
    # run main steps
    for i in 1:(N-n_eq)
        if length(walkers) > 30000
            j = floor(Int,i/n)
            returnEnergies[j] = ET
            return returnEnergies[1:j]
        end
        # Do FP step
        FPStep.(walkers)
        # Do Brancing
        walkers = branching(walkers, ET)
        ET = E₀ + a/Δτ*log(M,length(walkers))            
        next!(p)

        # push into return arrays 
        if i%n == 0
            # calc index
            j = floor(Int,i/n)
            returnEnergies[j] = ET
        end
        next!(p)
    end

    if walkersReturn
        return walkers
    end

    return returnEnergies

end

