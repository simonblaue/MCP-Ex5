using Statistics
using LinearAlgebra
using ProgressMeter


# define a two position struct
mutable struct heliumAtom
    r₁::Vector{Float64}
    r₂::Vector{Float64}
end

# Define the trial wave function
function trialWaveFunction(atom::heliumAtom, α::Float64, β::Float64, κ::Float64)
    # calculate distances
    r₁ = norm(atom.r₁)
    r₂ = norm(atom.r₂)
    r₁₂ = norm(atom.r₁-atom.r₂)

    #calculate single electron terms
    ψ1 = exp(-κ*r₁)
    ψ2 = exp(-κ*r₂)

    #calculate value of trial wave function
    ψ = ψ1 * ψ2 * exp(β*r₁₂/(1+α*r₁₂))

    return ψ
end

# Define the local Energy
function localEnergy(atom::heliumAtom, α::Float64, β::Float64, κ::Float64)
    # calculate distances
    r₁ = norm(atom.r₁)
    r₂ = norm(atom.r₂)
    r₁₂ = norm(atom.r₁-atom.r₂)

    # calculate u
    u = 1+α*r₁₂

    # calculate single electron parts
    E1 = (κ-2)/r₁
    E2 = (κ-2)/r₂

    #calculate interaction terms
    Eia = 1/r₁₂ * (1 - 2*β/u^2) + 2*β*α/u^3 - κ^2 - β^2/u^4 + κ*β/u^2* sum((atom.r₁/r₁ - atom.r₂/r₂).*((atom.r₁-atom.r₂)/r₁₂))

    #Local Energy
    E = E1 + E2 + Eia

    return E
end

# Define Quantum Force 
function quantumForce(atom::heliumAtom, α::Float64, β::Float64, κ::Float64)
    # calculate distances
    r₁ = norm(atom.r₁)
    r₂ = norm(atom.r₂)
    r₁₂ = norm(atom.r₁-atom.r₂)

    F1 = -κ/r₁ * atom.r₁ + β/(r₁₂*(1+α*r₁₂)) * (atom.r₁-atom.r₂) - α*β/(1+α*r₁₂)^2 * (atom.r₁-atom.r₂)
    F2 = -κ/r₂ * atom.r₂ - β/(r₁₂*(1+α*r₁₂)) * (atom.r₁-atom.r₂) + α*β/(1+α*r₁₂)^2 * (atom.r₁-atom.r₂)

    return [F1,F2]
end


function GreensFunction(atom1::heliumAtom, atom2::heliumAtom, Δτ::Float64)
    pref = 1/(2*π*Δτ)^3 
    # putting elctron positions in handable 6 vecs
    R1 = Vector{Float64}(undef, 6)
    R1[1:3] = atom1.r₁
    R1[4:6] = atom1.r₂
    R2 = Vector{Float64}(undef, 6)
    R2[1:3] = atom1.r₁
    R2[4:6] = atom1.r₂

    # Calcute exponent
    F = Vector{Float64}(undef, 6)
    qForce = quantumForce(atom1, α, β , κ)
    F[1:3] = qForce[1]
    F[4:6] = qForce[2]

    exponent = -sum((R2-R1-(Δτ/2).*F).^2)/(2*Δτ)

    # Calculate Greens Function Value
    G = pref * exp(exponent)

    return G
end

# Define random step
function randomStep(atom::heliumAtom, s::Float64, α::Float64, β::Float64, κ::Float64; FP=false, Δτ=0.5)

    # propose a step ( copy old one first)
    proposeAtom = heliumAtom(atom.r₁, atom.r₂)

    # Normal VMC
    if !FP
        randomStep = s.*(rand(Float64,3).-0.5)
        # choose one of the electrons randomly
        if rand()<0.5
            proposeAtom = heliumAtom(atom.r₁, atom.r₂)
            proposeAtom.r₁ += randomStep
        else
            proposeAtom = heliumAtom(atom.r₁, atom.r₂)
            proposeAtom.r₂ += randomStep
        end
    
    # Do Fokker Plank if nessecary
    else
        randomStep = sqrt(Δτ/2).*(rand(Float64,6))
        F = quantumForce(atom, α, β, κ)
        proposeAtom.r₁ += randomStep[1:3] + F[1]*Δτ/2
        proposeAtom.r₂ += randomStep[4:6] + F[2]*Δτ/2
    end

    # Calculate acceptance ratio VMC
    ψold = trialWaveFunction(atom, α, β, κ)
    ψpropose = trialWaveFunction(proposeAtom, α, β, κ)
    acceptance = ψpropose^2/ψold^2

    # Correct acceptance ratio for VMC-FP
    if FP
        wR2R1 = GreensFunction(proposeAtom, atom, Δτ)
        wR1R2 = GreensFunction(atom, proposeAtom, Δτ)
        acceptance *= wR2R1/wR1R2
        acceptance = min(1, acceptance)
    end

    # Take or leave the propsed new position
    if rand() < acceptance
         atom.r₁ = proposeAtom.r₁
         atom.r₂ = proposeAtom.r₂
    end
end

# Run the VMC(-FP) Algortihm
function runSimulation(;M=300,N=10000,n=1000,n_eq=0,s=0.1,α=0.15,β=0.5,κ=2.0, FP=false, Δτ=0.5, p::Progress)

    # Initalize M heliumAtoms
    walkers = Array{heliumAtom}(undef, M)
    for i in 1:M
        pos1 = rand(Float64,3).-0.5
        pos2 = rand(Float64,3).-0.5
        walkers[i] = heliumAtom(pos1,pos2)
    end

    # Initalize Energy and variance for each walker:
    E = zeros(M)
    # var =zeros(M)

    # run equilibration if any
    if n_eq>0
        for _ in 1:n_eq
            if FP
                randomStep.(walkers, s, α, β, κ, FP=true, Δτ=Δτ)
            else
                randomStep.(walkers, s, α, β, κ)
            end
            next!(p)
        end
        # when equilibrating I want mean over all steps
        n=N-n_eq
    end
    
    #Initalize return arrays
    returnEnergies = zeros(floor(Int,(N-n_eq)/n))
    returnStd = zeros(floor(Int,(N-n_eq)/n))
    
    # run main steps
    for i in 1:(N-n_eq)
        # update positions
        if FP
            randomStep.(walkers, s, α, β, κ, FP=true, Δτ=Δτ)
        else
            randomStep.(walkers, s, α, β, κ)
        end

        # Update energy and variance mean over all walkers
        E_new  = localEnergy.(walkers, α, β, κ)

        E .+= E_new

        # push into return arrays 
        if i%n == 0
            # calc index
            j = floor(Int,i/n)

            # first calc mean over i steps, then over walkers
            E ./= n

            # push the right values
            returnEnergies[j] = mean(E)
            returnStd[j] = √(var(E)/(M-1))

            #reset E 
            E .= 0.0
        end

        next!(p)
    end

    return returnEnergies,returnStd

end

