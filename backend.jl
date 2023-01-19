using Statistics
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


# Define Quantum Force 
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

function vmcStep(R::Vector{Vector{Float64}}; s=1.0, α=0.16, β=0.5, κ=2.0)
    # Propse position
    Rnew = copy(R)
    Rnew[rand(1:2)] +=  s*(rand(Float64,3).-0.5)
    # Calc Acceptance
    ψold = trialWaveFunction(R; α, β, κ)
    ψpropose = trialWaveFunction(Rnew; α, β, κ)
    if rand() < ψpropose^2/ψold^2
        return  Rnew
    end
    return R
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


function initWalkers(;M=300)
    walkers = [[rand(Float64,3).-0.5,rand(Float64,3).-0.5] for _ in 1:M]
    return walkers
end





# Run the Algortihm
function runSimulation(;M=300,N=10000,n=1000,n_eq=0,s=0.1,α=0.15,β=0.5,κ=2.0, FP=false, Δτ=0.01, walkersReturn=false, p::Progress)

    # Initalize M heliumAtoms
    walkers = initWalkers(;M)

    # run equilibration if any
    if n_eq>0
        for _ in 1:n_eq
            if FP
               walkers = FPStep.(walkers; α, β, κ, Δτ)
            else
                walkers = vmcStep.(walkers; s, α, β, κ)
            end
            next!(p)
        end
        # when equilibrating I want mean over all other steps
        n=N-n_eq
    end
    
    # Initalize Energy for each walker:
    E = zeros(M)
    #Initalize return arrays
    returnEnergies = zeros(floor(Int,(N-n_eq)/n))
    returnStd = zeros(floor(Int,(N-n_eq)/n))
    
    # run main steps
    for i in 1:(N-n_eq)
        # update positions
        if FP
            walkers = FPStep.(walkers; α, β, κ, Δτ)
        else
            walkers = vmcStep.(walkers; s, α, β, κ)
        end

        # Update energy and variance mean over all walkers
        E_new  = localEnergy.(walkers; α, β, κ)
        
        E .+= E_new
        
        # push into return arrays 
        if i%n == 0
            # calc index
            j = floor(Int,i/n)

            # first calc mean over n steps, then over walkers
            E ./= n

            # push the right values
            returnEnergies[j] = mean(E)
            returnStd[j] = √(var(E)/(M-1))

            #reset E 
            E .= 0.0
        end

        next!(p)
    end

    if walkersReturn
        return walkers
    end

    return returnEnergies,returnStd

end