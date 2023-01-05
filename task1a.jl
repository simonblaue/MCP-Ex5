include("vmc.jl")
using LaTeXStrings
using Plots

function calcAveandStds()
    M = 300
    N = 30000
    n = 1000

    α = 0.15
    β = 1/2
    κ = 2.

    sList = [0.1,1.0,10.0]

    avEnergies = []
    stdEnergies = []

    for s in sList
        walkers = initWalkers(M)
        avEnergy, stdEnergy = vmc(walkers, s, α,β,κ, N, n)
        push!(avEnergies, avEnergy)
        push!(stdEnergies, stdEnergy)
    end

    return avEnergies, stdEnergies
end


# @time avEnergies, stdEnergies = calcAveandStds()


labels = [L"s=1" L"s=0.1" L"s=10"]
p1 = plot(1:30, avEnergies, labels=labels, title="Average Energie")
p2 = plot(1:30, stdEnergies, labels=labels, title="Standard error Energie")

savefig(p1, "saves/figures/task1a.avEnergies.pdf")
savefig(p2, "saves/figures/task1a.avStd.pdf")
display(p1)
display(p2)

