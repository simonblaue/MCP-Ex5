include("vmc.jl")
using LaTeXStrings
using Plots
using Measurements
using ProgressBars

function calcElandStd()
    M = 300
    N = 10000
    n_equi = 2500
    s = 0.1
    κ = 2
    β = 1/2
    αs = 0:0.05:0.8

    avEnergies = []
    stdEnergies = []

    Threads.@threads for α in αs
        walkers = initWalkers(M)
        avEnergy, stdEnergy = vmc2(walkers, s, α, β, κ, N, n_equi)
        push!(avEnergies, avEnergy)
        push!(stdEnergies, stdEnergy)
    end

    return avEnergies, stdEnergies
end

# @time avEnergies, stdEnergies = calcElandStd()


p1 = plot( 0:0.05:0.8 ,avEnergies .± stdEnergies, title="Energy dependence on α", xlabel=L"α", ylabel=L"E_L", legend=false)

savefig(p1, "saves/figures/task1c.avEnergies.pdf")
# savefig(p2, "saves/figures/task1c.avStd.pdf")
display(p1)
