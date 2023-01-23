include("dmc.jl")
using ProgressMeter
using Plots
using LaTeXStrings
using Statistics


N = 30000
M = 300
α = 0.18
β = 0.5
κ = 2.0
NEQ = 10000
n = 100

walkersReturn = false
Δτ = 0.03

a = Δτ
E₀=-2.894

p = Progress(N)
# res = runSimulation(; M, N, n, n_eq=NEQ, α, β, κ, Δτ, E₀, p, walkersReturn)

meanEnergy = fill(mean(res[1]), length(res[1]))
meanNumWalkers = fill(mean(res[2]), length(res[1]))

l = @layout [a; b]

p1 = plot([res[1], meanEnergy],
    label=["" "Mean Energy at $(round(meanEnergy[1], digits=2))"])
p2 = plot([res[2], meanNumWalkers],
    label=["" "Mean #Walkers: $(round(meanNumWalkers[1], digits=2))"])


plt = plot(p1,p2, layout=l,
    xlabel=["" "Steps"],
    ylabel=[L"Energie $E_T$" "#Walkers"]
)

savefig(plt, "saves/task2b.Energies.pdf")
display(plt)
