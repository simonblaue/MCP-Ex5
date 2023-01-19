include("dmc.jl")
using ProgressMeter
using Plots
using LaTeXStrings

N = 30000
M = 300
α = 0.18
β = 0.5
κ = 2.0
NEQ = 10000
n = 1

walkersReturn = false
Δτ = 0.03

a = Δτ
E₀=-2.894

p = Progress(N)
# res = runSimulation(; M, N, n, n_eq=NEQ, α, β, κ, Δτ, E₀, p, walkersReturn)

plt = plot([res[1],res[2]],
    layout=(2,1), 
    label="", 
    xlabel=["" "Stpes"],
    ylabel=[L"Energie $E_T$" "# Walkers"]
    )


savefig(plt, "saves/task2b.Energies.pdf")
display(plt)
