include("dmc.jl")
using ProgressMeter
using Plots

N = 30000
M = 300
s = 1.0
α = 0.18
β = 0.38
κ = 1.85
NEQ = 1
n = 1

walkersReturn = false
Δτ = 0.03

a = Δτ
E₀=-3.3

p = Progress(N)
Ets = runSimulation(; M, N, n, n_eq=NEQ, s, α, β, κ, Δτ, E₀, p, walkersReturn)

plt = plot(Ets)
# display(plt)
