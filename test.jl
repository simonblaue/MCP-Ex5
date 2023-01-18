include("backend.jl")
using ProgressMeter
using Plots

# N = 10
# M = 300
# s = 1.0
# α = 0.18
# β = 0.38
# κ = 1.85
# NEQ = 10
# n = N-NEQ


# Δτ = 0.03

# a = Δτ
# E₀=-2.99


# Ets = runSimulation(p=p, E₀=-3.5)

# plot(Ets)


res = runSimulation(N=4, n=4);