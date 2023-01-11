include("backend.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit

N = 10000
M = 300
s = 0.1
α = 0.18
β = 0.38
κ = 1.85
NEQ = 4000
n = N-NEQ

p = Progress(N)
avE, std = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ,p=p)

print("Minimal Energy: $(avE[1]) ± $(std[1])")