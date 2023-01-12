include("backend.jl")
using Plots
using LaTeXStrings
using ProgressMeter
using Optim

# TODO: variier über alle drei parameter gleichzeitig whoop whoop

N = 8000
M = 300
s = 1.0
NEQ = 2500
n = N-NEQ

optim_steps = 100
startparmas = [0.18,0.38,1.85]

p = Progress(N*277)
function simulation(params)
    α = params[1]
    β = params[2]
    κ = params[3]
    resEnergy, resStd = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ,p=p)

    return resEnergy[1]
end

result = optimize(simulation, startparmas, iterations=optim_steps)

display(result)

optimalParams = Optim.minimizer(result)
foundMinimum = Optim.minimum(result)

println("α = $(optimalParams[1]) \nβ = $(optimalParams[2]) \nκ = $(optimalParams[3])")
println("Minimal Energy: $foundMinimum ")

p = Progress(N)
res = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=optimalParams[1],β=optimalParams[2],κ=optimalParams[3],p=p)
println("$res")

res = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=0.18,β=0.38,κ=1.85,p=p)