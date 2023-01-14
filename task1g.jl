include("backend.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit

N = 10000
M = 300
s = 1.0
α = 0.18
β = 0.38
κ = 1.85
NEQ = 2000
n = N-NEQ

FP = true
Δτs = [0.01, 0.05, 0.1, 0.2, 1.0]

function calculateResults(varyparam)
    stds = zeros(length(varyparam))
    avEs = zeros(length(varyparam))
    
    p = Progress(length(varyparam)*N)
    Threads.@threads for (i,Δτ) in collect(enumerate(varyparam))
        avE, stdE = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p, FP=FP, Δτ=Δτ)
        avEs[i] = avE[1]
        stds[i] = stdE[1]
    end

    return avEs, stds
end


@time avEs, stds = calculateResults(Δτs)


p1 = scatter(Δτs, avEs.±stds, label="", mc="darkblue")
xaxis!(:log)
xticks!(Δτs, string.(Δτs))
title!("Energy dependence on Δτ")
xlabel!(L"Δτ")
ylabel!(L"\bar{E}_L")

p2 = scatter(Δτs, stds, mc="darkblue", label="")
xaxis!(:log)
xticks!(Δτs, string.(Δτs))
title!("Std dependence on Δτ")
xlabel!(L"Δτ")
ylabel!(L"\bar{E}_L")

savefig(p1, "saves/task1g.avEnergies.pdf")
savefig(p2, "saves/task1g.avStd.pdf")
display(p1)
display(p2)