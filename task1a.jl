include("backend.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed

N = 30000
M = 300
α = 0.15
β = 0.5
κ = 2.0
NEQ = 0
n = 1000

s_list = [0.1, 1.1, 10.0]

function calculateResults(s_list)
    stds = zeros(3,30)
    avEs = zeros(3,30)
    
    p = Progress(3*N)
    Threads.@threads for (i,s) in collect(enumerate(s_list))
        avE, stdE = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p)
        avEs[i,:] = avE
        stds[i,:] = stdE
    end

    return avEs, stds
end

# @time avEs, stds = calculateResults(s_list)

display(stds[1][end])
display(stds[2][end])
display(stds[3][end])

labels = [L"s=0.1" L"s=1.0" L"s=10.0"]
p1 = plot(1:n:N-NEQ,avEs', label=labels)
xlabel!("Steps")
ylabel!(L"$\langle E_L^{1000} \rangle$")
title!("Mean energies")

p2 = plot(1:n:N-NEQ,stds', label=labels)
xlabel!("Steps")
ylabel!(L"$σ(\langle E_L^{1000} \rangle )$")
title!("Std error")

savefig(p1, "saves/task1a.energies.pdf")
savefig(p2, "saves/task1a.stds.pdf")

display(p1)
display(p2)