include("backend.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed

N = 10000
M = 300
s = 0.1
αs = [0., 0.1, 0.2, 0.3, 0.4, 0.5]
β = 0.5
κ = 2.0
NEQ = 0
n = 1000


function calculateResults(αs)
    stds = zeros(length(αs),10)
    avEs = zeros(length(αs),10)
    
    p = Progress(3*N)
    Threads.@threads for (i,α) in collect(enumerate(αs))
        avE, stdE = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p)
        avEs[i,:] = avE
        stds[i,:] = stdE
    end

    return avEs, stds
end

@time avEs, stds = calculateResults(αs)

labels = [L"\alpha = 0." L"\alpha = 0.1" L"\alpha = 0.2" L"\alpha = 0.3" L"\alpha = 0.4" L"\alpha = 0.5"]
p1 = plot(1:n:N-NEQ,avEs', label=labels)
vline!([4000], color="black", label=L"Equilbration seperation at $N=4000$")
xlabel!("Steps")
ylabel!("Mean Energy")
title!("Mean energies")

p2 = plot(1:n:N-NEQ,stds', label=labels)
xlabel!("Steps")
ylabel!("Mean Energy")
title!("Std error")

savefig(p1, "saves/task1b.energies.pdf")
savefig(p2, "saves/task1b.stds.pdf")

display(p1)
display(p2)