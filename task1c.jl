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
αs = 0:0.01:0.5
β = 0.5
κ = 2.0
NEQ = 4000
n = N-NEQ


function calculateResults(varyparam)
    stds = zeros(length(varyparam))
    avEs = zeros(length(varyparam))
    
    p = Progress(length(varyparam)*N)
    Threads.@threads for (i,α) in collect(enumerate(varyparam))
        avE, stdE = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p)
        avEs[i] = avE[1]
        stds[i] = stdE[1]
    end

    return avEs, stds
end

# @time avEs, stds = calculateResults(αs)

@. model(x,p) = p[1]*x^(1)+p[2]*x^(2)+p[3]*x^(3)+p[4]*x^(4)+p[5]
p0 = zeros(5)

fitEnergies = curve_fit(model, αs, avEs, p0)
fitStds = curve_fit(model, αs, stds, p0)

fittedEv = model(αs, coef(fitEnergies))
fittedStds = model(αs, coef(fitStds))

p1 = scatter(αs, avEs.±stds, mc="darkblue", label="")
plot!(αs, fittedEv, lc="red", label="4th order fit with min at α=$(αs[findmin(fittedEv)[2]])")
title!("Energy dependence on α")
xlabel!(L"α")
ylabel!(L"E_L")

p2 = scatter(αs, stds, label="", mc="darkblue")
plot!(αs,fittedStds, color="red", label="4th order fit with min at α=$(αs[findmin(fittedStds)[2]])")
title!("Std. of Energies depending on α")
xlabel!(L"α")
ylabel!(L"σ_{E_L}")


savefig(p1, "saves/task1c.avEnergies.pdf")
savefig(p2, "saves/task1c.avStd.pdf")
display(p1)
display(p2)