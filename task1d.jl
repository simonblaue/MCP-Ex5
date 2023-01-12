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
α = 0.16
β = 0.5
κs = 1.7:0.005:2.2
NEQ = 3000
n = N-NEQ


function calculateResults(varyparam)
    stds = zeros(length(varyparam))
    avEs = zeros(length(varyparam))
    
    p = Progress(length(varyparam)*N)
    Threads.@threads for (i,κ) in collect(enumerate(varyparam))
        avE, stdE = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p)
        avEs[i] = avE[1]
        stds[i] = stdE[1]
    end

    return avEs, stds
end

@time avEs, stds = calculateResults(κs)

@. model(x,p) = p[1]*x^(1)+p[2]*x^(2)+p[3]*x^(3)+p[4]
p0 = zeros(4)

fitEnergies = curve_fit(model, κs, avEs, p0)
fitStds = curve_fit(model, κs, stds, p0)

fittedEv = model(κs, coef(fitEnergies))
fittedStds = model(κs, coef(fitStds))

p1 = scatter(κs, avEs.±stds, mc="darkblue", label="")
plot!(κs, fittedEv, lc="red", label="3rd order fit with min at κ=$(κs[findmin(fittedEv)[2]])")
title!("Energy dependence on κ")
xlabel!(L"κ")
ylabel!(L"E_L")

p2 = scatter(κs, stds, label="", mc="darkblue")
plot!(κs,fittedStds, color="red", label="3rd order fit with min at κ=$(κs[findmin(fittedStds)[2]])")
title!("Std. of Energies depending on κ")
xlabel!(L"κ")
ylabel!(L"σ_{E_L}")


savefig(p1, "saves/task1d.avEnergies.pdf")
savefig(p2, "saves/task1d.avStd.pdf")
display(p1)
display(p2)