include("vmc.jl")
using LaTeXStrings
using Plots
using Measurements
using ProgressBars
using CurveFit

function calcElandStd(αs)
    M = 300
    N = 10000
    n_equi = 3000
    s = 0.1
    κ = 2
    β = 1/2

    avEnergies = zeros(length(αs))
    stdEnergies = zeros(length(αs))

    Threads.@threads for (i,α) in collect(enumerate(αs))
        walkers = initWalkers(M)
        # avEnergy, stdEnergy = vmc2(walkers, s, α, β, κ, N, n_equi)
        (avEnergies[i], stdEnergies[i]) = vmc2(walkers, s, α, β, κ, N, n_equi)
        # push!(stdEnergies, stdEnergy)
    end

    return avEnergies, stdEnergies
end


αs = 0:0.005:0.5
# @time avEnergies, stdEnergies = calcElandStd(αs)

fit = curve_fit(Polynomial, αs, avEnergies, 2)

fitα = 0:0.001:0.5
yb = fit.(fitα)

p1 = scatter(αs, avEnergies, title="Energy dependence on α", xlabel=L"α", ylabel=L"E_L", label="")
plot!(fitα,yb, color="red", label="Quadratic fit with min at α=$(fitα[findmin(yb)[2]])")
p2 = scatter(αs, stdEnergies, title="Std. of Energies depending on α", xlabel=L"α", ylabel=L"σ_{E_L}", legend=false)

savefig(p1, "saves/figures/task1c.avEnergies.pdf")
savefig(p2, "saves/figures/task1c.avStd.pdf")
display(p2)
display(p1)

