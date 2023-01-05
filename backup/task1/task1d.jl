include("vmc.jl")
using LaTeXStrings
using Plots
using Measurements
using ProgressBars
using CurveFit

function calcElandStd(κs)
    M = 300
    N = 10000
    n_equi = 3000
    s = 0.1
    α = 0.19
    β = 1/2

    avEnergies = zeros(length(κs))
    stdEnergies = zeros(length(κs))

    Threads.@threads for (i,κ) in collect(enumerate(κs))
        walkers = initWalkers(M)
        (avEnergies[i], stdEnergies[i]) = vmc2(walkers, s, α, β, κ, N, n_equi)
    end

    return avEnergies, stdEnergies
end


κs = 1.7:0.005:2.2
# @time avEnergies, stdEnergies = calcElandStd(κs)

fit = curve_fit(Polynomial, κs, avEnergies, 2)

fitκs =1.7:0.001:2.2
yb = fit.(fitκs)

p1 = plot(κs, avEnergies, ribbon=stdEnergies, lc="black", title="Energy dependence on κ", xlabel=L"κ", ylabel=L"E_L", label="")
scatter!(κs, avEnergies, label="")
plot!(fitκs,yb, color="red", label="Quadratic fit with min at κ=$(fitκs[findmin(yb)[2]])")
p2 = scatter(κs, stdEnergies, title="Std. of Energies depending on α", xlabel=L"κ", ylabel=L"σ_{E_L}", legend=false)

savefig(p1, "saves/figures/task1d.avEnergies.pdf")
savefig(p2, "saves/figures/task1d.avStd.pdf")
display(p2)
display(p1)

