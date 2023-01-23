include("dmc.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit
using StatsBase

N = 30000
M = 300
α = 0.18
β = 0.5
κ = 2.0
NEQ = 10000
n = 100

walkersReturn = false
Δτ = 0.03

a = Δτ
E₀=-2.894

p = Progress(N)
# res = runSimulation(; M, N, n, n_eq=NEQ, α, β, κ, Δτ, E₀, p, walkersReturn)

# meanEnergy = fill(mean(res[1]), length(res[1]))
# meanNumWalkers = fill(mean(res[2]), length(res[1]))

# l = @layout [a; b]

# p1 = plot([res[1], meanEnergy],
#     label=["" "Mean Energy at $(round(meanEnergy[1], digits=2))"])
# p2 = plot([res[2], meanNumWalkers],
#     label=["" "Mean #Walkers: $(round(meanNumWalkers[1], digits=2))"])


# plt = plot(p1,p2, layout=l,
#     xlabel=["" "Steps"],
#     ylabel=[L"Energie $E_T$" "#Walkers"]
# )

# savefig(plt, "saves/task2b.Energies.pdf")
# display(plt)

# walkers = runSimulation(; M, N, n, n_eq=NEQ, α, β, κ, Δτ, E₀, p, walkersReturn=true)
function calculatePositions(walkers)
    pos_list = zeros(length(walkers)*2)
    for (i,w) in enumerate(walkers)
        pos_list[i] = norm(w[1]) 
        pos_list[i+length(walkers)] =  norm(w[2])
    end

    return pos_list
end

res = calculatePositions(walkers)

h = fit(Histogram, res, nbins=60)
h = normalize(h, mode=:pdf)

x = h.edges[1][1:end-1]
x_fit = LinRange(x[1], x[end], 100)

gaussian(x, μ, σ) = exp(-(x - μ)^2 / 2σ^2) / (√(2π) * σ)
@. model(x, p) = 0.5 * gaussian(x, p[1], p[2])

p0 = [1.0, 1.0]
fit_model = curve_fit(model, x, h.weights, p0)

plt = plot(h, label="")

(μ0, σ0) = coef(fit_model)
plot!(x_fit, model(x_fit, coef(fit_model)), label="Two Gaussian Fit: μ₁=$(round(μ0, digits=2)), σ₁=$(round(σ0, digits=2))")


savefig(plt, "saves/task2b.density.pdf")

display(plt)