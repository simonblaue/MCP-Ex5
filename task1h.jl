include("backend.jl")
using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit
using StatsBase

N = 4000
M = 300
s = 1.0
α = 0.18
β = 0.38
κ = 1.85
NEQ = 2000
n = N-NEQ

FP = true
Δτ = 0.1

p = Progress(N)
@time walkers = runSimulation(M=M,N=N,n=n,n_eq=NEQ,s=s,α=α,β=β,κ=κ, p=p, FP=FP, Δτ=Δτ, walkersReturn=true)


function calculatePositions(walkers)
    pos_list = zeros(length(walkers)*2)
    for (i,w) in enumerate(walkers)
        pos_list[i] = norm(w[1]) 
        pos_list[i+length(walkers)] =  norm(w[2])
    end

    return pos_list
end

# function calculateProbDistri(walkers)
#     probs = zeros(length(walkers))
#     for (i,w) in enumerate(walkers)
#         probs[i] = trialWaveFunction(w, α, β, κ)^2
#     end

#     return probs
# end

res = calculatePositions(walkers)

# res = calculateProbDistri(walkers)

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


savefig(plt, "saves/task1h.density.pdf")

display(plt)


