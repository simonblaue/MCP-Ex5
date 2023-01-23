using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit
using LinearAlgebra
using StatsBase

include("pathInt.jl")

T(x, Δτ) =  1/2 * sum( ((x[2:end] - x[1:end-1])/Δτ).^2 )
V(x) = 1/2 * sum(((x[2:end] + x[1:end-1]) / 2).^2)


function integrate()
    t0 = 0
    tend = 100
    N = 400
    Δt = (tend-t0)/N
    steps = 20000

    x = zeros(N)
    
    h = zeros(steps)

    @showprogress for i in 1:steps
        x, changedx = step(x,Δt,N,h)
        h[i] = changedx
    end

    return h
end
###################################

###################################
res = integrate()

h = fit(Histogram, res, nbins=100)
h = normalize(h, mode=:pdf)

x = h.edges[1][1:end-1]
x_fit = LinRange(x[1], x[end], 100)

gaussian(x, μ, σ) = exp(-(x - μ)^2 / 2σ^2) / (√(2π) * σ)
@. model(x, p) = 0.5 * gaussian(x, p[1], p[2])

p0 = [0.0, 1.0]
fit_model = curve_fit(model, x, h.weights, p0)

plt = plot(h, label="")

(μ0, σ0) = coef(fit_model)
plot!(x_fit, model(x_fit, coef(fit_model)), 
    label="Gaussian Fit: μ₁=$(round(μ0, digits=2)), σ₁=$(round(σ0, digits=2))",
    lw=2)

savefig(plt, "saves/task3b.hist.pdf")
display(plt)
