using Plots
using LaTeXStrings
using Measurements
using ProgressMeter
using Distributed
using LsqFit
using LinearAlgebra
using StatsBase

include("pathInt.jl")
T(x, Δτ) =  1/2 * norm((x[2:end] - x[1:end-1])/Δτ)^2
V(x) = 1/2 * norm((x[2:end,:] + x[1:end-1,:])[:,1]/2 + (x[2:end,:]+ x[1:end-1,:])[:,2]*10 / 2)^2

wx = fill(1, 399)
wy = fill(10,399)

function integrate()
    t0 = 0
    tend = 100
    N = 400
    Δt = (tend-t0)/N
    steps = 200000
    p = Progress(2*steps)
    x = zeros(N,2)

    V_list = zeros(steps)
    T_list = zeros(steps)
    h = zeros(steps,2)

    for i in 1:steps
        V_list[i] =  V(x)
        T_list[i] = T(x, Δt)
        x, _ = step2d(x,Δt,N)
        next!(p)
    end

    for i in 1:steps
        x, h[i,:] = step2d(x,Δt,N)
        next!(p)
    end

    return V_list, T_list, h
end
###################################

###################################
using Plots
using ProgressMeter

Vs, Ts, h = integrate()

plt = plot(Vs, label="Potential")
plot!(Ts, label="Kinetic")

savefig(plt, "saves/task3c.energies.pdf")
display(plt)

lims = (min(h...), max(h...))

hist = histogram2d(h[:,1],h[:,2],
    bins=(100, 100), 
    normalize=:pdf, 
    xlims=lims,
    ylims=lims
    )
title!("Normalized 2D Histogram as squared wave function ")
xlabel!("x")
ylabel!("y")

savefig(hist, "saves/task3c.hist.pdf")
display(hist)

