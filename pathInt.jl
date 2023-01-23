V(x) = 1/2 * x^2

function S(x; Δt=1/4)
    res = 0
    N = length(x)
    for j in 2:N
        res += 1/2 * ((x[j] - x[j-1])/Δt)^2 + V((x[j]+x[j-1])/2)
    end
    return res
end


function step!(x, action, Δt)
    # Chose random pos in delta tau (not the ends)
    rng_idx = rand(2:N-1)
    xnew = copy(x)
    xnew[rng_idx] += (1-2*rand())*Δt

    actionNew = S(xnew; Δt)

    ΔS = abs(action-actionNew)
    
    # acceptance of random step 
    r = rand()
    p = exp(-Δt+ΔS)

    if r <= p
        action = actionNew
        x = xnew
    end
end

function integrate()
    t0 = 0
    tend = 100
    N = 400
    Δt = (tend-t0)/N

    x0 = 0
    xend = 1
    Δx = 0.01
    M = (xend-x0)/Δx

    lattice = rand(M, N)

    action = S(x; Δt)

    V_list = []

    for i in 1:10000
        push!(V_list, V.(x))
        step!(x,action,Δt)
    end

    return V_list
end
###################################

###################################

using Plots
using ProgressMeter
using Statistics

Vs = integrate()

plot(mean(Vs, dims=2))