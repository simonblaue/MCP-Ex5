
T(x, Δτ) =  1/2 * sum( ((x[2:end] - x[1:end-1])/Δτ).^2 )
V(x) = 1/2 * sum(((x[2:end] + x[1:end-1]) / 2).^2)
S(x, Δτ) = T(x, Δτ) + V(x)


function step(x, Δt, N)
    # Chose random pos in delta tau (not the ends)
    rng_idx = rand(2:N-1)
    xnew = copy(x)
    xnew[rng_idx] += (1-2*rand())*Δt

    # actionNew = S(xnew, Δt)

    ΔS = -(S(xnew, Δt) - S(x, Δt))
    
    # acceptance of random step 
    r = rand()
    p = exp(-Δt+ΔS)

    if r < p
        x = xnew
    end

    return x
end

function integrate()
    t0 = 0
    tend = 100
    N = 400
    Δt = (tend-t0)/N
    steps = 20000

    x = zeros(N)
    
    V_list = zeros(steps)
    T_list = zeros(steps)

    @showprogress for i in 1:steps
        V_list[i] =  V(x)
        T_list[i] = T(x, Δt)
        x = step(x,Δt,N)
    end

    return V_list, T_list, x
end
###################################

###################################
using Plots
using ProgressMeter

Vs, Ts, x = integrate()

plt = plot(Vs, label="Potential")
plot!(Ts, label="Kinetic")

hist = histogram(x, bins=100)

display(hist)
display(plt)