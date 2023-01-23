include("pathInt.jl")

T(x, Δτ) =  1/2 * sum( ((x[2:end] - x[1:end-1])/Δτ).^2 )
V(x) = 1/2 * sum(((x[2:end] + x[1:end-1]) / 2).^2)

function integrate()
    t0 = 0
    tend = 100
    N = 400
    Δt = (tend-t0)/N
    steps = 20000
    p = Progress(steps)
    x = zeros(N)

    V_list = zeros(steps)
    T_list = zeros(steps)

    @showprogress for i in 1:steps
        V_list[i] =  V(x)
        T_list[i] = T(x, Δt)
        x = step(x,Δt,N)
        next!(p)
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


savefig(plt, "saves/task3a.energies.pdf")
display(plt)