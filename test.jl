t0 = 0
tend = 100
N = 400
Δt = (tend-t0)/N

x0 = 0
xend = 1
Δx = 0.01
M = floor(Int,(xend-x0)/Δx)

lattice = rand(M,N)

lattice[:,1] .= 0
lattice[:,end] .= 0

display(lattice)