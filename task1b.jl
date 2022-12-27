include("vmc.jl")
using LaTeXStrings

function calcAveandStds()
    M = 300
    N = 10000
    n = 1000

    s = 0.1
    αs =[0., 0.1, 0.2, 0.3, 0.4, 0.5]


    avEnergies = []
    stdEnergies = []

    for α in αs
        walkers = initWalkers(M)
        avEnergy, stdEnergy = vmc(walkers, s, α, N, n)
        push!(avEnergies, avEnergy)
        push!(stdEnergies, stdEnergy)
    end

    return avEnergies, stdEnergies
end


@time avEnergies, stdEnergies = calcAveandStds()


labels = [L"\alpha = 0." L"\alpha = 0.1" L"\alpha = 0.2" L"\alpha = 0.3" L"\alpha = 0.4" L"\alpha = 0.5"]
p1 = plot(1:10, avEnergies, labels=labels, title="Average Energie")
p2 = plot(1:10, stdEnergies, labels=labels, title="Standard error Energie")

display(p1)
display(p2)

