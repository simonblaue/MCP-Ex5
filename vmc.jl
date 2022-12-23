using LinearAlgebra


const M = 300
const β = 1/2
const κ = 2

function localEnergy(r₁,r₂)
    r₁₂ = norm(r₁-r₂)
    u = 1 + α*r₁₂
    (κ-2)/norm(r₁) + (κ-2)/norm(r₂) + 1/r₁₂ * (1-2β/u^2) + 2*β*α/u^3 -κ^2 - β^2/u^4 + κ*β/u^2*(r₁/norm(r₁)- r₂/norm(r₂)) * (r₁-r₂)/r₁₂
end


function vmc()
    walkers = randn((M,2,2))
end

function initPos()
    
end


x = randn((3))
y = randn((3))

norm(x.*y)