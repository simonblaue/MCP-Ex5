
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

    return x, xnew[rng_idx]
end


function step2d(x, Δt, N)
    # Chose random pos in delta tau (not the ends)
    rng_idx = rand(2:N-1)
    xnew = copy(x)

    xnew[rng_idx,:] .+= (1 .- 2*rand(2))*Δt

    ΔS = -(S(xnew, Δt) - S(x, Δt))
    
    # acceptance of random step 
    r = rand()
    p = exp(-Δt+ΔS)
    if r < p
        x = xnew
    end

    return x, xnew[rng_idx,:]
end