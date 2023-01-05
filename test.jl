using Plots
using ProgressBars

# Threads.nthreads()

M = 1000000
cumbersome = zeros(M)

pusher = []

Threads.@threads for i in ProgressBar(1:M)
    res = 2*i+3/i
    cumbersome[i] = res
    # push!(pusher, res)
end
