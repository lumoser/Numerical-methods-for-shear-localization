using LinearAlgebra
using GLMakie
using DelimitedFiles

function compSlope(x, y)
    k = zeros(length(y)-1)
    for i in 2:length(y)
        k[i-1] = ((y[i] - y[i-1])/(x[i]-x[i-1])) 
    end
    return round(mean(k), digits = 1)
end


rTimesLocR4P = vec(readdlm("bTimes/tolRuntimR4P.txt")) .*1e-9
rTimesLocR5P = vec(readdlm("bTimes/tolRuntimR5P.txt")) .*1e-9

rTimesLocROCK2 = vec(readdlm("bTimes/rTimesLocROCK2.txt")) .*1e-9


ntTolR4P = vec(readdlm("bTimes/ntTolR4P.txt"))
ntTolR5P = vec(readdlm("bTimes/ntTolR5P.txt"))

ntTolROCK2 = vec(readdlm("bTimes/ntTolROCK2.txt"))

kR4P = compSlope(log10.(ntTolR4P), log10.(rTimesLocR4P))
kR5P = compSlope(log10.(ntTolR5P), log10.(rTimesLocR5P))
kROCK2 = compSlope(log10.(ntTolROCK2), log10.(rTimesLocROCK2))


fig1 = Figure(size = (800, 800), fontsize =32)
ax1 = Axis(fig1[1,1], xscale = log10, yscale = log10, xlabel ="Number of timesteps", ylabel="Runtime [s]", title = "Runtimes")

lines!(ax1, ntTolR4P, rTimesLocR4P , label ="Rodas4P, k ≈ $kR4P", linewidth=3)
lines!(ax1, ntTolR5P, rTimesLocR5P , label ="Rodas5P k ≈ $kR5P ", linewidth=3)
lines!(ax1, ntTolROCK2, rTimesLocROCK2 , label ="ROCK2, k ≈ $kROCK2", color=:purple, linewidth=3)

axislegend(ax1, merge = false, unique = false, position = :lt)

display(fig1)
#save("Figures/Figure_4_2.png", fig1)