using LinearAlgebra
using GLMakie
using DelimitedFiles
using Statistics


function compSlope(x, y)
    k = zeros(length(y)-1)
    for i in 2:length(y)
        k[i-1] = ((y[i] - y[i-1])/(x[i]-x[i-1])) 
    end
    return round(mean(k), digits = 1)
end


NyImplicit = [301, 351, 401, 451, 501, 601, 701, 801, 901, 1001]
NyExplicit  = [101, 151, 201, 227, 251, 277, 301, 327, 351]

timesR4P   = vec(readdlm("bTimes/spatialRuntimeR4P.txt"))
timesR5P   = vec(readdlm("bTimes/spatialRuntimeR5P.txt"))
timesROS2  = vec(readdlm("bTimes/spatialRuntimeROS2.txt"))
timesROCK2 = vec(readdlm("bTimes/spatialRuntimeROCK2.txt"))

ntR4P      = vec(readdlm("bTimes/spatialNTR4P.txt"))
ntR5P      = vec(readdlm("bTimes/spatialNTR5P.txt"))
ntROS2     = vec(readdlm("bTimes/spatialNTROS2.txt"))
ntROCK2    = vec(readdlm("bTimes/spatialNTROCK2.txt"))

kR4P       = compSlope(log10.(NyImplicit), log10.(timesR4P))
kR5P       = compSlope(log10.(NyImplicit), log10.(timesR5P))
kROS2      = compSlope(log10.(NyImplicit), log10.(timesROS2))
kROCK2     = compSlope(log10.(NyExplicit), log10.(timesROCK2))

fig1 = Figure(size = (1200, 600), fontsize =32)
ax1 = Axis(fig1[1, 1], xscale = log10, yscale = log10, xlabel= "Number of gridpoints", ylabel="Runtime [s]", title = "Runtimes" )
ax2 = Axis(fig1[1,2],xscale = log10, yscale = log10, ylabel = "Number of timesteps",  xlabel= "Number of gridpoints", yaxisposition = :right, title = "Number of timesteps")
lines!(ax1, NyImplicit ,timesR4P, linewidth = 5, label = "k ≈ $kR4P")
lines!(ax1, NyImplicit , timesR5P, linewidth = 5, label = "k ≈ $kR5P")
lines!(ax1, NyImplicit , timesROS2, linewidth = 5, label = "k ≈ $kROS2")
lines!(ax1, NyExplicit , timesROCK2, linewidth = 5, label = "k ≈ $kROCK2", color = :purple)

lines!(ax2, NyImplicit, ntR4P, label = "Rodas4P", linewidth = 5, )
lines!(ax2, NyImplicit, ntR5P, label = "Rodas5P", linewidth = 5)
lines!(ax2, NyImplicit, ntROS2, label = "ROS2", linewidth = 5)
lines!(ax2, NyExplicit, ntROCK2, label = "ROCK2", linewidth = 5, color = :purple)

axislegend(ax1, merge = false, unique = false, position = :lt)
axislegend(ax2, merge = true, unique = true, position=:lb)

display(fig1)
#save("Figures/Figure_4_1.png", fig1)


