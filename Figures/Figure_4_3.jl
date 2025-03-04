using LinearAlgebra
using GLMakie
using DelimitedFiles


reltol = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

rTimesLocR4P = vec(readdlm("bTimes/tolRuntimR4P.txt")) .*1e-9
rTimesLocR5P = vec(readdlm("bTimes/tolRuntimR5P.txt")) .*1e-9

dtLocR4P = vec(readdlm("bTimes/dtLocR4P.txt"))
dtLocR5P = vec(readdlm("bTimes/dtLocR5P.txt"))

reltolR = [5e-6, 1e-6, 1e-7, 1e-8 ,1e-9, 1e-10]

rTimesLocROCK2 = vec(readdlm("bTimes/rTimesLocROCK2.txt")) .*1e-9
dtLocROCK2     = vec(readdlm("bTimes/dtLocROCK2.txt"))

lims = (1e-4, 1e-10, nothing, nothing)

fig1 = Figure(size = (1200, 1200), fontsize =32)
ax1 = Axis(fig1[1,1], xscale = log10, yscale = log10,xlabel ="reltol", ylabel="Timestep at localization [s]", title = "Localization Δt implicit ", xreversed = true)
ax2 = Axis(fig1[1,2], xscale = log10,xlabel ="reltol", ylabel="Timestep at localization [s]", title = "Localization Δt explicit ", xreversed  =true)
lines!(ax1, reltol, dtLocR4P , label ="Rodas4P", linewidth=3)
lines!(ax1, reltol, dtLocR5P , label ="Rodas5P", linewidth=3)
lines!(ax2, reltolR, dtLocROCK2 , label ="ROCK2", color=:purple, linewidth=3)

axislegend(ax1, merge = false, unique = false, position = :rt)
axislegend(ax2, merge = false, unique = false, position = :rt)

ax3 = Axis(fig1[2,1], xscale = log10, xlabel ="reltol", ylabel="Runtime [s]", title = "Runtime implicit ", xreversed = true)
ax4 = Axis(fig1[2,2], xscale = log10, xlabel ="reltol", ylabel="Runtime [s]", title = "Runtime explicit ", xreversed  = true)

lines!(ax3, reltol, rTimesLocR4P, label = "Rodas4P", linewidth = 3)
lines!(ax3, reltol, rTimesLocR5P, label = "Rodas5P", linewidth = 3)

lines!(ax4, reltolR, rTimesLocROCK2, label = "ROCK2", linewidth = 3, color =:purple)

display(fig1)
#save("Figures/Figure_4_3.png", fig1)
