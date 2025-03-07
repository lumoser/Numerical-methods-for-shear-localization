using LinearAlgebra
using GLMakie
using DelimitedFiles
using Statistics

include("../benchRapidLoc.jl")

solτ    = vec(readdlm("FC/tauFCRapid301.txt"))
solT    = vec(readdlm("FC/TFCRapid1001.txt"))
solTime = vec(readdlm("FC/timesRapid301.txt"))

Ny      = 301
relTol  = [5e-7, 1e-7, 5e-8, 1e-8, 5e-9,1e-9]
secpyear    = 3600*24*365
HI = Int32((Ny-1)/2)
maxErr1 = zeros(length(relTol))
maxErr2 = zeros(length(relTol))
maxErr3 = zeros(length(relTol))

err1 = zeros(length(solTime))
err2 = zeros(length(solTime))
err3 = zeros(length(solTime))
for j in 1:length(relTol)
    sol1 , timer1 = mainSolver(Ny, relTol[j], "R4P")
    sol2 , timer2 = mainSolver(Ny, relTol[j], "R5P")
    sol3 , timer3 = mainSolver(Ny, relTol[j], "ROCK2")


    
    for i in 1:length(solTime)
        
        err1[i] = abs(solτ[i] - sol1(solTime[i])[Ny+HI])
        err2[i] = abs(solτ[i] - sol2(solTime[i])[Ny+HI])
        err3[i] = abs(solτ[i] - sol3(solTime[i])[Ny+HI])

    end
    maxErr1[j] = maximum(err1)
    maxErr2[j] = maximum(err2)
    maxErr3[j] = maximum(err3)
end
fig1 = Figure(size = (1200,800),fontsize = 34)
ax1 = Axis(fig1[1,1], xlabel="reltol",ylabel = "Error τ [Pa]",xscale = log10, yscale=log10)
empty!(ax1)
lines!(relTol[end-2:end], maxErr1[end-2:end], linewidth = 3, label ="Rodas4P" )
lines!(relTol[end-2:end], maxErr2[end-2:end], linewidth = 3, label ="Rodas5P")
lines!(relTol[end-2:end], maxErr3[end-2:end], linewidth = 3, label ="ROCK2", color =:purple )

axislegend(ax1, merge = false, unique = false, position = :lt)

display(fig1)
save("./Figures/Figure_4_5.png", fig1)

