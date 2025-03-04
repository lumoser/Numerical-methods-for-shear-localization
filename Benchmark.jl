using BenchmarkTools
using DelimitedFiles

include("benchRapidLoc.jl")


function spatialBmark()
    NyImplicit = [301, 351, 401, 451, 501, 601, 701, 801, 901, 1001]
    
    timesR4P   = zeros(length(NyImplicit))
    timesR5P   = zeros(length(NyImplicit))
    timesROS2  = zeros(length(NyImplicit))

    ntR4P      = zeros(length(NyImplicit))
    ntR5P      = zeros(length(NyImplicit))
    ntROS2     = zeros(length(NyImplicit))

    
    for i in 1:length(NyImplicit)

        sol1, timer1 = mainSolver(NyImplicit[i], 1e-6, "R4P")
        sol2, timer2 = mainSolver(NyImplicit[i], 1e-6, "R5P")
        sol3, timer3 = mainSolver(NyImplicit[i], 1e-6, "ROS2")

        ntR4P[i]    = length(sol1.t)
        ntR5P[i]    = length(sol2.t)
        ntROS2[i]   = length(sol3.t)

        timesR4P[i]  = timer1.times[1] *1e-9
        timesR5P[i]  = timer2.times[1] *1e-9
        timesROS2[i] = timer3.times[1] *1e-9

    end
    writedlm("bTimes/spatialRuntimeR4P.txt", timesR4P)
    writedlm("bTimes/spatialRuntimeR5P.txt", timesR5P)
    writedlm("bTimes/spatialRuntimeROS2.txt", timesROS2)

    writedlm("bTimes/spatialNTR4P.txt", ntR4P)
    writedlm("bTimes/spatialNTR5P.txt", ntR5P)
    writedlm("bTimes/spatialNTROS2.txt", ntROS2)
    

    NyExplicit  = [101, 151, 201, 227, 251, 277, 301, 327, 351]

    timesROCK2  = zeros(length(NyExplicit))
    ntROCK2     = zeros(length(NyExplicit))

    for i in 1:length(NyExplicit)

        sol4, timer4    = mainSolver(NyExplicit[i], 1e-6, "ROCK2")
        timesROCK2[i]   = timer4.times[1] *1e-9
        ntROCK2[i]      = length(sol4.t)

    end

    writedlm("bTimes/spatialRuntimeROCK2.txt", timesROCK2)
    writedlm("bTimes/spatialNTROCK2.txt", ntROCK2)

    return nothing

end

function tolTest()
    reltol = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
    reltolR = [5e-6, 1e-6, 1e-7, 1e-8 ,1e-9, 1e-10]

    rTimesLocR4P = zeros(length(reltol))
    rTimesLocR5P = zeros(length(reltol))

    dtLocR4P     = zeros(length(reltol)) 
    dtLocR5P     = zeros(length(reltol)) 


    for i in 1:length(reltol)
        sol1, timer1 = mainSolver(301, reltol[i],"R4P")
        sol2, timer2 = mainSolver(301, reltol[i],"R5P")

        rTimesLocR4P[i] = timer1.times[1]
        rTimesLocR5P[i] = timer2.times[1]

        dtLocR4P[i]     = round(minimum(compLocalMin(sol1.t)))
        dtLocR5P[i]     = round(minimum(compLocalMin(sol2.t)))


    end

    writedlm("bTimes/tolRuntimR4P.txt",rTimesLocR4P)
    writedlm("bTimes/tolRuntimR5P.txt",rTimesLocR5P)

    writedlm("bTimes/dtLocR4P.txt",dtLocR4P)
    writedlm("bTimes/dtLocR5P.txt",dtLocR5P)

    rTimesLocROCK2  = zeros(length(reltolR)) 
    dtLocROCK2  = zeros(length(reltolR)) 

    for i in 1:length(reltolR)
        sol1, timer1, = mainSolver(301, reltolR[i], "ROCK2")

        rTimesLocROCK2[i]   = timer1.times[1]
        dtLocROCK2[i]       = round(minimum(compLocalMin(sol1.t)))
    end

    writedlm("bTimes/rTimesLocROCK2.txt", rTimesLocROCK2)
    writedlm("bTimes/dtLocROCK2.txt", dtLocROCK2)

end

function highResSol()
    sol, timer = mainSolver()

    outT    = zeros(length(sol.u))
    outτ    = zeros(length(sol.u))
    for i in 1:length(sol.u)
        outT[i] = sol.u[i][HI]
        outτ[i] = sol.u[i][Ny+HI]
    end
     
    writedlm("FC/TFCRapid$Ny.txt", outT )
    writedlm("FC/tauFCRapid$Ny.txt", outτ )
    writedlm("FC/timesRapid$Ny.txt", solution.t)

end

