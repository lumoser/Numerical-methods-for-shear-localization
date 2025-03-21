using LinearAlgebra
using SparseArrays
using GLMakie
using DifferentialEquations
using DelimitedFiles
using BenchmarkTools
using ProfileView


include("../helperFunctions.jl")



GLMakie.activate!(title = "Benchmark", fxaa = false)

function mainSolver(Ny = 301, reltol = 1e-6, solver = "R5P" )


    secpyear    = 3600*24*365
    #----------------main switches-----------------  

    initTemp    = 550                           # Initial temperatureat BC
    Ny          = Ny                            # length of domain
    Δt          =  1e-4* secpyear               # Initial timestep
    dtmax       = 5e2*secpyear
    tmax        = 1e3*secpyear                  # end time of simulation
    ε_bg        = 1.3e-12                         # Background strainrate (Spang et al. 2024)4
    L           = -3000;                        # length of domain
    #vbc         = 1e-2 / secpyear              # absolute velocity at the boundary (only use if L = consant)
    vbc         = 2*ε_bg * L                    # Velocity as a function of L
    G           = 90e9                          # elastic Shear modulus
    σYield      = 2.4e9
    ηmin        = 1e19                          # regualrization viscosity
    η           = 1e23                          # background viscosity
    γ           = 0.1
    ρ           = 3300                          # density
    cp          = 1e3                           # heat capacity
    k           = 4.75                          # conductivity  

    
    #-------solver options-------------
    solvers = Dict(
    "R4P" => Rodas4P(autodiff = false),
    "R5P" => Rodas5P(autodiff = false),
    "ROS2" => ROS2(autodiff = false),
    "ROCK2" => ROCK2(),
    "ROCK4" => ROCK4()
    )

    chSolver    = solvers[solver] 

    reltol      = reltol        # default = 1e-3
    abstol      = 1e-7          # default = 1e-6
    gamma       = 0.9           # default = 0.9
    qmin        = 0.2           # default = 0.2
    qmax        = 13            # default = 10


    #--------init Condition---------------

    yu          = LinRange(0,1,Ny-1)                                    # unified grid (for gaussian pulse)
    Temp        =    initTemp .* ones(Ny-1)                             # temperature init
    Temp      .+= 50 .* exp.(-0.5 .* ((yu    .- 0.5) ./ 0.1) .^4)       # Gaussian pulse for temperature
    Tref        = ones(Ny-1) .*initTemp                                 # reference temperature for viscosity

    HI          = Int16((Ny-1)/2)                                       # Index for domain half for cell centers
    
    
    #-------setup halfspace cooling----------
    #Temp[1:HI]  = 1200 * erf.(abs.(yc[1:HI]) ./(2 *  sqrt(k1 * 2e7*secpyear))) 
    #Temp[HI:end]  = 1200 * erf.(abs.(yc[HI:end]) ./(2 *  sqrt(k2 * 2.5e8 *secpyear))) 

    #-----------------DiffEq Init------------------


    xc, D, h   = D2nd(Ny-1,L)                      # Matrix for 2nd derivative of vector (FD) and the coordinates to cell centers
    x          = range(0,L,Ny)                     # coordinates of the cell edges

  

    uInit   = zeros(2*(Ny -1))                      # Vector holding the initial values of τ and T
    uInit[1:Ny-1] = Temp                            # setting initial temperature, inital τ = 0

    vBC         = zeros(Ny)                         # velocity BC init
    vBC[1]      = vbc                               # velocity BC top
    vBC[end]    = 0 * vbc                           # velocity BC bottom

    vcm         = vbc *secpyear*1e2                 # velocity in cm/ year (for plotting)


    #----Init for other parametes-----
    ηs          = ones(Ny-1) .* η                                   # nonlinear viscosity
    A           = Tridiagonal(zeros(Ny-1), zeros(Ny), zeros(Ny-1))  # ceofficent matrix for velocity LSE
    ϵs          = zeros(Ny-1)                                       # strainrate
    vx          = zeros(Ny)                                         # velocity
    τII         = zeros(1)                                          # second invatriant of stress tensor
    pFlag       = ones(1)                                           # ceofficent for plasticity (can either be 1 or 0)
    TBuf        = copy(Temp)                                        # buffer for temperature
    τBuf        = similar(TBuf)                                     # buffer for stress

    P = (;Ny, ηs, vx, ϵs, A, σYield, ηmin, γ, G, D, h,  ρ, cp, k, vBC, η, Tref, τII, pFlag, HI, TBuf, τBuf)

    println("---------Started solver with Ny: $Ny -----------")

    function shearElastoPlastic!(dU,U,P,t)

        @views P.TBuf     .= U[1:P.Ny-1]
        P.ηs            .= (P.η .* exp.(-P.γ .*(TBuf .- P.Tref ))) .+ P.ηmin             # updating temperature dependent viscosity
       
        nLV!(P.A, P.Ny, P.ηs, P.h)                                                              # updating coefficent matrix for velocities
                               
        copyto!(P.vx, P.vBC)
        ldiv!(lu!(P.A),P.vx )                        # computing  velocities

        
        P.ϵs            .= 0.5 .* diff(P.vx) ./ P.h 
        @views P.τBuf   .= U[P.Ny:end]                                     # updating strainrate
        P.τII           .= sqrt(0.5 *(P.τBuf[P.HI]^2))                                        # updating second invariant of stress tensor
            
        dU[P.Ny:end]    .=  P.pFlag[1] .* ((2* P.G .* P.ϵs) .- (P.τBuf  .* P.G ./ P.ηs))   # solving dτ/dt = 2Gϵ - G* (τ/η) (elastic case, if plastic = 0)


        # solving ∂T/∂t = (κ ∂T²/∂y² + 2τϵ)/ρcₚ with ϵ = ϵ_total - ϵ_elastic, ϵ_elastic = 1/2G * dτ/dt  
             
        dU[1:Ny-1]      .= ((P.k .* *(P.D, P.TBuf)).+ (2 .*(dU[P.Ny:end] .+ P.τBuf).* (P.ϵs .-(1/2G .* dU[P.Ny:end])))) ./(P.ρ * P.cp) 

        dU[1]           = 0     # Dirichlet BC for temperature
        dU[Ny-1]        = 0     

    end


    function saveVals!(u, t, integrator)

        ηi         = P.ηs         
        vi         = P.vx  
        ϵi         = P.ϵs

        return ηi, ϵi, vi
    end


    function plasticityTrigger(u,t,integrator)

        integrator.p.σYield- sqrt(0.5* u[integrator.p.Ny + integrator.p.HI]^2) 

    end

    function addPlasticity!(integrator)

        println("Plasticity triggered @ ",round(integrator.t/secpyear), " years")
        integrator.p.pFlag[1] = 0

    end





    diffTerm = ODEProblem(shearElastoPlastic!, uInit, (0,tmax),P)                               # Init the ODEProblem
    
    plasticCallback = ContinuousCallback(plasticityTrigger, addPlasticity!)                     # Callback for plasticity   
    
    savedValues = SavedValues(Float64, Tuple{Vector{Float64},Vector{Float64},Vector{Float64}})  # callback for saving parameterEvolution  
    callback = SavingCallback((u,t,integrator) -> (saveVals!(u,t,integrator)), savedValues)

    cbs1     = CallbackSet(callback,plasticCallback)                                            # combining the callbacks to a set
    
    bb1 = () -> solve(
                    diffTerm, 
                    chSolver,                                      
                    reltol      = reltol,
                    abstol      = abstol,
                    dt          = Δt,
                    adaptive    = true, 
                    dtmax       = dtmax, 
                    callback    = cbs1,
                    gamma       = gamma,
                    qmin        = qmin,
                    qmax        = qmax
                    )
    #timer = @benchmark $bb1() 
    timer = 0
    
    println("method:", solver)  
    println("reltol:$reltol", ", abstol:$abstol")   

    
    
    solution = solve(
        diffTerm, 
        chSolver,                   
        reltol      = reltol,
        abstol      = abstol,
        dt          =Δt,
        adaptive    = true, 
        dtmax       = dtmax, 
        callback  = cbs1,
        gamma       = gamma,
        qmin        = qmin,
        qmax        = qmax
    )
    
      
    #----------plotting section----------
    evoPlt = false
    τplot  = true
    p2f    = false
    
    if evoPlt
        function evoPltAdaptive!(fig::Figure, ax1::Axis, ax2::Axis, ax3::Axis, ax4::Axis)
            empty!(ax1)
            empty!(ax2)
            empty!(ax3)
            empty!(ax4)
            lines!(ax1, Temp, xc/1e3, color=:blue,  linewidth = 3,label = "Init T")   # Initial condition
        
            #--------------DifferentialEquations--------------
            lines!(ax1,solution.u[i][1:Ny-1], xc/1e3, color =:red, linestyle=:solid, linewidth = 3, label = "T")
            lines!(ax4, log10.(savedValues.saveval[j][1])  , xc/1e3, color = :red, linewidth = 3, linestyle=:solid, label = "η")
            lines!(ax3, log10.(abs.(savedValues.saveval[j][2])), xc/1e3, color = :red, linewidth = 3,linestyle=:solid, label = "ϵ");
            lines!(ax2, savedValues.saveval[j][3] .*secpyear*1e2, x/1e3, color = :red, linewidth = 3, linestyle=:solid, label = "Vₓ");
            tyrs        = solution.t[i] /secpyear
            ax1.title   = rpad("Time=$(round(tyrs, digits=5)) yrs",20)
            
            sleep(0.0005)
            display(fig)
        end

        fig = Figure(size = (1450, 850), fontsize =34)
        tempLim = maximum(solution.u[end][HI] +50 )
        limits = ((500, tempLim,nothing, nothing))
        ax1 = Axis(fig[1, 1], xlabel="Temperature [°C]", ylabel="Depth y [km]", title="Time=$time myr",limits=limits)
        ax2 = Axis(fig[1, 2], xlabel="Velocity [cm/year]")  
        ax3 = Axis(fig[1, 3], xlabel="Log shear strain rate [1/s]")#, limits = ((-35,0, nothing,nothing)))  
        ax4 = Axis(fig[1, 4], xlabel="Log viscosity [Pa s]")
        i = 2           # 2 if plastic bc callbacks mess 
        j = 1
        println(length(solution.u))
        println(length(savedValues.saveval))
        for p in 1:length(solution.u)-1 
            evoPltAdaptive!(fig, ax1, ax2, ax3, ax4)
            i += 1
            j += 1           
        end
    
        axislegend(ax1, merge = true, unique = true)
        axislegend(ax2, merge = true, unique = true)
        axislegend(ax3, merge = true, unique = true)
        axislegend(ax4, merge = true, unique = true)
        if p2f
            save("Figures/ParamEvoPlastic1_0.png", fig)
        end

    elseif τplot 
        function stressPlt!(TmaxDE, timePlot, p2f::Bool)
            fig2 = Figure(size = (1600,800),fontsize = 34)
            tlim = tmax*1e-3 /secpyear 
            lims = ((0,tlim,0,3900))
            #tempLim = maximum(TmaxDE ) +40
            tempLim = 950
            ax7 = Axis(fig2[1,1], xlabel="Time [kyr]",ylabel = "Shear stress [MPa]", limits = lims)
            ax8 = Axis(fig2[1,1],ylabel = "Temperature [°C]",xlabel = "",yaxisposition = :right, limits = (0,tlim,550,tempLim))
            ax9 = Axis(fig2[1,2], ylabel = "Stepsize [yr]", xlabel = "Time[kyr]")  
        
            lines!(ax8, timePlot ./secpyear .*1e-3, TmaxDE, color = :red, linewidth = 3, label = "T")
            lines!(ax7, timePlot ./secpyear .*1e-3, tauplot .*1e-6, color = :blue, linewidth = 3 , label = "τxy")
        
            sSize = diff(solution.t) ./secpyear
            lines!(ax9, solution.t[1:end-1] ./secpyear .*1e-3, sSize,color = :black, label = solver , linewidth = 3)
            vpl     = Int64(abs(round(vcm)))
            dtmin   = Int64(round(minimum(diff(solution.t))))
            ndt     = length(solution.t)
            #yield   = Int64(round(σYield *1e-9))   
            yield   =  round(σYield *1e-9, digits =3)
            gT      = Int64(round(G *1e-9))
            bLine = 200
            posT     = 0.4
            
            
            dts = compLocalMin(solution.t)
            dtlocy = round( minimum(dts[1:end-1]) / secpyear  ,digits = 4)
            dtloc = round( minimum(dts[1:end-1])   ,digits = 1) 
            println("dtloc(sec)", minimum(dts[1:end-1]))
            println("dtloc(yrs)", dtlocy)
            #text!(ax7, "dtLoc = $dtlocy years", position= Point2f(posT, bLine +1200))
            #text!(ax7, "dtloc = $dtloc seconds", position= Point2f(posT, bLine +1200))
        
            #text!(ax7, "τyield = $yield GPa", position= Point2f(posT, bLine +1200))
            
          

            text!(ax7, "G = $gT GPa", position= Point2f(posT,bLine +800))
            text!(ax7, "vₓ = $vpl cm/year", position= Point2f(posT, bLine+1000))

            text!(ax7, "nSteps = $ndt ", position= Point2f(posT,bLine +1400))
            #text!(ax9, "reltol = $reltol" , position = Point2f(0,24) )
            #text!(ax9, "abstol = $abstol" , position = Point2f(0,21))
            axislegend(ax7, merge = false, unique = true, position =:lt)
            axislegend(ax8, merge = false, unique = true, position =:rt)
            axislegend(ax9, merge = false, unique = true, position =:rt)
            display(fig2)
            t_years = Int64(round(timePlot[end] /secpyear))
            if p2f
                save("Figures/Figure_3_4.png", fig2)
            end
        end

        plotsteps = 10000
        timePlot            = LinRange(0, tmax, plotsteps)
        tauplot          = zeros(plotsteps)
        TmaxDE           = zeros(plotsteps)
        for i in 1:length(timePlot)
            tauplot[i] = solution(timePlot[i])[Ny+HI]
            TmaxDE[i]  = solution(timePlot[i])[HI]  
        end
        stressPlt!(TmaxDE, timePlot, p2f)
    

     
    end
       
    println("adaptive nt: ", length(solution.t))   
    println("---------Solution done!-----------")

    return solution, timer

end




sol, timer =  mainSolver(501, 1e-4, "R4P")
display(timer)