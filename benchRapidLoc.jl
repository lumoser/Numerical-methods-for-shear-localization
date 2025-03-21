using LinearAlgebra
using SparseArrays
using GLMakie
using DifferentialEquations
using DelimitedFiles
using BenchmarkTools
using ProfileView
#=
this file contains the main solver with a parameter space that results in rapid localization
no plotting routines are implemented here as this is only for benchmarking purposes

=#
include("helperFunctions.jl")



GLMakie.activate!(title = "Benchmark", fxaa = false)

function mainSolver(Ny = 301, reltol = 1e-6, solver = "R5P" )


    secpyear    = 3600*24*365
    #----------------main switches-----------------  

    initTemp    = 550                           # Initial temperatureat BC
    Ny          = Ny                            # length of domain
    Δt          =  1e-4* secpyear               # Initial timestep
    dtmax       = 5e2*secpyear
    tmax        = 1e3*secpyear                  # end time of simulation
    ε_bg        = 1e-12                         # Background strainrate (Spang et al. 2024)4
    L           = -3000;                        # length of domain
    #vbc         = 1e-2 / secpyear              # absolute velocity at the boundary (only use if L = consant)
    vbc         = 2*ε_bg * L                    # Velocity as a function of L
    G           = 40e9                          # elastic Shear modulus
    σYield      = 2.4e9
    ηmin        = 1e15                          # regualrization viscosity
    η           = 1e23                          # background viscosity
    γ           = 0.1
    ρ           = 3300                          # density
    cp          = 1e3                           # heat capacity
    k           = 2.75                          # conductivity  

    
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
        P.ηs            .= (P.η .* exp.(-P.γ .*(P.TBuf .- P.Tref ))) .+ P.ηmin                  # updating temperature dependent viscosity
       
        nLV!(P.A, P.Ny, P.ηs, P.h)                                                              # updating coefficent matrix for velocities

        copyto!(P.vx, P.vBC)
        ldiv!(P.vx, lu!(P.A), P.vBC)                                                            # computing  velocities
        P.ϵs            .= 0.5 .* diff(P.vx) ./ P.h                                             # updating strainrate

        @views P.τBuf   .= U[P.Ny:end]   
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
        ct         = integrator.t / secpyear
        println("-----------current time: $ct years----------")
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
    timer = @benchmark $bb1() 
    #timer = 0
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
    
    
    println("method:", solver)  
    println("reltol:$reltol", ", abstol:$abstol")   

       
    println("adaptive nt: ", length(solution.t))   
    println("---------Solution done!-----------")

    return solution, timer

end

