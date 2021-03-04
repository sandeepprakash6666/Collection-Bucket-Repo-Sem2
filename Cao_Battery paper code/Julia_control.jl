
using DifferentialEquations
using Plots 
    plotlyJS()




mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end



function f(du,u,p,t)
    du[1] = -0.5*u[1] + u.f1
    du[2] = -0.5*u[2]
end

const tstop1 = [5.]
const tstop2 = [8.]


function condition(u,t,integrator)
    t in tstop1
  end
  
  function condition2(u,t,integrator)
    t in tstop2
  end




  function affect!(integrator)
    for c in full_cache(integrator)
      c.f1 = 1.5
    end
  end

  function affect2!(integrator)
    for c in full_cache(integrator)
      c.f1 = -1.5
    end
  end



  save_positions = (true,true)

  cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
  
  save_positions = (false,true)
  
  cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
  
  cbs = CallbackSet(cb,cb2)

  u0 = SimType([10.0;10.0], 0.0)
  prob = ODEProblem(f,u0,(0.0,10.0))

  const tstop = [5.;8.]
  sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)

  sol[1].f1
  
  plot(sol)




  ##*Using Parameterized Functions
  using DifferentialEquations
  using Plots 
    plotlyjs()

  #Define ODE eqns
    function f(du,u,p,t)
        du[1] = -0.5*u[1] + p
        du[2] = -0.5*u[2]
    end
    
    #
    mutable struct SimType{T} <: DEDataVector{T}    #?Ask evren here
        x::Array{T,1}
        f1::T
    end

    #nitial Conditions
    u0 = SimType([10.0;10.0], 0.0)
    p = 0.0
    tspan = (0.0, 10.0)
    prob = ODEProblem(f,u0,tspan,p)
    
    const tstop = [5.;8.]
    const tstop1 = [5.]
    const tstop2 = [8.]

        save_positions = (true,true)
        function condition(u,t,integrator)
            t in tstop1
        end
      
        function affect!(integrator)
            integrator.p = 1.5
        end

    cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
    
        save_positions = (false,true)
        function condition2(u,t,integrator)
            t in tstop2
        end
        function affect2!(integrator)
            integrator.p = -1.5
          end
    cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
    
    cbs = CallbackSet(cb,cb2)




    sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)



    plot(sol)


##*Parameterized Function : Bioreactor


    using DifferentialEquations
    using Plots 
      plotlyjs()
  
    #Define ODE eqns
    function f(du,u,p,t)
            #Parameters
            x2_f   = 4
            km     = 0.12
            k1     = 0.4545
            Y      = 0.4
            μ_max  = 0.53


        x1, x2 = u
        D = p
        μ = μ_max*x2/(km + x2 + k1*x2^2)

        du[1] = x1*(μ - D)
        du[2] = D*(x2_f - x2) - μ*x1/Y  

    end

    mutable struct SimType{T} <: DEDataVector{T}    #?Ask evren details
        x::Array{T,1}
        f1::T
    end

    #nitial Conditions [1, 1] @ time 0
    u0 = SimType([1.0;1.0], 0.0)
    p = 0.2916572110275656
    tspan = (0.0, 20.0)
    prob = ODEProblem(f,u0,tspan,p)


    #Callback Parameters
    const tstop = [5.;10.;15.]
    const tstop1 = [5.]
    const tstop2 = [10.]
    const tstop3 = [15.]

    #1st callback break
    save_positions = (true,true,true)               #?Ask evren details
    function condition(u,t,integrator)
        t in tstop1
    end
  
    function affect!(integrator)
        integrator.p = 0.3416572201748369
    end
cb = DiscreteCallback(condition, affect!, save_positions=save_positions)


    save_positions = (false,true,true)
    function condition2(u,t,integrator)
        t in tstop2
    end
    function affect2!(integrator)
        integrator.p = 0.3916572288927861
    end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)


    save_positions = (false,false, true)
    function condition3(u,t,integrator)
        t in tstop3
    end
    function affect3!(integrator)
        integrator.p = 0.33174142317075545
    end
cb3 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)

cbs = CallbackSet(cb,cb2,cb3)



sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)



plot(sol,ylim = [0.0, 2.0], title = "Sim as ODE : Callback_ParamFun "  )



