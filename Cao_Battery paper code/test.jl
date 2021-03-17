using Sundials



##*Single Parameter - VanderPol

  p = [1.0]
  function f!(du,u,p,t)

    du[1] = (1.2 - u[2]^2)*u[1] - u[2] + p[1]
    du[2] = u[1]    

    end

  #Initial Conditions
  u₀ = [0.01, 1.0]
  tspan = (0.0,10.0)

  prob = ODEProblem(f!,u₀,tspan,p)

  sol = solve(prob)
  plot(sol)



##*2 Inputs - Vanderpol

  using DifferentialEquations
  using Plots
    plotlyjs()
  using Interpolations


  #*To pass the correct MV to integrator
  function my_hold(u, t)
    
    n = round(t)
    Int_n = convert(Int, n) 

    return(u[Int_n])
  end

  # A = [1.0, 10.0]
  # itp = interpolate((A,), A, Gridded(Linear()))
  # itp = interpolate((A,), A, Gridded(Constant()))
  # itp
  # itp(4.5)
  # itp(8.5)

  # p2 = [1.0, 0.5]

  #*Defining model ODEs
  function f2!(du,u,p,t)

  # param = itp(t)
  param = my_hold(p, t)

  du[1] = (1.2 - u[2]^2)*u[1] - u[2] + param
  du[2] = u[1]  

  end

  #Initial Conditions
  u₀ = [0.01, 1.0]
  tspan = (1.0,10.0)
  p2 =     u = [1.0; 5*ones(5); 4*ones(4) ]


  prob = ODEProblem(f2!,u₀,tspan,p2)
  sol = solve(prob)
  plot(sol)


##*Input Profile: bioreactor (As ODE) 
  #!Implementation questionable at best 
  #Bcz of adhoc floor and if statements : confuses integrator timestepping algos
    
  using DifferentialEquations
  using Plots
    plotlyjs()

  tspan = (0.0, 20.0)
  function step_input(t)    #function to return correct step input
            # t = 0.0

    u_step = [0.049172880221869476
            0.0991728751856224
            0.14917284616952867
            0.19917261848265047
            0.21487352821766373
            0.22555820123690046
            0.27513080807784135
            0.3251307693218697
            0.3751307683754231
            0.4251307730448192
            0.4751307777520993
            0.4696865127630823
            0.31978417121221314
            0.2969992371566889
            0.2992734614688022
            0.30039102996976513
            0.3003726578943411
            0.3002515488631327
            0.3001784222320125
            0.30013583848070924
          ]    
      NFE = size(u_step)[1]

    nfe = floor(Int, t/(tspan[2] - tspan[1])*NFE +1)

    if nfe > NFE
      nfe = NFE
    elseif nfe < 1
      nfe = 1 
    end

    u_step[nfe]

    return u_step[nfe]

  end
  # step_input(0.256)

  function f_bioreactor!(du,u,p,t)
          # du = [0.0, 0.0]
          # u = [1.0, 1.0]
          # p = step_input
          # t = 0.26

          #Parameters
          x2_f   = 4
          km     = 0.12
          k1     = 0.4545
          Y      = 0.4
          μ_max  = 0.53


    x1, x2 = u
    D = p(t)
    μ = μ_max*x2/(km + x2 + k1*x2^2)

    du[1] = x1*(μ - D)
    du[2] = D*(x2_f - x2) - μ*x1/Y  

  end

  #Initial Conditions
  u₀ = [1.0, 1.0]
  p = step_input


  alg = Vern7() 
  prob = ODEProblem(f_bioreactor!,u₀,tspan,p)
  sol = solve(prob,  alg )
  plot(sol, ylim = [0.0, 2.0], title = "Sim as ODE : If Statements ? : $alg "  )




##*Input Profile: bioreactor (As DAE) 
  #!Implementation is Questionable at best

  using DifferentialEquations
  using Sundials
  using Plots
    plotlyjs()

  tspan = (0.0, 20.0)
  function step_input(t)
            # t = 0.0

    u_step = [0.049172880221869476
            0.0991728751856224
            0.14917284616952867
            0.19917261848265047
            0.21487352821766373
            0.22555820123690046
            0.27513080807784135
            0.3251307693218697
            0.3751307683754231
            0.4251307730448192
            0.4751307777520993
            0.4696865127630823
            0.31978417121221314
            0.2969992371566889
            0.2992734614688022
            0.30039102996976513
            0.3003726578943411
            0.3002515488631327
            0.3001784222320125
            0.30013583848070924
          ]

    NFE = size(u_step)[1]

    nfe = floor(Int, t/(tspan[2] - tspan[1])*NFE +1)

    if nfe > NFE
      nfe = NFE
    elseif nfe < 1
      nfe = 1 
    end

    u_step[nfe]

    return u_step[nfe]

  end
  step_input(0.256)


  function f_bioreactor!(out, du,u,p,t)
  #         out = [0.0, 0.0, 0.0]
  #         du = [0.0, 0.0, 0.0]
  #         u = [1.0, 1.0, 0.0]
  #         p = step_input
  #         t = 0.26

    #Parameters
    x2_f   = 4
    km     = 0.12
    k1     = 0.4545
    Y      = 0.4
    μ_max  = 0.53


    x1, x2, μ = u
    D = p(t)


    out[1] = x1*(μ - D)                     - du[1]
    out[2] = D*(x2_f - x2) - μ*x1/Y         - du[2]
    out[3] = μ_max*x2/(km + x2 + k1*x2^2)   - μ
  end

  #Initial State (alg state must be feasible)
  u₀ = [1.0, 1.0, 0.33661479834868213]
  du₀ = [0.0, 0.0, 0.0]
  tspan = (0.0,20.0)
  p = step_input
  differential_vars = [true,true,false]

  prob = DAEProblem(f_bioreactor!,du₀,u₀,tspan,p,differential_vars=differential_vars)

  sol = solve(prob,IDA())

  plot(sol, ylim = [0.0, 2.0], title = "Sim as DAE : If statements ?")




##*Naive approach : Calling integrator in for loop

  using DifferentialEquations
  using Plots
    plotlyjs()


  function f_bioreactor!(du,u,p,t)
    # du = [0.0, 0.0]
    # u = [1.0, 1.0]
    # p = step_input
    # t = 0.26

    #Parameters
    x2_f   = 4
    km     = 0.12
    k1     = 0.4545
    Y      = 0.4
    μ_max  = 0.53


    # x1 = u[1]
    # x2 = u[2]
    # D = p
    # μ = μ_max*u[2]/(km + u[2] + k1*u[2]^2)

    du[1] = u[1]*(      (μ_max*u[2]/(km + u[2] + k1*u[2]^2))   - p)
    du[2] = p*(x2_f - u[2]) - (   μ_max*u[2]/(km + u[2] + k1*u[2]^2)      )*u[1]/Y  


    # μ = μ_max*x2/(km + x2 + k1*x2^2)

    # du[1] = x1*(μ - D)
    # du[2] = D*(x2_f - x2) - μ*x1/Y 
  end

  u_step = [0.049172880221869476
            0.0991728751856224
            0.14917284616952867
            0.19917261848265047
            0.21487352821766373
            0.22555820123690046
            0.27513080807784135
            0.3251307693218697
            0.3751307683754231
            0.4251307730448192
            0.4751307777520993
            0.4696865127630823
            0.31978417121221314
            0.2969992371566889
            0.2992734614688022
            0.30039102996976513
            0.3003726578943411
            0.3002515488631327
            0.3001784222320125
            0.30013583848070924
          ]
    
  # u_step = [0.2916572110275656 
  #             0.3416572201748369 
  #             0.3916572288927861 
  #             0.33174142317075545
  #           ]

    NFE = size(u_step)[1]
    Tspan = 20 

    # x_end = zeros(1,2)
    x_end = [1.0, 1.0]
    
    x_list = []
    t_list = []
    
    for i in 1:NFE
      # i = 2

      global x_end
      global x_list, t_list
      
      #Initial Conditions
      u₀ = copy(x_end)
      
      # tspan = (Tspan/NFE*(i-1), Tspan/NFE*i)
      tspan = (0.0, 1.0) 
      p = u_step[i]
    
    
      alg = Vern7() 
      prob = ODEProblem(f_bioreactor!,u₀,tspan,p)
      sol = DifferentialEquations.solve(prob,  alg )

      # sol.u[end]
      # x_end_list
      # x_end_list[i+1,:]

      x_end = copy(sol.u[end])
      # x_end_list

      sol.t
      
      x_list = append!(x_list,  sol.u[:])
      t_list = append!(t_list,  (i-1).+sol.t )
    
    end

     
    x_list
    t_list

    N_list = size(t_list)[1]
      
    p1 = plot(t_list,  [x_list[i][1] for i in 1:N_list], label = "x1", ylim = (0.0, 2.0), title = "Naive")
    p1 = plot!(t_list, [x_list[i][2] for i in 1:N_list], label = "x2")

