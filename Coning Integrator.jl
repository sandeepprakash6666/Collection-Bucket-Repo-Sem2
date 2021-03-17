using DifferentialEquations
using Sundials
using Plots

##*
N_x = 15

W = 50          #m
dx = W/N_x 


function f_coning!(du,u,p,t)
        # du  = zeros(N_x)
        # u   = ones(N_x)*10
        # p   = 0
        # t   = 1



    h = u
    q0 = p


   #m


    x_grid    = collect(1 : 1 : N_x)
    indx_grid = collect(1 : 1 : N_x) 



    # h = ones(1+ NFE_x+1 +1)
    # h = [1; 2; 3; 4; 4]  

    
    #*Calculate using finite differences
    dh  = zeros(N_x)
    d²h = zeros(N_x)

    #LHS
    h_L = h[2] - q0*2*dx/(2*h[1])
    dh[1]  = (h[2] - h_L)           /(2*dx) 
    d²h[1] = (h[2] - 2*h[1] + h_L)  /(dx^2) 

    #middle
    for i in 2:(N_x-1)
        # i = 2
        # global dh, d²h      

        dh[i]  = (h[i+1]  - h[i-1])             /(2*dx)
        d²h[i] = (h[i+1] - 2*h[i] + h[i-1])    /(dx^2)

    end

    #RHS
    h_R      = h[N_x-1] 
    dh[N_x]  = (h_R - h[N_x - 1])               /(2*dx) 
    d²h[N_x] = (h_R - 2*h[N_x] + h[N_x - 1])    /(dx^2)


    # dh
    # d²h

    #*
     for i  in 1:N_x                   #Differentials wrt time
        
        du[i] = dh[i]^2 + h[i]*d²h[i]
    end

    # du
end


#Initial State (alg state must be feasible)
u₀ = ones(N_x)*10
du₀ = 0*copy(u₀)
tspan = (0.0,500.0)
p = 1


alg = Vern7() 
prob = ODEProblem(f_coning!,u₀,tspan,p)
sol = solve(prob,  alg )


# sol.u

# sol.t
# size(sol.t)[1]
# sol.u[1]
# [sol.u[i][1] for i in 1 : size(sol.t)[1]]

# plot(sol.t, [sol.u[i][1] for i in 1:size(sol.t)[1]])
# plot!(sol.t, [sol.u[i][2] for i in 1:size(sol.t)[1]])


x_grid = range(0, stop = W, length = N_x)
@gif for i = 1:size(sol.t)[1]
    plot(x_grid, sol.u[i], xlim = (-1,51), ylim = (0,11))
end



