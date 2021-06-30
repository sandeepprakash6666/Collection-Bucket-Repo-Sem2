using JuMP, Ipopt
##

#Continous time LTI system matirces
    Aₜ = [0 1; 
        -2 -3]

    Bₜ = [0; 1]

    Cₜ = [1 0]

    Dₜ = [0]

#Discrete LTI System Matrices
    A = [0.9909 0.0861; 
        -0.1722 0.7326]

    B = [0.0045; 0.0861]

    C = [1 0]

    D = [0]

N_x = 2
N_u = 1
N_y = 1

x00 = [0.4; 0]
(T0, Tf) = (0, 5)
dt = 0.1        #Fixed because discrete matrix made at this discretization
NFE = (Tf - T0)/dt

##Auxiliary Functions
function Discrete_input(from)
    input_res = 0.2
    floor(from/input_res)*input_res
end
function Sim_discrete_Plant(x0, u0)
    #Dummy variables for Debugging
        # x0 = [0.35 ;0]
        # u0 = [0.7]

        x1 = A*x0 + B.*u0
end
function Build_OCP(m, (Mₚ,Nₚ), dt)
    #Dummy variables for debugging
        # m = Model( with_optimizer(Ipopt.Optimizer))
        # (Mₚ,Nₚ) = (0,5)    #(0,5)
        # dt = 0.1
    #

    #Declaring Variables
    @variable(m, x[1:N_x,Mₚ:Nₚ])   
    @variable(m, u[1:N_u,Mₚ:Nₚ-1]) 

    #Initial Condition  
    @NLparameter(m,  x̄0[1:N_x] == 0.2)     
        m[:x0] = x̄0      
    @NLconstraint(m, Fix_Initial[nx in 1:N_x], x[nx,Mₚ] == x̄0[nx])

    #Difference Eqns
    @constraint(m, [nfe in Mₚ:Nₚ-1], x[1,nfe+1] == sum( A[1,nx]*x[nx,nfe]   for nx in 1:N_x) + B[1]*u[1,nfe]  )
    @constraint(m, [nfe in Mₚ:Nₚ-1], x[2,nfe+1] == sum( A[2,nx]*x[nx,nfe]   for nx in 1:N_x) + B[2]*u[1,nfe]  )

#Objective Term
@expression(m, SP_Obj, dt*sum(  2*(x[1,nfe+1] - 0.35)^2 + (x[2,nfe+1] - 0.0)^2 + (u[1,nfe] - 0.7)^2 for nfe in Mₚ:Nₚ-1))
# @NLexpression(m, SP_Obj, dt*sum(  (x[1,nfe+1] - 8)^2 + (x[2,nfe+1] - 4)^2 + (u[1,nfe] - 1)^2 for nfe in Mₚ:Nₚ-1))
@objective(m, Min, SP_Obj )        #Used directly for centralized
end


OCP = Model( with_optimizer(Ipopt.Optimizer, linear_solver = "mumps", max_cpu_time = 500.0, print_level = 4))

Build_OCP(OCP, (0, NFE), dt)
for nx in 1:N_x #Set Initial State
    set_value(OCP[:x0][nx], x00[nx])      #need 1+ because overlap defined as 0:NFE, but in array 1:NFE+1
end

##
xk = x00

            #for Plotting
                plot_t_plant = [0.0] 
                plot_x_plant = [x00]
                plot_u_plant = []
                plot_u_opt   = []  

##
for iter_MPC in 1:1
    JuMP.optimize!(OCP)
    # plot_OCP_sol(OCP)

    uk_opt = JuMP.value.(OCP[:u])[1,0]      #Output from MPC

    uk_plant = Discrete_input(uk_opt)       

    xk1 = Sim_discrete_Plant(xk, uk_plant)  #Simulating Plant


    #Updates for Next Iterarion
        for nx in 1:N_x
            set_value(OCP[:x0][nx], xk1[nx])  
        end

        xk = xk1

        #For Plotting
            push!(plot_t_plant, plot_t_plant[end] + dt)
            push!(plot_x_plant, xk1)
            push!(plot_u_plant, uk_plant)
            push!(plot_u_opt,   uk_opt)
end

plot_OCP_sol(OCP)

plot_closed_loop()


##Plotting
using Plots
plotly()

function plot_OCP_sol(m)

    t_grid = collect(T0 : dt: Tf)


    opt_x1 = JuMP.value.(m[:x])[1,:].data
    opt_x2 = JuMP.value.(m[:x])[2,:].data
    opt_u1 = JuMP.value.(m[:u])[1,:].data

    plt = plot(title = "OpenLoop", ylim = (-0.1, 1.1))
    plot!(plt, t_grid, opt_x1,  label = "x1")
    plot!(plt, t_grid, opt_x2,  label = "x2")
    plot!(plt, t_grid, [opt_u1; opt_u1[end]], label = "u1",    linetype = :steppost)

end

function plot_closed_loop()

    x1_plant = [plot_x_plant[i][1] for i in 1:size(plot_x_plant)[1] ]
    x2_plant = [plot_x_plant[i][2] for i in 1:size(plot_x_plant)[1] ]
    u1_plant = [plot_u_plant[i][1] for i in 1:size(plot_u_plant)[1] ]
    u1_opt   = [plot_u_opt[i][1]   for i in 1:size(plot_u_opt)[1]   ]


    plt = plot(title = "ClosedLoop", ylim = (-0.1, 1.1))
    plot!(plt, plot_t_plant, x1_plant, label = "x1_cl")
    plot!(plt, plot_t_plant, x2_plant, label = "x2_cl")
    plot!(plt, plot_t_plant, [u1_plant; u1_plant[end]], label = "u1_cl",        linetype = :steppost)
    plot!(plt, plot_t_plant, [u1_opt; u1_opt[end]],     label = "u1_MPC_cl",    linetype = :steppost, linestyle = :dash)

end



##
A
A - ones(2,2)

-B.*0.7

inv(A - ones(2,2))

inv(A - ones(2,2))*(    B.*0.7   )
inv(ones(2,2) - A)*(    B.*0.7   )


(1-0.9909)*(0.00861/0.01722)

0.00861/0.01722