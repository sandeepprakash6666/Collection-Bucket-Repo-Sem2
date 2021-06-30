
using Ipopt, JuMP
include("Collocation_Matrix.jl")
##

N_x = 3
N_u = 1

x0 = [0.0 0.0 0.0]



NFE = 30
NCP = 2

# T = 30
# dt = T/NFE


L = 300     #distance



##Forward Euler

    # @variable(m, x[1:N_x,0:NFE])
    # @variable(m, u[1:N_u,0:NFE-1])

    # @variable(m, dx[1:N_x,0:NFE-1])

    # @variable(m, dt)

    # for nfe in 0:NFE, nx in 1:N_x
    # # nx = 1
    # # nfe = 0
    # set_start_value(x[nx,nfe], 0)

    # end

    # #Initial Condition
    # @constraint(m, [nx in 1:N_x], x[nx,0] == x0[nx])



    #     #DAE Eqns
    #     @NLconstraint(m, [nfe in 0:NFE-1], dx[1,nfe] == x[2,nfe])
    #     @NLconstraint(m, [nfe in 0:NFE-1], dx[2,nfe] == u[1,nfe])
    #     @NLconstraint(m, [nfe in 0:NFE-1], dx[3,nfe] == 1)


    #     ##Inequality
    #     @NLconstraint(m, [nfe in 0:NFE-1], u[1,nfe] <= 1)
    #     @NLconstraint(m, [nfe in 0:NFE-1], u[1,nfe] >= -2)


    #     #Final Time constraint
    #     @constraint(m, x[1,NFE] == L)
    #     @constraint(m, x[2,NFE] == 0)

    #     #Forward Euler
    #     @constraint(m, [nx in 1:N_x,nfe in 0:NFE-1], x[nx,nfe+1] - x[nx,nfe] == dt*dx[nx,nfe] )


    #     @NLobjective(m, Min, x[3,NFE] )        #Used directly for centralized



##Collocation
    m = Model( with_optimizer(Ipopt.Optimizer, linear_solver = "mumps", max_cpu_time = 500.0, print_level = 5))

    @variable(m, x[1:N_x,0:NFE-1, 0:NCP])
    @variable(m, u[1:N_u,0:NFE-1])
    @variable(m, dx[1:N_x,0:NFE-1, 1:NCP])

    @variable(m, dt, start = 1)
    
    #Initial Condition
    @constraint(m, [nx in 1:N_x], x[nx,0,0] == x0[nx])



        #DAE Eqns
        @NLconstraint(m, [nfe in 0:NFE-1,ncp in 1:NCP], dx[1,nfe,ncp] == x[2,nfe,ncp])
        @NLconstraint(m, [nfe in 0:NFE-1,ncp in 1:NCP], dx[2,nfe,ncp] == u[1,nfe])
        @NLconstraint(m, [nfe in 0:NFE-1,ncp in 1:NCP], dx[3,nfe,ncp] == 1)

        ##Inequality
        @NLconstraint(m, [nfe in 0:NFE-1], u[1,nfe] <= 1)
        @NLconstraint(m, [nfe in 0:NFE-1], u[1,nfe] >= -2)

        #Final Time constraint
        @constraint(m, x[1,end,NCP] == L)
        @constraint(m, x[2,end,NCP] == 0)


        #
        @constraint(m, dt >= 0.1)

                    #region-Collocation Eqns
                    Pdotₘₐₜ, Pₘₐₜ = collocation_matrix(NCP, "Radau")

                    #Continuity Equations: Radau
                        @constraint(m, [nx in 1:N_x,nfe in 0:NFE-2], x[nx,nfe, NCP] == x[nx,nfe+1,0])

                    #Integration constraint : Only at internal Checkpoints
                        @constraint(m, [nx in 1:N_x,nfe in 0:NFE-1, ncp in 1:NCP], sum( x[nx,nfe,i]*Pdotₘₐₜ[i+1, ncp]  for i in 0:NCP   ) ==  dx[nx,nfe,ncp]*dt )
                #endregion



    @NLobjective(m, Min, x[3,end,end] )        #Used directly for centralized

##

JuMP.optimize!(m)


opt_x1 = [JuMP.value.(x[1,0,0]); JuMP.value.(x[1,:,NCP]).data]./100
opt_x2 = [JuMP.value.(x[2,0,0]); JuMP.value.(x[2,:,NCP]).data]./10
opt_x3 = [JuMP.value.(x[3,0,0]); JuMP.value.(x[3,:,NCP]).data]

opt_u = JuMP.value.(u[1,:]).data

opt_dt = JuMP.value(dt)

NFE*opt_dt
##


using Plots
plotly()

t_grid = collect(0:opt_dt:NFE*opt_dt)

plt = plot(legend = :topleft)
plot!(plt, t_grid, opt_x1)
plot!(plt, t_grid, opt_x2)
plot!(plt, t_grid[1:end-1], opt_u)

##

# raw_index(v::MOI.VariableIndex) = v.value


# JuMP.index.(u[1,2])
# JuMP.index.(m[:u][1,2])


# all_variables(m)
# JuMP.value.(all_variables(m))



1
## Checking Hessian





##
# using SparseArrays
# using LinearAlgebra

#     nlp = m
#     NLPBLOCK = JuMP.MOI.get(nlp, JuMP.MOI.NLPBlock())# get the models NLPblock through the MOI
#     NLPeval = NLPBLOCK.evaluator # get the NLP evaluator block to evaluate the gradients, jacobians and hessians of the model

#     #Hessian
#     hessStruct = MOI.hessian_lagrangian_structure(NLPeval)
#     hessLag = Array{Float64}(undef, length(hessStruct))

# hess_values = zeros(length(hessStruct))
# NLPeval.constraints
# JuMP.MOI.get(nlp, JuMP.MOI.NLPBlockDual())

# x_eval = value.(all_variables(nlp))

# σ = 1.0
# μ = ones(length(NLPeval.constraints))
# MOI.eval_hessian_lagrangian(NLPeval, hessLag, x_eval, σ, μ)

# HessianBig = sparse(Symmetric(sparse(first.(hessStruct), last.(hessStruct), hessLag), :L))



# H = sparse(
#     map(i -> i[1], hessStruct), 
#     map(i -> i[2], hessStruct),
#     hess_values,
#     length(x),
#     length(x),
# )




# m_eval = NLPEvaluator(m)
# MOI.initialize(m_eval, [:Hess])
# hess_struct = MOI.hessian_lagrangian_structure(m_eval)
# length(hess_struct)
# hess_values = zeros(length(hess_struct))
# x = value.(all_variables(m))
# σ = 1.0
# μ = ones(length(m_eval.constraints))
# MOI.eval_hessian_lagrangian(m_eval, hess_values, x, σ, μ)
# H = sparse(
#     map(i -> i[1], hess_struct), 
#     map(i -> i[2], hess_struct),
#     hess_values,
#     length(x),
#     length(x),
# )

# Array(H)

