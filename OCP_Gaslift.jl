# Oil well with gas lift, model equations & parameters taken from:
# Krishnamoorthy, D., Fjalestad, K. and Skogestad, S., 2019. Optimal operation of oil and gas production using simple feedback control structures. Control Engineering Practice, 91, p.104107.
# https://doi.org/10.1016/j.conengprac.2019.104107
# taken to Julia by evren, 1 well version

# using DifferentialEquations
using JuMP
using Ipopt

using Parameters
using Plots
plotlyjs()

include("Collocation_Matrix.jl")


##* Model Parameters
@with_kw struct WellPar
    # L, H, D = Length, Height, Diameter [m]
    
    # well
    L_w = 1500      # length well
    H_w = 1000      # height well
    D_w = 0.121     # diameter well
    
    # bottom hole
    L_bh = 500 
    H_bh = 500
    D_bh = 0.121
    
    # annulus
    L_a = L_w
    H_a = H_w
    D_a = 0.189

    rho_o  = 800     # density of oil, [kg/m3]
    mu_oil = 1*0.001

    C_iv = 1e-3     # injection valve characteristic,[m2]
    C_pc = 2e-3     # choke valve characteristic, [m2]

    GOR0 = 0.1      # Gas Oil Ratio
    PI0 = 2.2       # Productivity index, [kg/s/bar]
    
    Press_r = 150   # reservoir pressure [bar]
    p_m = 20        # manifold pressure  [bar]
    T_a = 28  + 273 # [K]
    T_w = 32  + 273 # [K]

    Mw = 20e-3      # ? units [g/mol]


end
# change parameters from the default here
par = WellPar()

# and unpacking all the parameters
@unpack_WellPar par

const A_bh  = π*D_bh^2/4
const LA_bh = L_bh*A_bh
const A_w   = π*D_w^2/4
const LwAw  = L_w*A_w
const g     = 9.81      # m/s2
const R     = 8.314     # J/mol K
const V_a   = L_a*(π*(D_a/2)^2 - π*(D_w/2)^2) 


#* Model Initial states and Guesses

#Initial state and Guesses
# x0 = [1.32; 0.8; 6.0]    

# dx_guess = 0*copy(x0)
# x_guess = copy(x0)
# # z_guess = [77.0; 47.0; 62.0; 240.0; 38.0; 4.4; 34.0; 61.0; 100.0; 1.0; 34.0; 3.4]
# u_guess = [1.0]

    #Initial State and Guesses - from SS soln
    x0 = [ 3.6902822709988854  1.9050192017979886  0.20234621601298794]
    dx_guess = [0.0 0.0 0.0]
    x_guess = [ 3.6902822709988854  1.9050192017979886  0.20234621601298794]
    z_guess = [216.23  100.199  1.7281  0.86515  52.6818  47.6234  5.05843  87.7611  127.007  47.1175  5.05843  5.05843] #From SS soln
    u_guess = [47.117526523816565]  




x_lb = [10e-3; 10e-3; 10e-3]
x_ub = [10e7;  10e7;  10e7]

z_lb = [1e-1; 1e-1;     1e-2; 1e-2;     1e-2; 1e-2; 1e-2;   0.1; 30.0;      1e-2; 1e-2; 1e-2]
z_ub = [150e4; 70e4;    900e4; 900e4;   50e4; 50e4; 50e4;   150e4; 150e4;   50e4; 50e4; 50e4]

u_lb = [5e-1]
u_ub = [50e4]

dx0 = 0*copy(x_guess)


#* OCP Parameters
T0  =  0.0  
Tf  =  1.0                      #*hrs
NFE =  30
NCP =  3

dt  =  (Tf - T0)*3600/NFE       #*sec

Nx = 3
Nu = 1
Nz = 12

#Steady State - Constant Disturbance
GOR = GOR0*ones(1,NFE) 
# PI  = PI0*ones( 1,NFE) 
    
    #Step Changes
    # GOR = hcat(     0.1*ones(1, convert(Int32, NFE/3)),      0.15*ones(1,convert(Int32, NFE/3)),     0.15*ones(1,convert(Int32, NFE/3)))
    PI  = hcat(     2.2*ones(1, convert(Int32, NFE/3)),      2.4*ones(1,convert(Int32, NFE/3)),      2.4*ones(1,convert(Int32, NFE/3)))


##
find_a_SS = 0
if NFE == 1 && NCP == 1
    global find_a_SS = 1
end
    
Sim_Mode = 1    #*To simulate using IpOPT 


    ##* Defining Solver
    model1 = Model(with_optimizer(Ipopt.Optimizer))
    set_optimizer_attributes(model1, "max_iter" => 1000)

    #region-> #*Define Variables and Objective

        #region-> Defining all variables
            ## Declare general OCP Variables
            @variable(model1, x[1:Nx, 1:NFE, 1:NCP])
            # @variable(model1, z[1:Nz, 1:NFE, 1:NCP])
            # @variable(model1, u[1:Nu, 1:NFE])
                                
            #unscaled variables
            # @variable(model1, q[    1:1,        1:NFE, 1:NCP])
            @variable(model1, dx_us[1:Nx,       1:NFE, 1:NCP])
            # @variable(model1, dq[   1:1,        1:NFE, 1:NCP])

            #region-> Auxiliary Variables - For writing model Equations clearly
                @variables(model1, begin
                #Differential Variables
                    x_lb[1] <=  mass_ga[nfe in 1:NFE, ncp in 1:NCP]     <=  x_ub[1]     , (start = x_guess[1]) # mass of gas in annulus
                    x_lb[2] <=  mass_gt[nfe in 1:NFE, ncp in 1:NCP]     <=  x_ub[2]     , (start = x_guess[2])# mass of gas in tubing
                    x_lb[3] <=  mass_ot[nfe in 1:NFE, ncp in 1:NCP]     <=  x_ub[3]     , (start = x_guess[3])# mass of oil in tubing
                
                #Algebraic variables
                    z_lb[1] <=  p_ai[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[1]      , (start = z_guess[1])
                    z_lb[2] <=  p_wh[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[2]      , (start = z_guess[2])

                    z_lb[3] <=  rho_ai[nfe in 1:NFE, ncp in 1:NCP]      <= z_ub[3]      , (start = z_guess[3])
                    z_lb[4] <=  rho_m[nfe  in 1:NFE, ncp in 1:NCP]      <= z_ub[4]      , (start = z_guess[4])
                    
                    z_lb[5] <=  w_pc[nfe   in 1:NFE, ncp in 1:NCP]      <= z_ub[5]      , (start = z_guess[5])# total flow through choke
                    z_lb[6] <=  w_pg[nfe   in 1:NFE, ncp in 1:NCP]      <= z_ub[6]      , (start = z_guess[6])# produced gas rate
                    z_lb[7] <=  w_po[nfe   in 1:NFE, ncp in 1:NCP]      <= z_ub[7]      , (start = z_guess[7])
                                
                    # pressure at well injection point
                    z_lb[8] <=  p_wi[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[8]      , (start = z_guess[8])
                    z_lb[9] <=  p_bh[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[9]      , (start = z_guess[9])
                    
                    z_lb[10] <= w_iv[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[10]     , (start = z_guess[10]) 
                    z_lb[11] <= w_ro[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[11]     , (start = z_guess[11])
                    z_lb[12] <= w_rg[nfe in 1:NFE, ncp in 1:NCP]        <= z_ub[12]     , (start = z_guess[12])

                    # Manipulated Inputs
                    u_lb[1] <=  w_gl[nfe in 1:NFE]                      <= u_ub[1]      , (start = u_guess[1])
                end)
            #endregion

            #region-> Mapping the auxilliary variables to general OCP variable vectors
                @constraints(model1, begin
                    #Differential states
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_ga[nfe, ncp]   == x[1,nfe,ncp]     # mass of gas in annulus
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_gt[nfe, ncp]   == x[2,nfe,ncp]     # mass of gas in tubing
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_ot[nfe, ncp]   == x[3,nfe,ncp]     # mass of oil in tubing                
                end)

            #endregion
        #endregion
    
        ## Objective
        if Sim_Mode == 0
            if find_a_SS == 0
                @NLobjective(model1, Min,  sum(    w_gl[nfe] -  w_po[nfe,1]      for nfe in 1:NFE )  )
            end
        end

    #endregion

    #region-> #*Define Constraints

        if find_a_SS == 0
            #fixing initial Point
            @NLconstraints(model1, begin
                Constr_x0[nx in 1:Nx],  x[nx,1,1] == x0[nx]
            end)
        end

        #Fixing the inputs if simulation mode
        if Sim_Mode == 1
            for nfe in 1:NFE
                fix(w_gl[nfe],  u_guess[1] ; force = true)
            end
        end


        #Defining the model ODEs in each line
        @NLconstraints(model1, begin
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx_us[1, nfe, ncp]      ==  (w_gl[nfe]      - w_iv[nfe, ncp])*1e-3 
            Constr_ODE2[nfe in 1:NFE, ncp in 1:NCP], dx_us[2, nfe, ncp]      ==  (w_iv[nfe, ncp] + w_rg[nfe, ncp]*1e-1 - w_pg[nfe, ncp])*1e-3 
            Constr_ODE3[nfe in 1:NFE, ncp in 1:NCP], dx_us[3, nfe, ncp]      ==  (w_ro[nfe, ncp] - w_po[nfe, ncp])*1e-3 
        end)

        #Defining Model Algebraic Equations in each line
        @NLconstraints(model1, begin
            #
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP],    p_ai[nfe,ncp]   == 1e-5*(((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga[nfe,ncp]*1e3) + (Mw/(R*T_a)*((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga[nfe,ncp]*1e3))*9.81*H_a)                         # annulus pressure    
            Constr_Alg2[nfe in 1:NFE, ncp in 1:NCP],    p_wh[nfe,ncp]   == 1e-5*(((R*T_w/Mw)*(mass_gt[nfe,ncp]*1e3/(L_w*A_w + L_bh*A_bh - mass_ot[nfe,ncp]*1e3/rho_o))) - ((mass_gt[nfe,ncp]*1e3+mass_ot[nfe,ncp]*1e3 )/(L_w*A_w))*9.81*H_w/2)  # wellhead pressure
            #
            Constr_Alg3[nfe in 1:NFE, ncp in 1:NCP],    rho_ai[nfe,ncp] == 1e-2*(Mw/(R*T_a)*p_ai[nfe,ncp]*1e5)  # gas, in annulus
            Constr_Alg4[nfe in 1:NFE, ncp in 1:NCP],    1e2*rho_m[nfe,ncp]*(mass_ot[nfe,ncp]*1e3*p_wh[nfe,ncp]*1e5*Mw + rho_o*R*T_w*mass_gt[nfe,ncp]*1e3)  == ((mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)*p_wh[nfe,ncp]*1e5*Mw*rho_o)     # mixture, in tubing
            #
            Constr_Alg5[nfe in 1:NFE, ncp in 1:NCP],    w_pc[nfe,ncp]^2     == C_pc^2*(rho_m[nfe,ncp]*1e2*(p_wh[nfe,ncp]*1e5 - p_m*1e5))                          # total flow through choke
            Constr_Alg6[nfe in 1:NFE, ncp in 1:NCP],    w_pg[nfe,ncp]*(mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)   == mass_gt[nfe,ncp]*1e3*w_pc[nfe,ncp]   # produced gas rate
            Constr_Alg7[nfe in 1:NFE, ncp in 1:NCP],    w_po[nfe,ncp]*(mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)   == mass_ot[nfe,ncp]*1e3*w_pc[nfe,ncp]   # produced oil rate
            #
            Constr_Alg8[nfe  in 1:NFE, ncp in 1:NCP],    p_wi[nfe,ncp]  == 1e-5*((p_wh[nfe,ncp]*1e5 + 9.81/(A_w*L_w)*(mass_ot[nfe,ncp]*1e3 + mass_gt[nfe,ncp]*1e3-rho_o*L_bh*A_bh)*H_w + 128*mu_oil*L_w*w_pc[nfe,ncp] /  (3.141*D_w^4*((mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)*p_wh[nfe,ncp]*1e5*Mw*rho_o) /  (mass_ot[nfe,ncp]*1e3*p_wh[nfe,ncp]*1e5*Mw + rho_o*R*T_w*mass_gt[nfe,ncp]*1e3)  )  ))
            Constr_Alg9[nfe  in 1:NFE, ncp in 1:NCP],    p_bh[nfe,ncp]  == 1e-5*(p_wi[nfe,ncp]*1e5 + rho_o*9.81*H_bh + 128*mu_oil*L_bh*w_po[nfe,ncp]/(3.14*D_bh^4*rho_o)) 
            #
            Constr_Alg10[nfe in 1:NFE, ncp in 1:NCP],    w_iv[nfe,ncp]^2  == C_iv^2*(rho_ai[nfe,ncp]*1e2*(p_ai[nfe,ncp]*1e5 - p_wi[nfe,ncp]*1e5))
            Constr_Alg11[nfe in 1:NFE, ncp in 1:NCP],    w_ro[nfe,ncp]  == PI[nfe]*1e-6*(Press_r*1e5 - p_bh[nfe,ncp]*1e5) 
            Constr_Alg12[nfe in 1:NFE, ncp in 1:NCP],    w_rg[nfe,ncp]  == 1e1*GOR[nfe]*w_ro[nfe,ncp] 

        end)

        #Additional constraints-Because of eqn reformulations
        @NLconstraints(model1, begin
            Constr_Alg5_1[ nfe in 1:NFE, ncp in 1:NCP],  (rho_m[nfe,ncp]*1e2*(p_wh[nfe,ncp]*1e5  - p_m*1e5))             >= 0
            Constr_Alg10_1[nfe in 1:NFE, ncp in 1:NCP],  (rho_ai[nfe,ncp]*1e2*(p_ai[nfe,ncp]*1e5 - p_wi[nfe,ncp]*1e5))   >= 0
        end)
    
        if isSS == 0
            #region-> #generic code -> Collocation Equation for Differential Equations AND Objective Function (No-scaling)
                            ## Creating a Radau collocation Matrix for NCP = 3
                            Pdotₘₐₜ, Pₘₐₜ = collocation_matrix(3, "Radau")
                            @NLconstraints(model1, begin
                                #Collocation Continuity
                                Constr_Coll_Cont_Diff[nx in 1:Nx, nfe in 1:NFE-1, ncp in 1:1], sum(Pₘₐₜ[i]*x[nx, nfe, i] for i in 1:NCP) == x[nx, nfe+1, ncp]
                                # Constr_Coll_Cont_quad[            nfe in 1:NFE-1, ncp in 1:1], sum(Pₘₐₜ[i]*q[1,  nfe, i] for i in 1:NCP) == q[1,  nfe+1, ncp]

                                #Collocation Integration
                                Constr_Coll_Int_Diff[nx in 1:Nx, nfe in 1:NFE, ncp in 2:NCP], sum(Pdotₘₐₜ[ncp-1, i]*x[nx, nfe, i] for i in 1:NCP)    == dt*dx_us[nx, nfe, ncp]
                                # Constr_Coll_Int_quad[            nfe in 1:NFE, ncp in 2:NCP], sum(Pdotₘₐₜ[ncp-1, i]*q[1, nfe, i]  for i in 1:NCP)    == dt*dq[1, nfe, ncp]

                            end)
                    #endregion
        else
            @NLconstraints(model1, begin
                Constr_SS[nx in 1:Nx, nfe in 1:NFE, ncp in 1:NCP], dx_us[nx,nfe,ncp] == 0    
            end)
              
        end
    #endregion

##* Solve Model 


model1

optimize!(model1)
JuMP.termination_status(model1)
JuMP.solve_time(model1::Model)



JuMP.value.(mass_ga[:, NCP])
JuMP.value.(mass_gt[:, NCP])
JuMP.value.(mass_ot[:, NCP])



if isSS == 1
    dx_guess = JuMP.value.(dx_us[:, NCP])
    x_guess = JuMP.value.(x[:,NFE,NCP])
    u_guess = JuMP.value.(w_gl[NFE,NCP])
    z_guess = hcat(
                    JuMP.value.(p_ai[NFE,NCP]),
                    JuMP.value.(p_wh[NFE, NCP]),

                    JuMP.value.(rho_ai[NFE, NCP]),
                    JuMP.value.(rho_m[NFE, NCP]), 
                                        
                    JuMP.value.(w_pc[NFE, NCP]),      
                    JuMP.value.(w_pg[NFE, NCP]),   
                    JuMP.value.(w_po[NFE, NCP]),      

                    JuMP.value.(p_wi[NFE, NCP]),        
                    JuMP.value.(p_bh[NFE, NCP]),        
                                        
                    JuMP.value.(w_iv[NFE, NCP]),       
                    JuMP.value.(w_ro[NFE, NCP]),       
                    JuMP.value.(w_rg[NFE, NCP])  
                   )

end



# star_w_po = JuMP.value.(w_po[:, NCP])
# star_w_gl = JuMP.value.(w_gl[:])



##* Plotting Solution of OCP

t_plot = collect(T0:  dt/3600:  Tf)    #Returns NFE+1 dimensional vector

p11 = plot(t_plot[1:end-1], GOR[:],                                 label = "GOR",  linetype = :steppost, linestyle = :dash)

p12 = plot(t_plot[1:end-1], JuMP.value.(w_gl[:]),                   label = "w_gl", linetype = :steppost)

p13 = plot(t_plot[1:end-1], JuMP.value.(w_po[:, NCP]),              label = "w_po")
p13 = plot!(t_plot[1:end-1], JuMP.value.(w_pg[:, NCP]),             label = "w_pg")

p14 = plot(t_plot[1:end-1], JuMP.value.(mass_ot[:, NCP]),           label = "m_ot")
p14 = plot!(t_plot[1:end-1], JuMP.value.(mass_gt[:, NCP]),          label = "m_gt")


##

p11
p12
p13
p14




