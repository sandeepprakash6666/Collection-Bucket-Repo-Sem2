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


##* Model constants
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
    L_a =  L_w
    H_a = H_w
    D_a = 0.189

    rho_o = 800     # density of oil, kg/m3
    mu_oil = 1*0.001

    C_iv = 1e-3     # injection valve characteristic, m2
    C_pc = 2e-3     # choke valve characteristic, m2

    GOR = 0.1       # Gas Oil Ratio
    PI = 2.2        # Productivity index, kg/s/bar
    
    Press_r = 150   # reservoir pressure
    p_m = 20        # manifold pressure
    T_a = 28  + 273 # K
    T_w = 32  + 273 # K

    Mw = 20e-3 # ? units g/mol


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


#* OCP Parameters
T0  =  0.0
Tf  =  30.0
dt  =  1000
NFE =  2
NCP =  3

Nx = 3
Nu = 1
Nz = 12

x0 = [1.32; 0.8; 6.0]   #initial state 

x_guess = copy(x0)
z_guess = [77.0; 47.0; 62.0; 240.0; 38.0; 4.4; 34.0; 61.0; 100.0; 1.0; 34.0; 3.4]
u_guess = [1.0]
dx0 = 0*x_guess

##
    ##* Defining Solver
    model1 = Model(with_optimizer(Ipopt.Optimizer))

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

            #Auxiliary Variables - For writing model Equations
                @variables(model1, begin
                #Differential Variables
                    mass_ga[nfe in 1:NFE, ncp in 1:NCP]     , (start = x_guess[1]) # mass of gas in annulus
                    mass_gt[nfe in 1:NFE, ncp in 1:NCP]     , (start = x_guess[2])# mass of gas in tubing
                    mass_ot[nfe in 1:NFE, ncp in 1:NCP]     , (start = x_guess[3])# mass of oil in tubing
                
                #Algebraic variables
                    p_ai[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[1])
                    p_wh[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[2])

                    rho_ai[nfe in 1:NFE, ncp in 1:NCP]      , (start = z_guess[3])
                    rho_m[nfe  in 1:NFE, ncp in 1:NCP]      , (start = z_guess[4])
                    
                    w_pc[nfe   in 1:NFE, ncp in 1:NCP]      , (start = z_guess[5])# total flow through choke
                    w_pg[nfe   in 1:NFE, ncp in 1:NCP]      , (start = z_guess[6])# produced gas rate
                    w_po[nfe   in 1:NFE, ncp in 1:NCP]      , (start = z_guess[7])
                                
                    # pressure at well injection point
                    p_wi[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[8])
                    p_bh[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[9])
                    
                    w_iv[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[10]) 
                    w_ro[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[11])
                    w_rg[nfe in 1:NFE, ncp in 1:NCP]        , (start = z_guess[12])

                # Manipulated Inputs
                    w_gl[nfe in 1:NFE]                      , (start = u_guess[1])
            end)

            #Set Initial Guesses for variables
            for nx in 1:Nx, nfe in 1:NFE, ncp in 1:NCP
                set_start_value(x[nx, nfe, ncp],    x_guess[nx])
                # set_start_value(z[nz, nfe, ncp],    z_guess[nz])
                set_start_value(dx_us[nx, nfe, ncp],dx0[nx])
                # set_start_value(u[nu, nfe],         u_guess[nu])
                # set_start_value(q[1, nfe, ncp],     q0)
                # set_start_value(dq[1, nfe, ncp],    dq0)


          end



            #region-> Mapping the auxilliary variables to general OCP variable vectors
                @constraints(model1, begin
                    #Differential states
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_ga[nfe, ncp]   == x[1,nfe,ncp]     # mass of gas in annulus
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_gt[nfe, ncp]   == x[2,nfe,ncp]     # mass of gas in tubing
                        [nfe in 1:NFE, ncp in 1:NCP],    mass_ot[nfe, ncp]   == x[3,nfe,ncp]     # mass of oil in tubing
                    
                    #Algebraic states
                    #     [nfe in 1:NFE, ncp in 1:NCP],   p_ai[nfe, ncp]      == z[1,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   p_wh[nfe, ncp]      == z[2,nfe,ncp]

                    #     [nfe in 1:NFE, ncp in 1:NCP],   rho_ai[nfe, ncp]    == z[3,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   rho_m[nfe, ncp]     == z[4,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_pc[nfe, ncp]      == z[5,nfe,ncp]     # total flow through choke
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_pg[nfe, ncp]      == z[6,nfe,ncp]     # produced gas rate
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_po[nfe, ncp]      == z[7,nfe,ncp]
                                        
                    #     # pressure at well injection point
                    #     [nfe in 1:NFE, ncp in 1:NCP],   p_wi[nfe , ncp]     == z[8,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   p_bh[nfe, ncp]      == z[9,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_iv[nfe, ncp]      == z[10,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_ro[nfe, ncp]      == z[11,nfe,ncp]
                    #     [nfe in 1:NFE, ncp in 1:NCP],   w_rg[nfe , ncp]     == z[12,nfe,ncp]

                    # # Manipulated Inputs
                    #     [nfe in 1:NFE],                 w_gl[nfe]           == u[1,nfe]
                
                end)

            #endregion
        #endregion
    
        ## Objective
        @NLobjective(model1, Min,  sum(    w_gl[nfe] -  w_po[nfe,1]      for nfe in 1:NFE )  )
    
    #endregion

    nfe = 1 
    ncp = 1

##
    #region-> #*Define Constraints

        #fixing initial Point
        @NLconstraints(model1, begin
            Constr_x0[nx in 1:Nx],  x[nx,1,1] == x0[nx]
        end)

        #Defining the model ODEs in each line
        @NLconstraints(model1, begin
            Constr_ODE1[nfe in 1:NFE, ncp in 1:NCP], dx_us[1, nfe, ncp]      ==  (w_gl[nfe]      - w_iv[nfe, ncp])*1e-3 
            Constr_ODE2[nfe in 1:NFE, ncp in 1:NCP], dx_us[2, nfe, ncp]      ==  (w_iv[nfe, ncp] + w_rg[nfe, ncp]*1e-1 - w_pg[nfe, ncp])*1e-3 
            Constr_ODE3[nfe in 1:NFE, ncp in 1:NCP], dx_us[3, nfe, ncp]      ==  (w_ro[nfe, ncp] - w_po[nfe, ncp])*1e-3 
        end)

        #Defining Model Algebraic Equations in each line
        @NLconstraints(model1, begin
            #
            Constr_Alg1[nfe in 1:NFE, ncp in 1:NCP],    p_ai[nfe,ncp]   == 1e-5#*(((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga[nfe,ncp]*1e3) + (Mw/(R*T_a)*((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga[nfe,ncp]*1e3))*9.81*H_a)                         # annulus pressure    
            Constr_Alg2[nfe in 1:NFE, ncp in 1:NCP],    p_wh[nfe,ncp]   == 1e-5#*(((R*T_w/Mw)*(mass_gt[nfe,ncp]*1e3/(L_w*A_w + L_bh*A_bh - mass_ot[nfe,ncp]*1e3/rho_o))) - ((mass_gt[nfe,ncp]*1e3+mass_ot[nfe,ncp]*1e3 )/(L_w*A_w))*9.81*H_w/2)  # wellhead pressure
            #
            Constr_Alg3[nfe in 1:NFE, ncp in 1:NCP],    rho_ai[nfe,ncp] == 1e-2#*(Mw/(R*T_a)*p_ai[nfe,ncp]*1e5)  # gas, in annulus
            Constr_Alg4[nfe in 1:NFE, ncp in 1:NCP],    rho_m[nfe,ncp]  == 1e-2#*(((mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)*p_wh[nfe,ncp]*1e5*Mw*rho_o)/(mass_ot[nfe,ncp]*1e3*p_wh[nfe,ncp]*1e5*Mw + rho_o*R*T_w*mass_gt[nfe,ncp]*1e3))     # mixture, in tubing
            #
            Constr_Alg5[nfe in 1:NFE, ncp in 1:NCP],    w_pc[nfe,ncp]^2   == C_pc^2*(rho_m[nfe,ncp]*1e2*(p_wh[nfe,ncp]*1e5 - p_m*1e5))                          # total flow through choke
            Constr_Alg6[nfe in 1:NFE, ncp in 1:NCP],    w_pg[nfe,ncp]*(mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)   == mass_gt[nfe,ncp]*1e3*w_pc[nfe,ncp]   # produced gas rate
            Constr_Alg7[nfe in 1:NFE, ncp in 1:NCP],    w_po[nfe,ncp]*(mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)   == mass_ot[nfe,ncp]*1e3*w_pc[nfe,ncp]   # produced oil rate
            #
            Constr_Alg8[nfe  in 1:NFE, ncp in 1:NCP],    p_wi[nfe,ncp]  == 1e-5#*((p_wh[nfe,ncp]*1e5 + 9.81/(A_w*L_w)*(mass_ot[nfe,ncp]*1e3 + mass_gt[nfe,ncp]*1e3-rho_o*L_bh*A_bh)*H_w + 128*mu_oil*L_w*w_pc[nfe,ncp] /  (3.141*D_w^4*((mass_gt[nfe,ncp]*1e3 + mass_ot[nfe,ncp]*1e3)*p_wh[nfe,ncp]*1e5*Mw*rho_o) /  (mass_ot[nfe,ncp]*1e3*p_wh[nfe,ncp]*1e5*Mw + rho_o*R*T_w*mass_gt[nfe,ncp]*1e3)  )  ))
            Constr_Alg9[nfe  in 1:NFE, ncp in 1:NCP],    p_bh[nfe,ncp]  == 1e-5#*(p_wi[nfe,ncp]*1e5 + rho_o*9.81*H_bh + 128*mu_oil*L_bh*w_po[nfe,ncp]/(3.14*D_bh^4*rho_o)) 
            #
            Constr_Alg10[nfe in 1:NFE, ncp in 1:NCP],    w_iv[nfe,ncp]  == C_iv#*sqrt(rho_ai[nfe,ncp]*1e2*(p_ai[nfe,ncp]*1e5 - p_wi[nfe,ncp]*1e5))
            Constr_Alg11[nfe in 1:NFE, ncp in 1:NCP],    w_ro[nfe,ncp]  == PI#*1e-6*(Press_r*1e5 - p_bh[nfe,ncp]*1e5) 
            Constr_Alg12[nfe in 1:NFE, ncp in 1:NCP],    w_rg[nfe,ncp]  == 1e1#*GOR*w_ro[nfe,ncp] 

        end)
    

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


    #endregion


##* Solve Model

model1

optimize!(model1)
JuMP.termination_status(model1)
JuMP.solve_time(model1::Model)







# function gas_well!(dx,x,p,t)

    #* mass flow rates
    mass_ga = x[1] # mass of gas in annulus
    mass_gt = x[2] # mass of gas in tubing
    mass_ot = x[3] # mass of oil in tubing
    
    #* Parameters
    w_gl = p
    _GOR = GOR
    #* NL calculations

    #Pressures
    p_ai = 1e-5*(((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga*1e3) + (Mw/(R*T_a)*((R*T_a/(V_a*Mw) + 9.81*H_a/V_a)*mass_ga*1e3))*9.81*H_a) # annulus pressure
    p_wh = 1e-5*(((R*T_w/Mw)*(mass_gt*1e3/(L_w*A_w + L_bh*A_bh - mass_ot*1e3/rho_o))) - ((mass_gt*1e3+mass_ot*1e3 )/(L_w*A_w))*9.81*H_w/2) # wellhead pressure

    # densities
    rho_ai = 1e-2*(Mw/(R*T_a)*p_ai*1e5)  # gas, in annulus
    rho_m = 1e-2*(((mass_gt*1e3 + mass_ot*1e3)*p_wh*1e5*Mw*rho_o)/(mass_ot*1e3*p_wh*1e5*Mw + rho_o*R*T_w*mass_gt*1e3))  # mixture, in tubing

    w_pc = C_pc*sqrt(rho_m*1e2*(p_wh*1e5 - p_m*1e5)) # total flow through choke
    w_pg = (mass_gt*1e3/(mass_gt*1e3+mass_ot*1e3))*w_pc # produced gas rate
    w_po = (mass_ot*1e3/(mass_gt*1e3+mass_ot*1e3))*w_pc # produced oil rate


    # pressure at well injection point
    p_wi = 1e-5*((p_wh*1e5 + 9.81/(A_w*L_w)*(mass_ot*1e3+mass_gt*1e3-rho_o*L_bh*A_bh)*H_w + 128*mu_oil*L_w*w_pc/(3.141*D_w^4*((mass_gt*1e3 + mass_ot*1e3)*p_wh*1e5*Mw*rho_o)/(mass_ot*1e3*p_wh*1e5*Mw + rho_o*R*T_w*mass_gt*1e3)))) 
    p_bh = 1e-5*(p_wi*1e5 + rho_o*9.81*H_bh + 128*mu_oil*L_bh*w_po/(3.14*D_bh^4*rho_o)) 
    w_iv = C_iv*sqrt(rho_ai*1e2*(p_ai*1e5 - p_wi*1e5)) 
    w_ro = (PI)*1e-6*(Press_r*1e5 - p_bh*1e5) 
    w_rg = 1e1*_GOR*w_ro 

    #* Differential Equations
    dx[1] = (w_gl - w_iv)*1e-3 
    dx[2] = (w_iv + w_rg*1e-1 - w_pg)*1e-3 
    dx[3] = (w_ro - w_po)*1e-3 


# end

#* Initial conditions

x0 = zeros(3)
# 
x0[1:3]= [1.32, 0.8, 6.0]

#* Parameters
w_gl = 4.0 # kg/s

p = w_gl # GOR can also be an uncertain parameter

tspan = (0.0, 2000.0)
prob = ODEProblem(gas_well!, x0, tspan, p)
sol = solve(prob)

plot(sol, label = ["gas annulus" "gas tubing" "oil tubing"])



