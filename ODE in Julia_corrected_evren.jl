# Oil well with gas lift, model equations & parameters taken from:
# Krishnamoorthy, D., Fjalestad, K. and Skogestad, S., 2019. Optimal operation of oil and gas production using simple feedback control structures. Control Engineering Practice, 91, p.104107.
# https://doi.org/10.1016/j.conengprac.2019.104107
# taken to Julia by evren, 1 well version

using DifferentialEquations
using Parameters
using Plots
plotlyjs()

## Declaring some constants
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
    Press_r = 150   # reservoir pressure
    PI = 2.2        # Productivity index, kg/s/bar
    p_m = 20        # manifold pressure
    T_a = 28 + 273  # K
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


##
function gas_well!(dx,x,p,t)

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
end

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



