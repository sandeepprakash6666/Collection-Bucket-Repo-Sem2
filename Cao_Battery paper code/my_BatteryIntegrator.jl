using Parameters
using DifferentialEquations
using Ipopt
using JuMP
using Interpolations
using Sundials
using Plots

##*Parameters

    @with_kw struct BatteryParam



        
        # A. Number of node points
        N1 = 20     #?
        N2 = 20     #?
        
        Ncp = 2     #? number of concentration at the positive side
        Ncn = 2     #? number of concentration at the positive side
        Nsei = 3    #?
        Ncum = 4    #?

    
        #=
        B. Parameters
        1: cathode, 2: anode
        Dp,Dn: Solid phase diffusivity (m2/s),
        kp, kn: Rate constant for lithium intercalation reaction (m2.5/(mol0.5s))
        cspmax,csnmax: Maximum solid phase concentration at positive (mol/m3),
        lp,ln : Region thickness (m),
        ap,an: Particel surface area to volume (m2/m3)
        Rpp ,Rpn: particle radius at positive (m),
        ce: Electrolyte Concentration (mol/m3)
        M[sei]: Molecular weight of SEI (Kg/mol)
        Kappa[sei]: SEI ionic conductivity (S/m)
        rho[sei]: SEI density (Kg/m3)
        ksei: rate constant of side reaction (C m/s mol)
        F : Faraday's constant (C/mol) , R : Ideal gas constant (J/K/mol), T : Temperature (K) 
        =#

        F = 96487
        R = 8.3143
        T = 298.15
        M_sei = 0.073
        Kappa_sei = 5e-6
        rho_sei = 2.1e3
        ksei = 1.5e-12                            # 1e-8
        Urefs = 0.4
        Rsei = 0.01

        
        area    = 1  # 0.3108
        cspmax  = 10350
        csnmax  = 29480
        lp      = 6.521e-5
        lnn     = 2.885e-5
        Rpp     = 1.637e-7
        Rpn     = 3.596e-6
        ce      = 1042
        ep      = 1 - 0.52
        en      = 1 - 0.619
        ap      = 3 * ep / Rpp
        an      = 3 * en / Rpn
        Sp      = area * lp * ap
        Sn      = area * lnn * an
        kp      = 1.127e-7 / F
        kn      = 8.696e-7 / F
        Dn      = 8.256e-14
        Dp      = 1.736e-14
        TC      = 2.3 / 0.3108
        Qmax    = TC
        P_nominal = TC * 3.1


    end

    #Change parameters from the default here
     Battery_par = BatteryParam()

    #Unpacking all the parameters
     @unpack_BatteryParam Battery_par

    #Specifying if variable is Differential or Algebraic
        differential_vars = trues(Ncp + Ncn + 4 + Nsei + Ncum)
        differential_vars[Ncp]              = false
        differential_vars[Ncp+Ncn]          = false
        differential_vars[Ncp+Ncn+1]        = false
        differential_vars[Ncp+Ncn+2]        = false
        differential_vars[Ncp+Ncn+3]        = false
        differential_vars[Ncp+Ncn+4]        = false
        differential_vars[Ncp+Ncn+5]        = false
        differential_vars[Ncp+Ncn+6]        = false
    differential_vars

    #*Initial Differential state of Battery
     soc = 0.6
        csp_avg0 = cspmax * 1 - csnmax * soc * lnn * en / lp / ep            # 49503.111
        csn_avg0 = csnmax * soc
     delta_sei0 = 1e-10

        u0  = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        du0 = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        u0[1:Ncp]               .= csp_avg0
        u0[(Ncp+1):(Ncp+Ncn)]   .= csn_avg0
        u0[Ncp+Ncn + 7]         = delta_sei0    # delta_sei
        u0[Ncp+Ncn + 8]         = 0             # cm
        u0[Ncp+Ncn + 9]         = 0             # cp
        u0[Ncp+Ncn + 10]        = 0             # Q
        u0[Ncp+Ncn + 11]        = 0             # cf
    u0
    du0

    #todo Time and Power Profile Change Manually 
    tspan = (0.0, 1*3600.0)
    TIME_FR_segment = 0:2:tspan[2]
    
    P_FR_segment             = NaN*ones(1801)
    P_FR_segment[1:900]     .= 20
    P_FR_segment[901:end]  .= -20
    # P_FR_segment[1201:end]  .= 0
    P_FR_segment

     #Linear Interpolation for parameter profile
     itp = interpolate((TIME_FR_segment,), P_FR_segment, Gridded(Linear()))


##*DAE Function for Battery Electrochemistry Equations
function f_common(out, du, u, p, t)
    #some default values for debugging
        # du  = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        # # u   = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        # u   = [5173.778657414506, 5173.778657414506, 14740.0, 14740.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-10, 0.0, 0.0, 0.0, 0.0]
        # out = zeros(Ncp + Ncn + 4 + Nsei + Ncum)

    #Unpacking variables

            csp         = u[1 : Ncp]
            csn         = u[Ncp+1 : Ncp+Ncn]

        csp_avg     = csp[1]
        csp_s       = csp[2]
        csn_avg     = csn[1]
        csn_s       = csn[2] 

        iint        = u[Ncp+Ncn + 1]
        phi_p       = u[Ncp+Ncn + 2] 
        phi_n       = u[Ncp+Ncn + 3]
        pot         = u[Ncp+Ncn + 4]
        
        it          = u[Ncp+Ncn + 5]
        isei        = u[Ncp+Ncn + 6]
        delta_sei   = u[Ncp+Ncn + 7] 

        cm          = u[Ncp+Ncn + 8]
        cp          = u[Ncp+Ncn + 9]
        Q           = u[Ncp+Ncn + 10]
        cf          = u[Ncp+Ncn + 11]    

    #Calculating Intermediate Variables
        #Positive
            theta_p = csp_s/cspmax

            Up = 7.49983 - 13.7758 * theta_p .^ 0.5 + 21.7683 * theta_p - 12.6985 * theta_p .^ 1.5 +
                0.0174967 ./ theta_p - 0.41649 * theta_p .^ (-0.5) -
                0.0161404 * exp.(100 * theta_p - 97.1069) +
                0.363031  * tanh.(5.89493 * theta_p - 4.21921)


            jp = 2 * kp * ce^(0.5) * 
                (cspmax - csp[Ncp])^(0.5) * csp[Ncp]^(0.5) *
                sinh(0.5 * F / R / T * (phi_p - Up)) 

        #Negative
            theta_n = csn_s/csnmax

            Un = 9.99877 - 9.99961 * theta_n .^ 0.5 - 9.98836 * theta_n + 8.2024 * theta_n .^ 1.5 +
                0.23584 ./ theta_n - 2.03569 * theta_n .^ (-0.5) -
                1.47266 * exp.(-1.14872 * theta_n + 2.13185) -
                9.9989  * tanh.(0.60345 * theta_n - 1.58171)

            jn = 2 *kn *ce^(0.5) *
                (csnmax - csn[Ncn])^(0.5) * csn[Ncn]^(0.5) *
                sinh(0.5 * F / R / T * (phi_n - Un + (Rsei + delta_sei / Kappa_sei) * it / an / lnn))

    #*C1.Governing Equations 
        #Positive 
        out[1] = -3*jp/Rpp - du[1]
        out[2] = 5 * (csp_s - csp_avg) + Rpp * jp / Dp

        #Negative
        out[Ncp+1] = -3 * jn / Rpn - du[Ncp+1]
        out[Ncp+2] = 5 * (csn_s - csn_avg) + Rpn * jn / Dn

    #*C2.Additional Equaions
        out[Ncp+Ncn+1] = (jp - (it   / ap / F / lp))      #?
        out[Ncp+Ncn+2] = (jn + (iint / an / F / lnn))
        out[Ncp+Ncn+3] = pot - phi_p + phi_n

    #*C3.SEI Layer Equations
        out[Ncp+Ncn+4] = - iint + it - isei      #?
        out[Ncp+Ncn+5] = - isei 
                        + an*lnn *ksei*exp( -1 * F / R / T * 
                                                (phi_n - Urefs  + it / an / lnn * (delta_sei / Kappa_sei + Rsei)))
        
        out[Ncp+Ncn+6] = isei * M_sei / F / rho_sei / an / lnn - du[Ncp+Ncn+7]   # d delta_sei/dt

    #*C4.Charge Stored
        out[Ncp+Ncn+7]      = iint     / 3600 - du[Ncp+Ncn+8]       # dcm/dt
        out[Ncp+Ncn+5+Nsei] = it * pot / 3600 - du[Ncp+Ncn+9]       # dcp/dt
        out[Ncp+Ncn+6+Nsei] = it       / 3600 - du[Ncp+Ncn+10]      # dQ/dt
        out[Ncp+Ncn+7+Nsei] = isei     / 3600 - du[Ncp+Ncn+11]      # dcf/dt

end

##*Function to get feasible initial Algebraic states
function getinitial(csp_avg, csn_avg, delta_sei, power)
    #Some default values for debugging
        # csp_avg   = 5173.778657414506
        # csn_avg   = 14740.0
        # delta_sei = 1.0e-10
        # power     = 20.0

    #Calculating Initial Guesses for variables
        theta_p_guess = min(0.9, csp_avg / cspmax)
        theta_n_guess = min(0.9, csn_avg / csnmax)

        Un_guess = 9.99877 - 9.99961 * theta_n_guess .^ 0.5 - 9.98836 * theta_n_guess +
                   8.2024 * theta_n_guess .^ 1.5 +
                   0.23584 ./ theta_n_guess - 2.03569 * theta_n_guess .^ (-0.5) -
                   1.47266 * exp.(-1.14872 * theta_n_guess + 2.13185) -
                   9.9989 * tanh.(0.60345 * theta_n_guess - 1.58171)
        
        Up_guess = 7.49983 - 13.7758 * theta_p_guess .^ 0.5 + 21.7683 * theta_p_guess -
                   12.6985 * theta_p_guess .^ 1.5 + 0.0174967 ./ theta_p_guess -
                   0.41649 * theta_p_guess .^ (-0.5) -
                   0.0161404 * exp.(100 * theta_p_guess - 97.1069) +
                   0.363031 * tanh.(5.89493 * theta_p_guess - 4.21921)

        pot0 = max(min(Up_guess - Un_guess, 3.3), 2.0)
        it0 = power/pot0
        isei0 = an *lnn *ksei * exp( -1 * F / R / T *(Un_guess - Urefs + delta_sei / Kappa_sei * it0 / an / lnn))

    #*Model for Optimization
     m = Model(with_optimizer(Ipopt.Optimizer))
     #Variables
        #Algebraic variables
        @variable(m,        csp_s,              start = csp_avg)
        @variable(m,        csn_s,              start = csn_avg)

        @variable(m,        iint,               start = it0)
        @variable(m,        phi_p,              start = Up_guess)
        @variable(m,        phi_n,              start = Un_guess)

        @variable(m,        it,                 start = it0)
        @variable(m,        isei,               start = isei0)

        #Intermediate Variables
        @variable(m, 0 <=   theta_p     <= 1,   start = theta_p_guess)
        @variable(m, 0 <=   theta_n     <= 1,   start = theta_n_guess)
        @variable(m,        Up,                 start = Up_guess)
        @variable(m,        Un,                 start = Un_guess)

     #Constraints

        #Intermediate Variables
            @constraint(m, theta_p * cspmax == csp_s)
            @constraint(m, theta_n * csnmax == csn_s)

            @NLconstraint(m, Up == 7.49983 - 13.7758 * theta_p^0.5 + 21.7683 * theta_p - 12.6985 * theta_p^1.5 +
                                   0.0174967 / theta_p - 0.41649 * theta_p^(-0.5) -
                                   0.0161404 * exp(100 * theta_p - 97.1069) +
                                   0.363031  * tanh(5.89493 * theta_p - 4.21921)
                        )

            @NLconstraint(m, Un == 9.99877 - 9.99961 * theta_n^0.5 - 9.98836 * theta_n + 8.2024 * theta_n^1.5 +
                                   0.23584 / theta_n - 2.03569 * theta_n^(-0.5) -
                                   1.47266 * exp(-1.14872 * theta_n + 2.13185) -
                                   9.9989 * tanh(0.60345 * theta_n - 1.58171)
                        )

            @NLconstraint(m, (cspmax - csp_s)^(0.5) * csp_s^(0.5) * sinh(0.5 * F / R / T * (phi_p - Up)) - 
                              (it / ap / F / lp / (2 * kp * ce^(0.5))) == 0
                        )   

            @NLconstraint(m, (csnmax - csn_s)^(0.5) * csn_s^(0.5) * sinh( 0.5 * F / R / T * (phi_n - Un  + (delta_sei / Kappa_sei + Rsei) * it / an / lnn),) + 
                              iint / an / F / lnn / (2 * kn * ce^(0.5)) == 0
                        )

        #Governing Equations
            @constraint(m, 5 * (csp_s - csp_avg) + Rpp * it   / F / Dp / ap / lp  == 0)
            @constraint(m, 5 * (csn_s - csn_avg) - Rpn * iint / F / Dn / an / lnn == 0)
        
        #SEI Layer Equations     
            @constraint(m, -iint + it - isei == 0)   
            @NLconstraint(m, 1e4 * (-isei + 
                                     an *lnn *ksei *exp( -1 * F / R / T * (phi_n - Urefs + (delta_sei / Kappa_sei + Rsei) * it / an / lnn))
                                   ) == 0
                         )

        #External connection Equations
            @constraint(m, it * (phi_p - phi_n) == power)
        
    #*Solve Model
        JuMP.optimize!(m)

        csp_s0  = getvalue(csp_s)
        csn_s0  = getvalue(csn_s)
        iint0   = getvalue(iint) 
        phi_p0  = getvalue(phi_p)
        phi_n0  = getvalue(phi_n)
        pot0    = phi_p0 - phi_n0
        it0     = getvalue(it)  
        isei0   = getvalue(isei)

    return csp_s0, csn_s0, iint0, phi_p0, phi_n0, pot0, it0, isei0

end  

##*DAE Function for Battery System Equations
function f_FR(out, du, u, param, t)
    #some default values for debugging
        # du  = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        # # u   = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        # u   = [5173.778657414506, 5173.778657414506, 14740.0, 14740.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-10, 0.0, 0.0, 0.0, 0.0]
        # out = zeros(Ncp + Ncn + 4 + Nsei + Ncum)
        # t = 1800.5

    #Unpacking variables
     p_to_battery = itp[t]
    # p_to_battery = power_to_battery0
    phi_p  = u[Ncp+Ncn + 2]
    phi_n  = u[Ncp+Ncn + 3]
    it     = u[Ncp+Ncn + 5]

    #Electrochemistry Equations
    f_common(out, du, u, p_to_battery, t)
    out

    #Battery System Equations
    out[end] = (phi_p - phi_n) * it - p_to_battery

end

##*Simulate Battery using Integrator
    
    #Get feasible algebraic states for the current differential states
    csp_s0, csn_s0, iint0, phi_p0, phi_n0, pot0, it0, isei0 = 
        getinitial(csp_avg0, csn_avg0, delta_sei0, P_FR_segment[1])

    #Collecting initial algebraic states into vector
        u0[Ncp]         = csp_s0
        u0[Ncp+Ncn]     = csn_s0
        u0[Ncp+Ncn+1]   = iint0                  # iint
        u0[Ncp+Ncn+2]   = phi_p0                 # phi_p
        u0[Ncp+Ncn+3]   = phi_n0                 # phi_n
        u0[Ncp+Ncn+4]   = pot0                   # pot
        u0[Ncp+Ncn+5]   = it0                    # it
        u0[Ncp+Ncn+6]   = isei0                  # isei
        
    u0

##  

    prob = DAEProblem(f_FR, du0, u0, tspan, differential_vars = differential_vars)
    sol  = DifferentialEquations.solve(prob, IDA())

##*Plotting everything
    plotlyjs()
    N_plot = size(sol.t)[1]

    [sol.u[i][1] for i in 1:N_plot]
    p1 = Plots.plot(sol.t,  [sol.u[i][1]  for i in 1:N_plot], label = "csp_avg")
    p1 = Plots.plot!(sol.t, [sol.u[i][2]  for i in 1:N_plot], label = "csp_s")
    p1 = Plots.plot!(sol.t, [sol.u[i][3]  for i in 1:N_plot], label = "csn_avg")
    p1 = Plots.plot!(sol.t, [sol.u[i][4]  for i in 1:N_plot], label = "csn_s")

    p2 = Plots.plot(sol.t,  [sol.u[i][5]  for i in 1:N_plot], label = "I_int")
    p2 = Plots.plot!(sol.t, [sol.u[i][6]  for i in 1:N_plot], label = "phi_p")
    p2 = Plots.plot!(sol.t, [sol.u[i][7]  for i in 1:N_plot], label = "phi_n")
    p2 = Plots.plot!(sol.t, [sol.u[i][8]  for i in 1:N_plot], label = "V or pot")
    p2 = Plots.plot!(sol.t, [sol.u[i][9]  for i in 1:N_plot], label = "It")
    p2 = Plots.plot!(sol.t, [sol.u[i][10] for i in 1:N_plot], label = "Isei")

    p3 = Plots.plot(sol.t,  [sol.u[i][11] for i in 1:N_plot], label = "delta_sei")
    p3 = Plots.plot!(sol.t, [sol.u[i][15] for i in 1:N_plot], label = "Cf")


    p4 = Plots.plot(sol.t,  [sol.u[i][12] for i in 1:N_plot], label = "cm")
    p4 = Plots.plot!(sol.t, [sol.u[i][13] for i in 1:N_plot], label = "cp")
    p4 = Plots.plot!(sol.t, [sol.u[i][14] for i in 1:N_plot], label = "Q")


    p1
##
    p2
    p3
    p4



##

E =  NaN*ones(1801)
E[1] = 0

for i in 1:1800

    E[i+1] = E[i] + P_FR_segment[i]*2/3600 

end


rang

t_plot = collect(range(0, stop = 3600, length = 1801))

plot!(t_plot, E)

E
P_FR_segment