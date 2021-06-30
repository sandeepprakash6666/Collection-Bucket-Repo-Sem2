


#*Calculating P and Pdot matrices automatically

    #region #*Calculating RADAU Pdot matrices automatically
    #Code adapter from Carol

    # using Polynomials    

    # #*NCP = 1
    #     NCP = 1
    #     t0 = 0.0
    #     t1 = 1.0

    #     P0 = fromroots([t1])/((t0 - t1))
    #     P1 = fromroots([t0])/((t1 - t0))

    #     P0_dot = derivative(P0)
    #     P1_dot = derivative(P1)


    #     P_dot = zeros(NCP+1, NCP)
    #     tau = [t1]

    #     for i in 1:NCP
    #         # i = 2
    #         P_dot[1,i] = P0_dot(tau[i])
    #         P_dot[2,i] = P1_dot(tau[i])
    #     end

    #     P = [0.0; 1.0] 
    #     P_dot

    # #*NCP = 2
    #     NCP = 2
    #     t0 = 0.0
    #     t1 =  0.333333
    #     t2 = 1.0

    #     P0 = fromroots([t1, t2])/((t0 - t1)*(t0 - t2))
    #     P1 = fromroots([t0, t2])/((t1 - t0)*(t1 - t2))
    #     P2 = fromroots([t0, t1])/((t2 - t0)*(t2 - t1))

    #     P0_dot = derivative(P0)
    #     P1_dot = derivative(P1)
    #     P2_dot = derivative(P2)


    #     P_dot = zeros(NCP+1, NCP)
    #     tau = [t1, t2]

    #     for i in 1:NCP
    #         # i = 2
    #         P_dot[1,i] = P0_dot(tau[i])
    #         P_dot[2,i] = P1_dot(tau[i])
    #         P_dot[3,i] = P2_dot(tau[i])
    #     end

    #     P = [0.0; 0.0; 1.0] 
    #     P_dot
    
    # #*NCP = 3
    #     NCP = 3
    #     t0 = 0.0
    #     t1 =  0.155051
    #     t2 = 0.644949
    #     t3 = 1.0

    #     P0 = fromroots([t1, t2, t3])/((t0 - t1)*(t0 - t2)*(t0 - t3))
    #     P1 = fromroots([t0, t2, t3])/((t1 - t0)*(t1 - t2)*(t1 - t3))
    #     P2 = fromroots([t0, t1, t3])/((t2 - t0)*(t2 - t1)*(t2 - t3))
    #     P3 = fromroots([t0, t1, t2])/((t3 - t0)*(t3 - t1)*(t3 - t2))

    #     P0_dot = derivative(P0)
    #     P1_dot = derivative(P1)
    #     P2_dot = derivative(P2)
    #     P3_dot = derivative(P3)


    #     P_dot = zeros(NCP+1, NCP)
    #     tau = [t1, t2, t3]

    #     for i in 1:NCP
    #         # i = 2
    #         P_dot[1,i] = P0_dot(tau[i])
    #         P_dot[2,i] = P1_dot(tau[i])
    #         P_dot[3,i] = P2_dot(tau[i])
    #         P_dot[4,i] = P3_dot(tau[i])
    #     end

    #     P = [0.0; 0.0; 0.0; 1.0] 
    #     P_dot

    #endregion


#*P and Pdot matrices Hardcoded for efficiency
function collocation_matrix(NCP, Poly)
    
    #todo - error handling to be properly implemented. 
    # struct incomplete_code <: Exception end

    #*Radau Polynomials
    if Poly == "Radau"

        if NCP == 1
            # t[k,i] = [0, 1]

            Pdot = [
                -1.0
                1.0
                ]

            P = [
                0.0  
                1.0
                ]

        elseif NCP == 2
            #t[k,i] = [0, 0.333333, 1]

            Pdot = [
                    -2.0        2.0
                    1.5       -4.5
                    0.499999   2.5
                    ]

            P = [
                0.0  
                0.0   
                1.0
                ] 

        elseif NCP == 3
            #t[k,1] = [0, 0.155051, 0.644949, 1]

            Pdot = [
                -4.13939    1.73939   -3.0
                3.22475   -3.56784    5.53197
                1.16784    0.775255  -7.53197
                -0.253197   1.0532     5.0
                ]	

            P = [
                0.0  
                0.0   
                0.0 
                1.0
                ]

        else
            # throw(incomplete_code())
        end



    #*Legendre Polynomials
    elseif Poly == "Legendre"
        #todo - Complete matrices for Legendre
        if NCP == 1
            # throw(incomplete_code())
        elseif NCP == 2
            # throw(incomplete_code())
        elseif NCP == 3
            # throw(incomplete_code())
        else
            # throw(incomplete_code())
        end
    end






    return Pdot, P 
end



# myPdot, myP = collocation_matrix(3, "Radau")
# myPdot
# myP