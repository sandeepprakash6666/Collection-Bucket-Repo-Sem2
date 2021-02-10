


function collocation_matrix(NCP, Poly)
    
    # NCP = 3
    # Poly = "Radau"

    #todo - error handling properly implemented. 
    # struct incomplete_code <: Exception end



    #*Radau Polynomials
    if Poly == "Radau"

        if NCP == 2
            # t[k,i] = [0, 1]

            Pdot = [ -1.0  1.0]

            P = [0  1]

        elseif NCP == 3
            #t[k,i] = [0, 0.333333, 1]

            Pdot = [   -2.0   1.5   0.499999		
                        2.0  -4.5   2.5		]

            P = [0  0   1] 

        elseif NCP == 4
            #t[k,1] = [0, 0.155051, 0.644949, 1]

            Pdot = [    -4.1394     3.2247      1.1678      -0.2532		
                        1.7394      -3.5678     0.7753      1.0532	
                        -3.0000     5.5320      -7.5320     5.0000  ]	

            P = [-1.78E-15  0   8.88E-16 1]

        # elseif NCP = 5 #!not completed

        #     Pdot = [    -7.1556     5.6441      1.9235      -0.5859     0.1739
        #                 2.5082      -5.0492     1.2211      1.7547      -0.4348
        #                 -1.9649     3.4925      -3.9845     0.6348      1.8221
        #                 4.0000      -6.9235     6.5953      -12.1718    8.5000  ]	

        #     P = [0.0    0.0   0.0   0.0   1.0]   

        else
            # throw(incomplete_code())
        end



    #*Legendre Polynomials
    elseif Poly == "Legendre"
        #todo - Complete matrices for Legendre
        if NCP == 2
            # throw(incomplete_code())
        elseif NCP == 3
            # throw(incomplete_code())
        elseif NCP == 4
            # throw(incomplete_code())
        else
            # throw(incomplete_code())
        end
    end






    return Pdot, P 
end



# myPdot, myP = collocation_matrix(4, "Radau")
# myPdot
# myP