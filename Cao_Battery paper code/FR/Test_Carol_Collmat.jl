
import DifferentialEquations
using JuMP, Ipopt, Polynomials

# c = 3 # No of col points
    # Dmin = 0
    # Dmax = 1
    # mumax = 0.4
    # x1max = 4.5
    # np = hor[1]
    # nm = hor[2]
    # R = 0.5
    # Q = 1
    # h = dt
    # #Collocation Points Using Radau Roots 3rd degree polynomial
    # t0 = 0
    # t1 = 0.155051
    # t2 = 0.644949
    # t3 = 1.000000
    # # Lagrange polynomials
    # l0 = fromroots([t1,t2,t3])/((t0-t1)*(t0-t2)*(t0-t3))
    # l1 = fromroots([t0,t2,t3])/((t1-t0)*(t1-t2)*(t1-t3))
    # l2 = fromroots([t0,t1,t3])/((t2-t0)*(t2-t1)*(t2-t3))
    # l3 = fromroots([t0,t1,t2])/((t3-t0)*(t3-t1)*(t3-t2))
    # # 1st derivatives
    # dl0 = derivative(l0)
    # dl1 = derivative(l1)
    # dl2 = derivative(l2)
    # dl3 = derivative(l3)
    # # Collocation matrix: 1st derivatives evaluated at the
    # # collocation points
    # adot = zeros(4,4)
    # tau = [t0,t1,t2,t3]
    # for i = 1:c+1
    # adot[1,i] = dl0(tau[i])
    # adot[2,i] = dl1(tau[i])
    # adot[3,i] = dl2(tau[i])
    # adot[4,i] = dl3(tau[i])
    # end





# if c == 3
    #Collocation Points Using Radau Roots 3rd degree polynomial
    t0 = 0
    t1 = 0.155051
    t2 = 0.644949
    t3 = 1.000000
    # Lagrange polynomials
    l0 = fromroots([t1,t2,t3])/((t0-t1)*(t0-t2)*(t0-t3))
    l1 = fromroots([t0,t2,t3])/((t1-t0)*(t1-t2)*(t1-t3))
    l2 = fromroots([t0,t1,t3])/((t2-t0)*(t2-t1)*(t2-t3))
    l3 = fromroots([t0,t1,t2])/((t3-t0)*(t3-t1)*(t3-t2))
    # 1st derivatives
    dl0 = derivative(l0)
    dl1 = derivative(l1)
    dl2 = derivative(l2)
    dl3 = derivative(l3)
    # Collocation matrix: 1st derivatives evaluated at the
    # collocation points
    adot = zeros(4,4)
    tau = [t0,t1,t2,t3]
    for i = 1:c+1
    adot[1,i] = dl0(tau[i])
    adot[2,i] = dl1(tau[i])
    adot[3,i] = dl2(tau[i])
    adot[4,i] = dl3(tau[i])
    end


adot
    # return adot