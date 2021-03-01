using DifferentialEquations
using Sundials
using Plots
using Interpolations



#Single Parameter

p = [1.0]
function f!(du,u,p,t)

  du[1] = (1.2 - u[2]^2)*u[1] - u[2] + p[1]
  du[2] = u[1]    

  end

#Initial Conditions
u₀ = [0.01, 1.0]
tspan = (0.0,10.0)

prob = ODEProblem(f!,u₀,tspan,p)

sol = solve(prob)
plot(sol)



#2 Inputs

A = [1.0, 10.0]

function my_hold(t)
  
  u = [1.0; 5*ones(5); 4*ones(4) ]
  n = round(t)
  Int_n = convert(Int, n) 


  return(u[Int_n])
end

# itp = interpolate((A,), A, Gridded(Linear()))
# itp = interpolate((A,), A, Gridded(Constant()))
# itp
# itp(4.5)
# itp(8.5)

# p2 = [1.0, 0.5]

function f2!(du,u,p,t)

# param = itp(t)
param = my_hold(t)

du[1] = (1.2 - u[2]^2)*u[1] - u[2] + param
du[2] = u[1]  

end

#Initial Conditions
u₀ = [0.01, 1.0]
tspan = (1.0,10.0)

prob = ODEProblem(f2!,u₀,tspan,p2)
sol = solve(prob)
plot(sol)

