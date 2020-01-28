# # Code to run time series for the spatial Lotka-Volterra approximation
# # with a fixed threshold
using DifferentialEquations, Distributions, LinearAlgebra, Plots, RecursiveArrayTools, SpecialFunctions, Sundials

# parameters for the growth curves
μᵤ = 9.0 # mean growth rate for non-delaying
μᵥ = 6.0 # mean growth rate for delaying
σᵤ = 2.0 # std dev growth rate for non-delaying
σᵥ = 2.0 # std dev growth rate for delaying
θ = 6 # time at which the next population is produced
f = 100*pdf.(TruncatedNormal(μᵤ, σᵤ, 0, 50),θ)/θ # growth rate for non-delaying
g = 100*pdf.(TruncatedNormal(μᵥ, σᵥ, 0, 50),θ)/θ # growthe rate for delaying
# note that I multiplied the above by 100, since they would be very low
# otherwise (10^-2)
K = 1000.0 # carrying capacity

# parameters and variables for space
D = 0.001 # diffusion parameter
N = 100 # discretization
X = reshape([i for i in 1:100 for j in 1:100],N,N) # length
Y = reshape([j for i in 1:100 for j in 1:100],N,N) # width

# Laplacian operator
Mx = Matrix(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
My = copy(Mx)
# Do the reflections
Mx[2,1] = 2.0
Mx[end-1,end] = 2.0
My[1,2] = 2.0
My[end,end-1] = 2.0

# Define the initial condition as normal arrays
U = 10*rand(N,N); V = 10*rand(N,N);
u0 = ArrayPartition((U,V))

# Define the discretized reaction-diffusion equation (RDE)
function RDE(du,u,p,t)
  U,V = u.x
  dU,dV = du.x
  DU = D*(My*U.*(U+V) + (U.*(U+V))*Mx)
  DV = D*(My*(V.*(U+V)) + (V.*(U+V))*Mx)
  @. dU = DU + f*U - U.*(U+V)/K
  @. dV = DV + g*V - V.*(U+V)/K
end

# Solve the discretized RDE
prob = ODEProblem(RDE,u0,(0.0,1000.0))
sol = solve(prob,BS3(),progress=true,save_everystep=false,save_start=false)
# Plot the frequency of uninfected over space. Given that θ is fixed, one
# or the other type will be dominant; there will be no bistability.
heatmap(sol[end].x[1]./(sol[end].x[1] + sol[end].x[2]))
