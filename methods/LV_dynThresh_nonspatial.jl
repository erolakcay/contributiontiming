# Code to run time series for the Lotka-Volterra approximation
# with variable threshold at only one isolated patch
using DifferentialEquations, Distributions, Plots, SpecialFunctions, Sundials

# finds the roots of a function, used to solve for the threshold, θ
function my_bisection_method(Q::Function, a::Real, b::Real)

    ## redefine to be self-contained
    close_enough(x, y) = abs(x - y) < sqrt(eps())
    close_enough(x) = close_enough(x, 0)

    ## check endpoints
    if Q(a) * Q(b) > 0
      error("Function values at endpoints must have different signs: " ,Q(a) , " " , Q(b))
    end

    ## endpoints aren't already zeroes?
    if close_enough(Q(a))
      return(a)
    elseif close_enough(Q(b))
        return(b)
    end

    ## begin
    update_mid(a, b) = (a + b)/2
    mid = update_mid(a, b)
    ctr = 1

    while ( !close_enough(Q(mid)) )
        ## update a or b and then mid
        Q(a) * Q(mid) < 0 ? (b = mid) : (a = mid)
        mid = update_mid(a, b)
        ctr = ctr + 1
    end

    #println("Converged to $mid in $ctr steps")

    return(mid)
end

# parameters for the growth curves
μᵤ = 9.0 # mean growth rate for non-delaying
μᵥ = 6.0 # mean growth rate for delaying
σᵤ = 2.0 # std dev growth rate for non-delaying
σᵥ = 2.0 # std dev growth rate for delaying
ρ = 0.25 # worker to queen threshold
K = 1000.0 # carrying capacity

f = function(t) # non-delaying production curve
    return pdf.(TruncatedNormal(μᵤ, σᵤ, 0, 50),t)
end
F = function(t) # the cumulative non-delaying production
    return cdf(TruncatedNormal(μᵤ, σᵤ, 0, 50),t)
end
g = function(t) # delaying production curve
    #return pdf.(Normal(MV,sv),t)
    return pdf.(TruncatedNormal(μᵥ, σᵥ, 0, 50),t)
end
G = function(t) # the cumulative delaying production curve
    #return 0.5*(1 + erf.((t-MV)/(sv*sqrt(2))))
    return cdf(TruncatedNormal(μᵥ, σᵥ, 0, 50),t)
end

# The differentia-algebraic equation (DAE) we wish to solve
function DAE(out,du,u,p,t)
         out[1] = (f(u[3])*u[1]/u[3]) - u[1]*(u[1]+u[2])/K - du[1] # eq for uninfected
         out[2] = (g(u[3])*u[2]/u[3]) - u[2]*(u[1]+u[2])/K - du[2] # eq for infected
         out[3] = u[1]*F(u[3]) + u[2]*G(u[3]) - ρ*(u[1] + u[2]) # eq for threshold, θ
end

# Initial conditions
u₀ = [0.95, 0.05, 0.0] # initial uninfected, infected, and θ (to be solved)
q(θ) = u₀[1]*F(θ) + u₀[2]*G(θ) - ρ*(u₀[1] + u₀[2]) # algebraic condition for θ₀ given u₀[1] and u₀[2]
u₀[3] = my_bisection_method(q, 0.0, 1000.0) # solve for θ₀
du₀ = [(f(u₀[3])*u₀[1]/u₀[3]) - u₀[1]*(u₀[1]+u₀[2])/K, (g(u₀[3])*u₀[2]/u₀[3]) - u₀[2]*(u₀[1]+u₀[2])/K, 0.0] # initial condition for du
tspan = (0.0,10000.0) # time span

# Solve DAE
differential_vars = [true,true,false] # u[1] and u[2] are differential, u[3] is not
prob = DAEProblem(DAE,du₀,u₀,tspan,differential_vars=differential_vars)
sol = solve(prob,IDA())

# Plot infected, uninfected, and θ over time in three panels
plot(sol, tspan=(0, 2000), layout=(3,1), label=["Uninfected" "Infected" "Threshold \\theta"])
