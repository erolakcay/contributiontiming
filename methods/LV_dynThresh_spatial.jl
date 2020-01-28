## Code to run time series for the spatial Lotka-Volterra approximation with a dynamic threshold

# Parameters
Dᵤ = D # diffusion for non-delayers
Dᵥ = D # diffusion for delayers
num_sims = numberofsims # number of sims
N = sidelength # discretization (side length)
K = carryingcapacity # carrying capacity (equal across space)
tspan = (0.0,maxtime) # time span
MassMatrix = [Matrix(1.0I, 2*N^2, 2*N^2) zeros(2*N^2,N^2); zeros(N^2,3*N^2)] # mass matrix
mtime = Int(maxtime+1)
output = [zeros(mtime*num_sims,3) α*ones(mtime*num_sims,1) δ*ones(mtime*num_sims,1) zeros(mtime*num_sims,1)]

# Production functions
μᵤ = 10 # mean for the non-delayers, u
σᵤ = 2 # std dev for the non-delayers, u
μᵥ = μᵤ*(1+δ) # mean for the delayers, v
σᵥ = 2 # std dev for the delayers, v
f(t) = pdf.(Normal(μᵤ, σᵤ),t)
F(t) = cdf.(Normal(μᵤ, σᵤ),t)
g(t) = pdf.(Normal(μᵥ, σᵥ),t)
G(t) = cdf.(Normal(μᵥ, σᵥ),t)

# Define the discretized reaction-diffusion differential-algebraic equation (RDAE)
function RDAE(du,u,p,t)
    for m = 2:N-1
        for n = 2:N-1
            count = n + (m-1)*N
            du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count-N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
            du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
            du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        end
    end

    for m = 2:N-1
        # left, n=1
        count = (m-1)*N+1
        du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count-N]+u[count+N-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count+N-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        # right, n=N
        count = m*N
        du[count] = Dᵤ*(u[count+N]+u[count+1-N]+u[count-N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1-N]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    end

    for n = 2:N-1
        # top, n=1
        count = n
        du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count+(N-1)*N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
        du[N^2+count] = Dᵥ*u[count]Dᵥ*(u[count]+u[N^2+count])*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count+(N-1)*N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        # bottom, n=N
        count = n+N*(N-1)
        du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1]+u[count-N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    end

    # top left
    count = 1
    du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count+(N-1)*N]+u[count+N-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count+(N-1)*N]+u[N^2+count+N-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # top right
    count = N
    du[count] = Dᵤ*(u[count+N]+u[count+1-N]+u[count+(N-1)*N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1-N]+u[N^2+count+(N-1)*N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # bottom left
    count = N*(N-1)+1
    du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1]+u[count-N]+u[count+N-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count+N-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # bottom right
    count = N^2
    du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1-N]+u[count-N]+u[count-1]-4*u[count]) + (f(u[2*N^2+count])/u[2*N^2+count])*u[count]*(1 - (u[count]+u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1-N]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (g(u[2*N^2+count])/u[2*N^2+count])*u[N^2+count]*(1 - (u[count]+u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
end

for sim = 1:num_sims
    # Determine initial conditions
    U₀ = 5*rand(Beta(1,1),N,N) # Beta(0.5,0.5), Beta(1,1), Beta(1,4), Beta(4,1), Beta(4,4) # initial number not-delaying over space
    V₀ = 5*ones(N,N) - U₀ #5*rand(Beta(0.5,0.5),N,N) # initial number delaying over space
    Θ₀ = zeros(N,N)
    for m = 1:N
        for n = 1:N
            # solve for Θ₀[m,n]
            Θ₀[m,n] = find_zero(θ -> U₀[m,n]*F(θ) + V₀[m,n]*G(θ) - α*(U₀[m,n] + V₀[m,n]), (0, 50))
        end
    end
    u₀ = [reshape(U₀', N^2);reshape(V₀', N^2);reshape(Θ₀', N^2)] # initial conditions

    # Solve as an ODE with a singular mass matrix
    RDAEfunc = ODEFunction(RDAE;mass_matrix = MassMatrix)
    prob = ODEProblem(RDAEfunc,u₀,tspan)
    sol = solve(prob,Rodas4P(),saveat=1)

    # Calculate average proportion of delayers
    avg_delay = (mean(sol[N^2+1:2*N^2,:],dims=1)./(mean(sol[1:N^2,:],dims=1) + mean(sol[N^2+1:2*N^2,:],dims=1)))'

    # Calculate Moran's I for the Von Neumann neighbour.
    w = 1/4 # Assume that the weights don't vary across focal nodes
    MoransI = zeros(length(sol.t[:]),1)
    for m = 1:length(sol.t[:])
        P = reshape(sol[N^2+1:2*N^2,m],(N,N))./(reshape(sol[1:N^2,m],(N,N)) + reshape(sol[N^2+1:2*N^2,m],(N,N)))
        ΔP = P - mean(P)*ones(size(P)) # P - mean P
        M = w*circshift(ΔP,(0,1)) + w*circshift(ΔP,(0,-1)) + w*circshift(ΔP,(1,0)) + w*circshift(ΔP,(-1,0))
        MoransI[m,1] = sum(ΔP.*M)/sum(ΔP.^2)
    end

    colj = (sim-1)*mtime+1
    colk = sim*mtime
    global output[colj:colk,1:3] = [sol.t[:] avg_delay MoransI]
    global output[colj:colk,6] = sim*ones(mtime,1)
end

pat=string(paths, "/output/");
if !isdir(pat)
  mkdir(pat)
end

# outpath = string(pat,"LV_dynThresh_spatial_rout_",queue);

alpha=α;
delay=δ;
diffusion=D;

#@rput final_delay alpha delay diffusion output;
@rput alpha delay diffusion output;

R"""
#library(ggplot2)
#library(cowplot); theme_set(theme_cowplot())
#library(viridis)

#final_delay <- as.data.frame(final_delay)
#names(final_delay) <- c("X","Y","Delaying")
output <- as.data.frame(output)
names(output) <- c("time","delaying","MoransI","alpha","delta","sim")

#save(final_delay,file=paste("/home/bmorsky/wolbachia/output/",paste("finalspace_alpha", alpha*10,"delay",delay*100,sep="_"),".Rda", sep=""))
#save(avg_delay,file=paste("/home/bmorsky/wolbachia/output/",paste("timeseries_alpha", alpha*10,"delay",delay*100,sep="_"),".Rda", sep=""))
#save(MoransI,file=paste("/home/bmorsky/wolbachia/output/",paste("MoransI_alpha", alpha*10,"delay",delay*100,sep="_"),".Rda", sep=""))

save(output,file=paste("/home/bmorsky/wolbachia/output/",paste("output", alpha*10,"delay",delay*100,sep="_"),".Rda", sep=""))
#q1 <- ggplot() + scale_fill_viridis(option="magma",limits = c(0,1)) + geom_raster(data=final_delay,aes(x=X,y=Y,fill=Delaying)) + theme_void() + theme(legend.position="none")

#q2 <- ggplot() + geom_line(data=avg_delay,aes(x=Time,y=Delaying),size=1) + scale_y_continuous(limits = c(0, 1)) + ggtitle("Time series") + theme(legend.title=element_blank()) + xlab(expression(paste("Time"))) + ylab(expression(paste("Proportion delaying, ", p)))

#q3 <- ggplot() + geom_line(data=MoransI,aes(x=Time,y=MoransI),size=1) + scale_y_continuous(limits = c(-1, 1)) + ggtitle("Moran's I") + theme(legend.title=element_blank()) + xlab(expression(paste("Time"))) + ylab(expression(paste("Moran's ", I)))

#save_plot(q1,filename=paste("/home/bmorsky/wolbachia/output/",paste("alpha", alpha*10,"delay",delay*100,sep="_"),".png", sep=""))

#save_plot(q2,filename=paste("/home/bmorsky/wolbachia/output/",paste("avg_delay_alpha", alpha*10,"delay",delay*100,sep="_"),".png", sep=""))

#save_plot(q3,filename=paste("/home/bmorsky/wolbachia/output/",paste("MoransI_alpha", alpha*10,"delay",delay*100,sep="_"),".png", sep=""))
"""
