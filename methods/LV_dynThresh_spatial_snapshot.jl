## Code to run time series for the spatial Lotka-Volterra approximation with a dynamic threshold

# Parameters
Dᵤ = 0.25*D # diffusion for non-delayers
Dᵥ = 0.25*D # diffusion for delayers
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
            du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count-N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
            du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
            du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        end
    end

    for m = 2:N-1
        # left, n=1
        count = (m-1)*N+1
        du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count-N]+u[count+N-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count+N-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        # right, n=N
        count = m*N
        du[count] = Dᵤ*(u[count+N]+u[count+1-N]+u[count-N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1-N]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    end

    for n = 2:N-1
        # top, n=1
        count = n
        du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count+(N-1)*N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[N^2+count] = Dᵥ*u[count]Dᵥ*(u[count]+u[N^2+count])*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count+(N-1)*N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
        # bottom, n=N
        count = n+N*(N-1)
        du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1]+u[count-N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
        du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    end

    # top left
    count = 1
    du[count] = Dᵤ*(u[count+N]+u[count+1]+u[count+(N-1)*N]+u[count+N-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1]+u[N^2+count+(N-1)*N]+u[N^2+count+N-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # top right
    count = N
    du[count] = Dᵤ*(u[count+N]+u[count+1-N]+u[count+(N-1)*N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count+N]+u[N^2+count+1-N]+u[N^2+count+(N-1)*N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # bottom left
    count = N*(N-1)+1
    du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1]+u[count-N]+u[count+N-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1]+u[N^2+count-N]+u[N^2+count+N-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
    # bottom right
    count = N^2
    du[count] = Dᵤ*(u[count-(N-1)*N]+u[count+1-N]+u[count-N]+u[count-1]-4*u[count]) + (u[count]/u[2*N^2+count])*(f(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[N^2+count] = Dᵥ*(u[N^2+count-(N-1)*N]+u[N^2+count+1-N]+u[N^2+count-N]+u[N^2+count-1]-4*u[N^2+count]) + (u[N^2+count]/u[2*N^2+count])*(g(u[2*N^2+count]) - (f(u[2*N^2+count])*u[count]+g(u[2*N^2+count])*u[N^2+count])/K)
    du[2*N^2+count] = u[count]*F(u[2*N^2+count]) + u[N^2+count]*G(u[2*N^2+count]) - α*(u[count] + u[N^2+count])
end

# Determine initial conditions
using GaussianRandomFields

an = AnisotropicExponential([2000 0; 0 2000])
gau = Gaussian(0.3)
cov = CovarianceFunction(2,gau)
pts = 0:0.01:1
grf = GaussianRandomField(cov,CirculantEmbedding(),pts,pts)
YY = sample(grf)

@rput N
R"""
library(gstat)
xy <- expand.grid(1:N, 1:N)
names(xy) <- c("x","y")
g.dummy <- gstat(formula=z~1, locations=~x+y,dummy=T,beta=1,model=vgm(psill=0.025, range=5, model='Per'),nmax=20)
yy <- predict(g.dummy, newdata=xy, nsim=1)
YY <- matrix((yy[,3] - min(yy[,3]))/max(yy[,3]- min(yy[,3])),ncol=N)
"""
@rget YY
N = 101
x = kron(collect(range(1,stop=N,length=N)),ones(N,1))
y = kron(ones(N,1),collect(range(1,stop=N,length=N)))
YY = [x y reshape(YY,N^2,1)]

@rput YY
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)

YY <- as.data.frame(YY)
names(YY) <- c("X","Y","Init")

q <- ggplot() + scale_fill_viridis() + geom_raster(data=YY,aes(x=X,y=Y,fill=Init)) + theme_void() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + coord_equal() + labs(x=NULL, y=NULL)
ggsave(q,filename="/Users/brycemorsky/Desktop/Gauss.png", width = 10, height = 10)
"""

U₀ = 5*YY #rand(Beta(4,4),N,N) # Beta(0.5,0.5), Beta(1,1), Beta(1,4), Beta(4,1), Beta(4,4) # initial number not-delaying over space
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
sol = solve(prob,Rodas4P())

# Format output for final state
final_delay = reshape(sol[end][N^2+1:2*N^2],(N,N))./(reshape(sol[end][1:N^2],(N,N)) + reshape(sol[end][N^2+1:2*N^2],(N,N)))
x = kron(collect(range(1,stop=N,length=N)),ones(N,1))
y = kron(ones(N,1),collect(range(1,stop=N,length=N)))
final_delay = [x y reshape(final_delay,N^2,1)]

pat=string(paths, "/output/");
if !isdir(pat)
  mkdir(pat)
end

# outpath = string(pat,"LV_dynThresh_spatial_rout_",queue);

alpha=α;
delay=δ;
diffusion=D;

@rput alpha delay diffusion final_delay;

R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)

final_delay <- as.data.frame(final_delay)
names(final_delay) <- c("X","Y","Delaying")

save(final_delay,file=paste("/home/bmorsky/wolbachia/outputs/",paste("corr_snapshot_D",diffusion,"alpha",alpha*10,"delay",delay*100,sep="_"),".Rda", sep=""))

q <- ggplot() + scale_fill_viridis(option="magma",limits = c(0,1)) + geom_raster(data=final_delay,aes(x=X,y=Y,fill=Delaying)) + theme_void() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + coord_equal() + labs(x=NULL, y=NULL)
ggsave(q,filename=paste("/home/bmorsky/wolbachia/outputs/",paste("corr_snapshot_D",diffusion,"alpha", alpha*10,"delay",delay*100,sep="_"),".png", sep=""), width = 10, height = 10)
"""
