# Draw bifurcation diagrams for the manuscript and supplementary information
using Distributions, Roots

# Production functions
δ = 0.25 # delay parameter, i.e. a 25% delay
μ = 10 # mean
σ = 2 # standard deviation
f(t) = pdf.(Normal(μ, σ),t) # non-delayed production curve
F(t) = cdf.(Normal(μ, σ),t) # the cumulative non-delayed production
fδ(t) = pdf.(Normal(μ*(1+δ), σ),t) # delayed production curve
Fδ(t) = cdf.(Normal(μ*(1+δ), σ),t) # the cumulative delayed production
q1(t) = f(t) - fδ(t)

# For the manuscript
θ = find_zero(q1,(0.01,40)) # find λ = 1

α₁ = (Fδ(θ)+0.001):0.001:0.999 # α > Fδ, delayers stable
α₂ = 0.001:0.001:(F(θ)-0.001) # α < F(θ), non-delayers stable
α = Fδ(θ):0.001:F(θ) # Fδ < α < F(θ), bistability

stable_p1 = [ones(length(α₁),1) α₁ ones(length(α₁),1)] # region delayers are stable
stable_p0 = [zeros(length(α₂),1) α₂ zeros(length(α₂),1)] # region non-delayers are stable
unstable = [ones(998-length(α₁),1) 0.001:0.001:Fδ(θ) 0.5*ones(998-length(α₁),1);
[(α - F(θ)*ones(length(α),1))/(Fδ(θ) - F(θ)) α 0.5*ones(length(α),1)];
zeros(998-length(α₂),1) F(θ):0.001:0.999 0.5*ones(998-length(α₂),1)]
length(0.001:0.001:Fδ(θ))
bif = [stable_p0; unstable; stable_p1]

# For the SI
𝚯 = zeros(1001)

𝚯[1] = find_zero(q1,(0.01,40))

for m = 1:1:1000
    τ = 0.15*m
    q2(t) = F(t+τ) - F(t) - Fδ(t+τ) + Fδ(t)
    𝚯[m+1] = find_zero(q2,(0.01,20+τ/2))
end

α₁ = Fδ(𝚯)
α₂ = F(𝚯)
bif_app = [0.15*collect(0:1:1000) zeros(1001) α₁ α₂ ones(1001)]

using RCall
@rput bif bif_app
R"""
library(ggplot2)
library(cowplot)
library(extrafont)
loadfonts()
library("reshape2")
theme_set(theme_cowplot())

bif <- as.data.frame(bif)
bif_app <- as.data.frame(bif_app)

p <- ggplot(bif,aes(x = V2, y = V1, group = V3)) + geom_line(aes(linetype=factor(V3))) +
        scale_linetype_manual(labels = c("stable","unstable","stable"), values=c("solid","dashed","solid")) +
        theme(legend.title = element_blank(), legend.position = "right", legend.direction = "vertical") +
        xlab(expression(paste("Total contribution per capita, ", alpha))) +
        ylab(expression(paste("Proportion delaying, ", p))) + ggtitle("Bifurcation diagram")

p_app <- ggplot(data = bif_app) +
        geom_ribbon(aes(x=V1,ymin=V3,ymax=V4), fill=rgb(244/256,157/256,55/256)) +
        geom_ribbon(aes(x=V1,ymin=V2,ymax=V3), fill=rgb(226/256,26/256,57/256)) +
        geom_ribbon(aes(x=V1,ymin=V4,ymax=V5), fill=rgb(17/256,119/256,184/256)) +
        scale_y_continuous(expand = c(0, 0), limits=c(0,1.1)) +
        scale_x_continuous(expand = c(0, 0), limits=c(0,15.5)) +
        ylab(expression(paste("Total contribution per capita, ", alpha))) +
        xlab(expression(paste("Competition duration, ", tau))) + ggtitle("Regions of monostability and bistability") +
        coord_flip()

save_plot(p,filename="/Users/brycemorsky/Desktop/bif_diagram.pdf",base_height = 4,base_width = 6)
save_plot(p_app,filename="/Users/brycemorsky/Desktop/app_alphaVtau.pdf",base_height = 4,base_width = 6)
"""
