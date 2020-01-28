# Draw a bifurcation diagram
using Distributions, Roots

# Production functions
δ = 0.25 # delay parameter, i.e. a 25% delay
μᵤ = 10
σᵤ = 2
μᵥ = μᵤ*(1+δ)
σᵥ = 2
f = function(t) # non-delayed production curve
    return pdf(Normal(μᵤ, σᵤ),t)
end
F = function(t) # the cumulative non-delayed production
    return cdf(Normal(μᵤ, σᵤ),t)
end
g = function(t) # delayed production curve
    return pdf(Normal(μᵥ, σᵥ),t)
end
G = function(t) # the cumulative delayed production curve
    return cdf(Normal(μᵥ, σᵥ),t)
end


q(t) = f(t)-g(t)
T = find_zero(q,(0.1,50))

L = round(Integer,F(T)*1000)
stable_p0 = [zeros(L,1) 0.001:0.001:L/1000 zeros(L,1)]
stable_p1 = [ones(L,1) (1-L/1000):0.001:0.999 ones(L,1)]
unstable = [ones(1000-L,1) 0.001:0.001:(1-L/1000) 0.5*ones(1000-L,1);
(collect((1-L/1000):0.001:L/1000) - F(T)*ones(2*L-999,1))/(G(T) - F(T)) (1-L/1000):0.001:L/1000 0.5*ones(2*L-999,1);
zeros(1000-L,1) (L/1000):0.001:0.999 0.5*ones(1000-L,1)]

Output = [stable_p0; unstable; stable_p1]

# Figure 1, varying production curves
Fig1a = [collect(0:0.01:20) f(0:0.01:20) zeros(length(0:0.01:20),1); collect(0:0.01:20) g(0:0.01:20) ones(length(0:0.01:20),1)]
Fig1b = [collect(0:0.01:20) f(0:0.01:20)./(f(0:0.01:20) + g(0:0.01:20))]

using RCall
@rput Output Fig1a Fig1b;
R"""
library(ggplot2)
library(cowplot)
library(extrafont)
loadfonts()
library("reshape2")
theme_set(theme_cowplot())

Output <- as.data.frame(Output)
Fig1a <- as.data.frame(Fig1a)
Fig1b <- as.data.frame(Fig1b)

p <- ggplot(Output,aes(x = V2, y = V1, group = V3)) + geom_line(aes(linetype=factor(V3))) +
     scale_linetype_manual(labels = c("stable","unstable","stable"), values=c("solid","dashed","solid")) +
     theme(legend.title = element_blank(), legend.position = "right", legend.direction = "vertical") +
     xlab(expression(paste("Total contribution per capita, ", alpha))) + ylab(expression(paste("Proportion delaying, ", p))) + ggtitle("Bifurcation diagram")

qa <- ggplot(Fig1a,aes(x = V1, y = V2, group = V3)) + geom_line(aes(linetype=factor(V3))) + scale_linetype_manual(labels = c("non-delayed","delayed"), values=c("solid","dashed")) + theme(legend.key.width = unit(5,"cm"),legend.title = element_blank(), legend.position = "bottom", legend.direction = "vertical") + xlab(expression(paste("Time, ", t))) + ylab(expression(paste("Contribution"))) + ggtitle("Contribution functions")

qb <- ggplot(Fig1b,aes(x = V1, y = V2)) + geom_line() + theme(legend.title = element_blank()) + xlab(expression(paste("Time, ", t))) + ylab(expression(paste("Contribution"))) + ggtitle("Contribution functions")

save_plot(p,filename="/Users/brycemorsky/Desktop/bif_diag.pdf",base_height = 4,base_width = 6)
save_plot(qa,filename="/Users/brycemorsky/Desktop/contribution_functions_a.png",base_height = 3,base_width = 6)
save_plot(qb,filename="/Users/brycemorsky/Desktop/contribution_functions_b.png",base_height = 3,base_width = 6)
"""
