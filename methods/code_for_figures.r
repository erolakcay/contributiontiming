library(ggplot2)
library(cowplot) 
theme_set(theme_cowplot())
library(dplyr)
library(colorspace)
# library(viridis)

# setwd("~/Desktop/D = 0.001, maxtime = 1000, numsims = 40, N = 40/")
setwd("~/Desktop/D=e-3, maxtime=1000, B(0.5,0.5)/")
# setwd("~/Desktop/D = 0.0001, maxtime = 10,000, numsims = 40, N = 40/")
f <- list.files(pattern = "Rda")
dat <- do.call("rbind",
               lapply(f, function(file){
                 load(file)
                 return(output)
               }))

dat[dat$alpha%in%c(.4,.5,.6) & dat$delta%in%c(.1,.3,.5),] %>% group_by(., alpha, delta, time) %>% summarise(meanDelaying=mean(delaying), meanMoran=mean(MoransI)) -> res

p <- ggplot(res) + 
  xlab("Time") + 
  facet_wrap("alpha", nrow=1, labeller=label_both) + 
  theme(strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        axis.line=element_blank()) +
  panel_border(size=1, colour="black") +
  scale_color_discrete_sequential(palette="BluGrn", name=bquote(delta)) +
  scale_x_continuous(breaks=c(0, max(res$time)/2, max(res$time))) +
  # scale_x_continuous(breaks=c(0, max(res$time)/2, max(res$time)), labels=c(0,5,10)) +
  scale_y_continuous(breaks=c(0,.5,1), limits=c(0,1))

# m <- p + 
#   ylab("Moran's index") +
#   geom_line(aes(x=time, y=meanMoran, col=factor(delta)), size=1.2)

d <- p +
  ylab("Delaying") +
  geom_line(aes(x=time, y=meanDelaying, col=factor(delta)), size=1.2)
save_plot(d, filename = paste("~/Desktop/average_proportion_delayers_moransi_ushaped_",max(res$time),".pdf",sep=""), base_height = 4, base_width = 10)
# 
# q <- plot_grid(d,m,ncol=1, labels="AUTO")
# save_plot(q, filename = paste("~/Desktop/average_proportion_delayers_moransi_",max(res$time),".pdf",sep=""), base_height = 6, base_width = 10)
