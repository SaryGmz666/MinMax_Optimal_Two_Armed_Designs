library(ggplot2)
library(reshape2)

ytext <- 25
xtext <- 18 # xtext <- 25
legtext <- 20
psize <- 8
xtitle <- 20
ytitle <- 20
legtitle <- 20
plot.theme <- theme(axis.text.y  = element_text(size=ytext, colour = 'black'), axis.text.x  = element_text(size=xtext, colour = 'black'),
                    axis.title.y  = element_text(size=ytitle, vjust=0.35), axis.title.x  = element_text(size=xtitle,vjust=0),
                    panel.background = element_rect(fill = "white", colour = 'black'),
                    legend.text = element_text(size = legtext), legend.title = element_text(size=legtitle),
                    plot.title = element_text(lineheight=.8, face="bold", size = 20,hjust = 0.5),
                    strip.text = element_text(size=20), strip.text.x = element_text(size=20, face="bold"), 
                    #legend.position="none", 
                    panel.grid.major = element_line(size = 0.5, linetype = 'dashed',
                                                    colour = "grey"), 
                    panel.grid.minor = element_line(size = 0.25, linetype = 'dashed',
                                                    colour = "grey"))

fds.plot <- function(X, H, Zs){
  Sigmainvs <- list()
  Sigma.approx <- list()
  for(j in 1:ncol(X)) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    Sigmainvs[[j]] <- solve(Sigma)
    Sigma.approx[[j]] <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
  }
  zs <- matrix(0, 0, ncol(X))
  for(i in 1:nrow(Zs)) {
    z <- Zs[i, ]
    zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
  }
  
  zs.sort <- matrix(0, ncol = ncol(zs), nrow = nrow(zs))
  for (i in 1:ncol(zs.sort)){
    zs.sort[,i] <- sort(zs[,i])
  }
  colnames(zs.sort) <- colnames(X)
  dat.plot <- melt(zs.sort)
  fds <- ggplot(data = dat.plot, aes(x = Var1, y = value, group = Var2, color = Var2)) + geom_line()
  fds <- fds + xlab("") + ylab("Scaled prediction variance") + plot.theme
  fds <- fds + scale_color_discrete(name = "Design")
  return(fds)
}


