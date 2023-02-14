library("lars")
data(diabetes)

# Trajectory ----
diabetes$y2 <- diabetes$y - mean(diabetes$y)

myx <- diabetes$x[,1:4]

myx2 <- myx*sqrt(441)

lm1 <- lm(y2~ myx2 - 1, data=diabetes)

lars1 <- with(diabetes, lars( myx2[,1:4], y2 ), intercept=FALSE )
plot(lars1)

# ingredients for classic version
mu <- coef(lm1)
Vinv <- solve( vcov(lm1) )

# ingredients for actually-Bayes version
n <- nrow(diabetes) #442
g <- 0.2*n
s2 <- 60^2 #summary(lm1)$sigma^2
X <- myx2[,1:4]

post.mean.beta <- (g/(g+1)) * coef(lm1)
post.var.beta  <- (g/(g+1)) * s2 * solve(t(X)%*%X)

mu <- post.mean.beta
Vinv <- solve(post.var.beta)

g <- function(d, lam){
  (t(d-mu)%*% Vinv %*%(d-mu))[1,1] + lam*sum(abs(d)) 
}

llvals <- seq(0,1,l=11)
sapply( llvals, function(lam){
  nlm1 <- nlm(function(d){g(c(0,0,d,0),lam)}, 0)
  nlm2 <- nlm(function(d){g(c(0,0,d),  lam)}, rep(0,2))
  nlm3 <- nlm(function(d){g(c(0,d),    lam)}, rep(0,3))
  nlm4 <- nlm(function(d){g(d,lam)}, coef(lm1))
  print(order(c(nlm1$minimum,nlm2$minimum,nlm3$minimum,nlm4$minimum)))
  list(nlm1, nlm2, nlm3, nlm4)[order(c(nlm1$minimum,nlm2$minimum,nlm3$minimum,nlm4$minimum))[1]]
  invisible()
})

# optimum solves
# 
# 2 Vinv %*% (d-mu) + lam*sign(d) == 0, for non-zero d_i elements
# 
# e.g. for 3, 4 this is
# 
# 2 * Vinv[3,] %*% (c(0,0,d3,d4)-mu) + lam*sign(d) == 0
# 2 * Vinv[4,] %*% (c(0,0,d3,d4)-mu) + lam*sign(d) == 0
# 
# which is
# 
# 2 * Vinv[3,]%*%-mu + Vinv[3,1:2]%*%c(0,0) + 2* Vinv[3,3:4]%*%d[3:4] + lam*sign(d)[3] == 0
# 2 * Vinv[4,]%*%-mu + Vinv[4,1:2]%*%c(0,0) + 2* Vinv[4,3:4]%*%d[3:4] + lam*sign(d)[4] == 0

sol1 <- function(lam){
  d <- solve( 2*Vinv[3,3], -lam + 2*Vinv[3,]%*%mu)[,1]
  c(0,0,d,0)
}

sol2 <- function(lam){
  d <- solve( 2*Vinv[3:4,3:4], -lam*c(1,1) + 2*Vinv[3:4,]%*%mu)[,1]
  c(0,0,d)
}

sol3 <- function(lam){
  d <- solve( 2*Vinv[2:4,2:4], -lam*c(-1,1,1) + 2*Vinv[2:4,]%*%mu)[,1]
  c(0,d)
}

sol4 <- function(lam){
  solve( 2*Vinv[1:4,1:4], -lam*c(1,-1,1,1) + 2*Vinv[1:4,]%*%mu)[,1]
}

crit12 <- uniroot( function(lam){ g(sol1(lam),lam) - g(sol2(lam),lam) } , c(5,10))$root
sol1(crit12)
sol2(crit12)
crit23 <- uniroot( function(lam){ g(sol2(lam),lam) - g(sol3(lam),lam) } , c(0,5))$root
sol2(crit23)
sol3(crit23)
crit34 <- uniroot( function(lam){ g(sol3(lam),lam) - g(sol4(lam),lam) } , c(0,1))$root
sol3(crit34)
sol4(crit34)

betamat.bayes <- rbind(
  rep(0,4),
  sol1(crit12),
  sol2(crit23),
  sol3(crit34),
  mu)

f <- function(d, lam){
  (t(d-mu)%*% Vinv %*%(d-mu)) + lam*sum(abs(d)) 
}

my.interp <- function(prop, betamat){
  n <- nrow(betamat)
  prop.vals <- (0:(n-1))/(n-1)
  a1 <- approxfun(prop.vals, betamat[,1])
  a2 <- approxfun(prop.vals, betamat[,2])
  a3 <- approxfun(prop.vals, betamat[,3])
  a4 <- approxfun(prop.vals, betamat[,4])
  c(a1(prop), a2(prop), a3(prop), a4(prop) )
}


get.lmin <- function(lam, betamat, ngrid=101){
  ppvals <- seq(0,1,l=ngrid)
  interp.beta <- sapply(ppvals, function(p){my.interp(p, betamat)})
  interp.L <- apply(interp.beta, 2, function(d){f(d, lam=lam)})
  pp.bounds <- ppvals[order(interp.L)[1:2]]
  opt1 <- optimize(function(pp){f(my.interp(pp, betamat), lam)}, pp.bounds[1:2])
  c(my.interp(opt1$minimum, betamat), opt1$minimum, opt1$objective)
}


lamvals <- seq(0,12, l=1001)
nlmsols <- sapply(lamvals, function(lam){ 
  nlm1 <- nlm( function(d){f(d, lam)}, mu )
  nlm0 <- nlm( function(d){f(d, lam)}, rep(0, 4) )
  n.list <- list(nlm1, nlm0)[order(c(nlm1$min, nlm0$min))]
  c(n.list[[1]]$est, n.list[[1]]$min)
})

mysols <- sapply(phivals, function(lam){get.lmin(lam, betamat.bayes)})
allsols <- t(mysols)[,1:4]
allsols[ t(nlmsols)[,5] < t(mysols)[,6],] <- t(nlmsols)[t(nlmsols)[,5] < t(mysols)[,6],1:4]


lam <- c(0,lamvals[c(which(allsols[,1]==0)[1],which(allsols[,2]==0)[1],
                     which(allsols[,4]==0)[1],which(allsols[,3]==0)[1])])

plot(x=lam, y=rev(lars1$beta[,3]), xlab=expression(lambda), ylab="Estimates", type="n", main = "Diabetes", ylim=range(lars1$beta),xlim=c(0,12),axes=FALSE)
lines(lam, rev(lars1$beta[,1]), col = "black", lwd = 1)
lines(lam, rev(lars1$beta[,2]), col = "red", lwd = 1)
lines(lam, rev(lars1$beta[,3]), col = "blue", lwd = 1)
lines(lam, rev(lars1$beta[,4]), col = "green", lwd = 1)
lines(lam, rev(betamat.bayes[,1]), col = "black", lwd = 1, lty=2)
lines(lam, rev(betamat.bayes[,2]), col = "red", lwd = 1, lty=2)
lines(lam, rev(betamat.bayes[,3]), col = "blue", lwd = 1, lty=2)
lines(lam, rev(betamat.bayes[,4]), col = "green", lwd = 1, lty=2)
abline(v=lam[-5], lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
axis(side=1, at=3)
axis(side=2, at=c(-10,0,10,20,30,40))
axis(side=1, at=round(lam,2))
axis(side=1, at=3)
box()
legend("topright", lty=1, col=c("black","red","blue","green"),
       c("age", "sex","BMI","map"),
       bty="n", horiz=FALSE)


# prior - volume ----
# interval volume -- regret method
myx <- diabetes$x[,1:4]
lm1 <- lm(I(y-mean(y))~myx-1, data=diabetes)

plot.data <- function(vers,inf, N.mc)
{
  d_grid=seq(from=-500,to=500,by=0.01)
  n=1000
  alpha=0.05

  # ingredients for actually-Bayes version
  n <- nrow(diabetes)
  g <- inf*n
  s2 <- 60^2 #summary(lm1)$sigma^2
  X <- myx

  post.mean.beta <- (g/(g+1)) * coef(lm1)
  post.var.beta  <- (g/(g+1)) * s2 * solve(t(X)%*%X)

  mu <- post.mean.beta
  V <- post.var.beta
  Vinv <- solve(post.var.beta)
  D_pos <- mvrnorm(n,mu=mu,Sigma=V)
  D1_pos <- D_pos[,1]
  D2_pos <- D_pos[,2]
  D3_pos <- D_pos[,3]
  D4_pos <- D_pos[,4]



  if(vers=="phi") {
    f <- function(d, phi){
      ( log(1 + (t(d-mu)%*% Vinv %*%(d-mu)) ) + phi*sum(abs(d)) )
    }
    paravalues <- seq(0,1,0.01)
  }

  if(vers=="lambda") {
    f <- function(d, lambda){
      ( (t(d-mu)%*% Vinv %*%(d-mu)) + lambda*sum(abs(d)) )
    }
    paravalues <- seq(0,1,0.01)
  }


  area.allpara <- vector()

  for (k in 1:length(paravalues)) {

    set.seed(2022)
    lossvalue <- vector()

    for (i in 1:n) {
      lossvalue[i] <- f(c(D1_pos[i], D2_pos[i], D3_pos[i],
                          D4_pos[i]),paravalues[k])
    }

    gamma <- quantile(lossvalue, probs = 1-alpha);

    mc_loss <- vector()
    mc_d1 <- sample(d_grid,N.mc,replace = TRUE)
    mc_d2 <- sample(d_grid,N.mc,replace = TRUE)
    mc_d3 <- sample(d_grid,N.mc,replace = TRUE)
    mc_d4 <- sample(d_grid,N.mc,replace = TRUE)


    for(i in 1:length(mc_d1))
    {
      mc_loss[i] <- f(c(mc_d1[i], mc_d2[i], mc_d3[i],
                        mc_d4[i]),paravalues[k])
    }

    index.within <- which(mc_loss <= gamma)
    area.allpara[k]  <- length(index.within)
  }

  datasave <- cbind(paravalues,area.allpara)
  return(datasave)

}

vol_g1 <- plot.data(vers="lambda",inf = 0.2, N.mc=10^12)
vol_g2 <- plot.data(vers="lambda",inf = 0.01, N.mc=10^12)

plot(vol_g1[,1],log10(vol_g1[,2]),type = "l",ylab = expression(log[10](volume)), xlab=expression(lambda), lwd=1,ylim = c(9,12))lines(vol_g2[,1],vol_g2[,2],type = "l", col = "blue", lwd=1)
lines(vol_g2[,1],log10(vol_g2[,2]),type = "l", col = "blue", lwd=1)
