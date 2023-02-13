library(lattice)

# Ridge ----
rho <- 0
s1 <- 0.5
s2 <- 0.5
m1 <- 3
m2 <- 2
V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
Vinv <- solve(V)
m <- c(m1,m2)
f <- function(d, phi=1){log(1 + t(m-d) %*% Vinv %*% (m-d)) + phi*sum(d^2)}

get.ridge <- function(lambda){
  as.vector(solve( diag(rep(1,2)) + lambda*V ) %*% c(m1,m2) )}

quick.opt.ridge <- function(phix){
  gridlam <- seq(0,100, 0.1)
  fx <- function(lam){f(get.ridge(lam), phi=phix)}
  gridfx <- sapply(gridlam, fx)
  optim1 <- optimize(fx, sort( gridlam[order(gridfx)[1:2]] ) )
  c(get.ridge(optim1$minimum), lam=optim1$minimum)
}

get.lambda.ridge <- function(dvec){
  #	c(m1,m2) - dvec = lambda*V %*% dvec (come from ?)
  (c(m1,m2)-dvec)/((V %*% dvec)[,1])
}


phivals.ridge <- seq(0,2, l=1001)
alld.ridge <- sapply(phivals.ridge, function(phi){
  
  nlm1 <- nlm(function(d){f(d, phi=phi)}, c(m1,m2))
  nlm2 <- nlm(function(d){f(d, phi=phi)}, c(0,0))
  
  if(nlm1$min < nlm2$min){nlm.best <- nlm1} else { nlm.best <- nlm2 }
  qq <- quick.opt.ridge(phi)
  
  if(f(qq[1:2], phi) < nlm.best$min){
    nlm.best$estimate <- qq[1:2]
    nlm.best$minimum <- f(qq[1:2], phi)
  }
  
  c(nlm.best$estimate, nlm.best$min, phi=phi)
})


all.lam.lambda.ridge <- sapply( phivals.ridge, function(phi){get.lambda.ridge(quick.opt.ridge(phi)[1:2])})
all.lam.ridge <- colMeans(all.lam.lambda.ridge)

par(mar=c(4,5,0,0)+1.5)
# phi version plot ----
mycut <- 0.303
plot(phivals.ridge, alld.ridge[1,], ylim=c(0,3), xlim=c(0,1),
     xlab=expression(phi), ylab="Estimate, d", type="n", main = "Ridge",axes=FALSE)
abline(v=mycut, lty=2, col="grey50")
abline(h=quick.opt.ridge(0.302)[1:2], lty=3, col=1:2)
abline(h=quick.opt.ridge(0.303)[1:2], lty=3, col=1:2)
lines(phivals.ridge[phivals.ridge<mycut], alld.ridge[1,phivals.ridge<mycut], type="l", col=1)
lines(phivals.ridge[phivals.ridge>mycut], alld.ridge[1,phivals.ridge>mycut], type="l", col=1)
lines(phivals.ridge[phivals.ridge<mycut], alld.ridge[2,phivals.ridge<mycut], type="l", col=2)
lines(phivals.ridge[phivals.ridge>mycut], alld.ridge[2,phivals.ridge>mycut], type="l", col=2)
axis(side=1, at=mycut)
axis(side=1, at=c(0,0.5, 1))
axis(side=2, at=round(quick.opt.ridge(0.302)[1:2],2), las=1)
axis(side=2, at=round(quick.opt.ridge(0.303)[1:2],2), las=1)
box()

# lambda version plot ----
f <- function(d, phi=1){t(m-d) %*% Vinv %*% (m-d) + phi*sum(d^2)}
lamvals.ridge <- seq(0, 30, by=0.01)
lam.ridge <- sapply(lamvals.ridge, function(phi){
  
  nlm1 <- nlm(function(d){f(d, phi=phi)}, c(m1,m2))
  nlm2 <- nlm(function(d){f(d, phi=phi)}, c(0,0))
  
  if(nlm1$min < nlm2$min){nlm.best <- nlm1} else { nlm.best <- nlm2 }
  c(nlm.best$estimate, nlm.best$min, phi=phi)
})

plot(all.lam.ridge, alld.ridge[1,], ylim=c(0,3), xlab=expression(lambda),
     ylab="Estimate, d", type="n", xlim=c(0,20), axes=FALSE, main = "Ridge")
abline(h=c(2.67,1.78), lty=3, col=1:2)
abline(h=c(1.21,0.81), lty=3, col=1:2)
abline(v=0.49, lty=2, col="gray50")
abline(v=5.9, lty=2, col="gray50")
lines(lam.ridge[4,lam.ridge[4,]<0.49], lam.ridge[1,lam.ridge[4,]<0.49], type="l", col=1)
lines(lam.ridge[4,lam.ridge[4,]>5.9], lam.ridge[1,lam.ridge[4,]>5.9], type="l", col=1)
lines(lam.ridge[4,lam.ridge[4,]<0.49], lam.ridge[2,lam.ridge[4,]<0.49], type="l", col=2)
lines(lam.ridge[4,lam.ridge[4,]>5.9], lam.ridge[2,lam.ridge[4,]>5.9], type="l", col=2)

lines(lam.ridge[4,lam.ridge[4,]>=0.49 & lam.ridge[4,]<=5.9], lam.ridge[1,lam.ridge[4,]>=0.49 & lam.ridge[4,]<=5.9], lty=2, col=1)

lines(lam.ridge[4,lam.ridge[4,]>=0.49 & lam.ridge[4,]<=5.9], lam.ridge[2,lam.ridge[4,]>=0.49 & lam.ridge[4,]<=5.9], lty=2, col=2)

axis(side=1, at=0.49, labels=paste("        ",0.49))
axis(side=1, at=5.9)
axis(side=1, at=c(0,10,15,20))
axis(side=2, at=c(2.67,1.78,1.21,0.81), las=1)
box()
legend("topright", lty=1, col=1:2, c(expression(d[1]), expression(d[2])),
       bty="n", horiz=TRUE)





# Lasso ----
rho <- 0
s1 <- 0.5
s2 <- 0.5
m1 <- 3
m2 <- 2
V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
Vinv <- solve(V)
m <- c(m1,m2)

f <- function(d, phi=1){
  ( log(1 + (t(d-m)%*% Vinv %*%(d-m)) ) + phi*sum(abs(d)) )
}

get.lasso <- function(lam){
  solve( 2*Vinv, -lam*c(1,1) + 2*Vinv%*%c(m1,m2))[,1]
}

quick.opt.lasso <- function(phix){
  gridlam <- seq(0,100, 0.1)
  fx <- function(lam){f(get.lasso(lam), phi=phix)}
  gridfx <- sapply(gridlam, fx)
  optim1 <- optimize(fx, sort( gridlam[order(gridfx)[1:2]] ) )
  c(get.lasso(optim1$minimum), lam=optim1$minimum)
}

get.lambda.lasso <- function(dvec){
  # lambda*sign(d)-2*Vinv%*%c(m1,m2)+2*Vinv%*%dvec=0 
  (2*Vinv%*%c(m1,m2)-2*Vinv%*%dvec)
}




phivals.lasso <- seq(0,2, l=1001)
alld.lasso <- sapply(phivals.lasso, function(phi){
  
  nlm1 <- nlm(function(d){f(d, phi=phi)}, c(m1,m2))
  nlm2 <- nlm(function(d){f(d, phi=phi)}, c(0,0))
  
  if(nlm1$min < nlm2$min){nlm.best <- nlm1} else { nlm.best <- nlm2 }
  qq <- quick.opt.lasso(phi)
  
  if(f(qq[1:2], phi) < nlm.best$min){
    nlm.best$estimate <- qq[1:2]
    nlm.best$minimum <- f(qq[1:2], phi)
  }
  
  c(nlm.best$estimate, nlm.best$min, phi=phi)
})

all.lam.lambda.lasso <- sapply( phivals.lasso, function(phi){get.lambda.lasso(quick.opt.lasso(phi)[1:2])})
all.lam.lasso <- colMeans(all.lam.lambda.lasso)



par(mar=c(4,5,0,0)+1.5)
# phi version plot ----
mycut <- 0.812
plot(phivals.lasso, alld.lasso[1,], ylim=c(0,3), xlim=c(0,1.6), 
     xlab=expression(phi), ylab="Estimate, d", type="n", axes=FALSE, main = "Lasso")
abline(h=0, lty=3, col=1:2)
abline(v=mycut, lty=2, col="grey50")

abline(h=quick.opt.lasso(0.811)[1:2], lty=3, col=1:2)
abline(h=quick.opt.lasso(0.812)[1:2], lty=3, col=1:2)

lines(phivals.lasso[phivals.lasso<mycut], alld.lasso[1,phivals.lasso<mycut], type="l", col=1)
lines(phivals.lasso[phivals.lasso>mycut], alld.lasso[1,phivals.lasso>mycut], type="l", col=1)

lines(phivals.lasso[phivals.lasso<mycut], alld.lasso[2,phivals.lasso<mycut], type="l", col=2)
# The same as d1 and lines are overlapped
lines(phivals.lasso[phivals.lasso>mycut], alld.lasso[2,phivals.lasso>mycut], type="l", col=2)

axis(side=1, at=mycut)
axis(side=1, at=c(0,0.5,1.5))
axis(side=2, at=0, las=1)
axis(side=2, at=round(quick.opt.lasso(0.811)[1:2],2), las=1)
# The same as 0.811 (same value)
axis(side=2, at=round(quick.opt.lasso(0.812)[1:2],2), las=1)
box()

# lambda version plot ----
plot(all.lam.lasso, alld.lasso[1,], ylim=c(0,3), xlim=c(0,30), xlab=expression(lambda), 
     ylab="Estimate, d", type="n", axes=FALSE, main = "Lasso")

abline(h=c(2.89,1.89), lty=3, col=1:2)
abline( v=24, lty=2,  col="gray50")
abline( v=0.89, lty=2, col="gray50")
abline(h=0, lty=3, col=1)

lines(c(0,0.89), c(3,2.89), type="l", col=1)
lines(c(0.89,24), c(2.89,0), type="l", col=1,lty=2)
lines(c(24,30), c(0,0), type="l", col=2)

lines(c(0,0.89), c(2,1.89), type="l", col=2)
lines(c(0.89,16), c(1.89,0), type="l", col=2,lty=2)
lines(c(16,24), c(0,0), type="l", col=2,lty=2)

axis(side=1, at=0.89, labels=paste("        ",0.89))
axis(side=1, at=c(0,16,24,30), las=1)
axis(side=2, at=c(1.89,2.89), las=1)
axis(side=2, at=0, las=1)
box()
legend("topright", lty=1, col=1:2, c(expression(d[1]), expression(d[2])),
       bty="n", horiz=TRUE)
