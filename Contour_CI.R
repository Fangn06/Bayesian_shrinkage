library(ggplot2)
library(tidyverse)
library(palmerpenguins)

# Ridge ------
rho <- 0
s1 <- 0.5
s2 <- 0.5
m1 <- 3
m2 <- 2
V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
Vinv <- solve(V)
m <- c(m1,m2)

# Bayesian estimates ----
Bayes_est <- function(lam){
  
  f <- function(d, lam=1){( (t(d-m)%*% Vinv %*%(d-m)) + lam*sum(d^2) )}
  
  nlm1 <- nlm(function(d){f(d, lam=lam)}, c(m1,m2))
  nlm2 <- nlm(function(d){f(d, lam=lam)}, c(0,0))
  
  if(nlm1$min < nlm2$min){nlm.best <- nlm1} else { nlm.best <- nlm2 }
  
  c(nlm.best$estimate, nlm.best$min, lam=lam)
}

lam_cand <- c(0.5,10,16,30)

estd0_ridge <- Bayes_est(0)
estd1_ridge <- Bayes_est(lam_cand[1])
estd2_ridge <- Bayes_est(lam_cand[2])
estd3_ridge <- Bayes_est(lam_cand[3])
estd4_ridge <- Bayes_est(lam_cand[4])



# Plots ------
plot.data <- function(phi=0,lam)
{
  rho <- 0
  s1 <- 0.5
  s2 <- 0.5
  m1 <- 3
  m2 <- 2
  V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
  Vinv <- solve(V)
  m <- c(m1,m2)
  alpha <- 0.05
  n <- 1000
  
  # functions
  f_phi <- function(d, phi=1){log(1 + t(m-d) %*% Vinv %*% (m-d)) + phi*sum(d^2)}
  f_lam <- function(d, lam=1){t(m-d) %*% Vinv %*% (m-d) + lam*sum(d^2)}
  
  
  
  D1_pos <- vector()
  D2_pos <- vector()
  
  for (i in 1:n) {
    D1_pos[i] =rnorm(1, m1, 0.5);
    D2_pos[i] =rnorm(1, m2, 0.5);
  }
  
  
  # grid simulation
  d1=seq(from=-8,to=8,by=0.1)
  d2=seq(from=-8,to=8,by=0.1)
  f_new=vector()
  D1=vector()
  D2=vector()
  for (i in 1:length(d1)) {
    for (j in 1:length(d2)){
      D1=c(D1,d1[i]);
      D2=c(D2,d2[j]);
    }
  }
  
  
  lossvalue_phi <- vector()
  lossvalue_lam <- vector()
  
  for (i in 1:n) {
    lossvalue_phi[i] <- f_phi(c(D1_pos[i], D2_pos[i]),phi)
    lossvalue_lam[i] <- f_lam(c(D1_pos[i], D2_pos[i]),lam)
  }
  
  gamma_phi <- quantile(lossvalue_phi, probs = 1-alpha);
  gamma_lam <- quantile(lossvalue_lam, probs = 1-alpha);
  
  
  f_new_phi <- vector()
  f_new_lam <- vector()
  for (i in 1:length(D1)) {
    f_new_phi[i]=f_phi(c(D1[i],D2[i]),phi)
    f_new_lam[i]=f_lam(c(D1[i],D2[i]),lam)
  }
  
  F_new_phi=cbind(D1,D2,f_new_phi)
  F_new_phi <- data.frame(F_new_phi)
  colnames(F_new_phi) <- c("d1","d2","value")
  break_lines_phi <- c(min(F_new_phi$value),gamma_phi,max(F_new_phi$value))
  
  F_new_lam=cbind(D1,D2,f_new_lam)
  F_new_lam <- data.frame(F_new_lam)
  colnames(F_new_lam) <- c("d1","d2","value")
  break_lines_lam <- c(min(F_new_lam$value),gamma_lam,max(F_new_lam$value))
  
  return(list(F_new_lam=F_new_lam,break_lines_lam=break_lines_lam,
              F_new_phi=F_new_phi,break_lines_phi=break_lines_phi))
  
  
  
}



cand0_ridge <- plot.data(lam = 0)
cand1_ridge <- plot.data(lam = lam_cand[1])
cand2_ridge <- plot.data(lam = lam_cand[2])
cand3_ridge <- plot.data(lam = lam_cand[3])
cand4_ridge <- plot.data(lam = lam_cand[4])
cand8_ridge <- plot.data(lam = 1000)

shrinkage1 <- c("0.5","10","16","30","Inf")
shrinkage <- "0"


figure.ridge <-  ggplot() + 
  geom_point(aes(x=estd0_ridge[1],y=estd0_ridge[2]), shape=16)+
  geom_contour(data=cand0_ridge$F_new_phi,mapping=aes(x=d1,y=d2,z=value,colour= shrinkage),
               breaks =cand0_ridge$break_lines_phi) + 
  
  geom_point(aes(x=estd1_ridge[1],y=estd1_ridge[2],colour=shrinkage1[2]), shape=16)+
  geom_contour(data=cand1_ridge$F_new_lam, mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[2]),
               breaks =cand1_ridge$break_lines_lam) + 
  
  geom_point(aes(x=estd2_ridge[1],y=estd2_ridge[2],colour=shrinkage1[3]),shape=16)+
  geom_contour(data=cand2_ridge$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[3]),
               breaks =cand2_ridge$break_lines_lam)+
  
  geom_point(aes(x=estd3_ridge[1],y=estd3_ridge[2],colour=shrinkage1[4]),shape=16)+
  geom_contour(data=cand3_ridge$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[4]),
               breaks =cand3_ridge$break_lines_lam)+
  
  geom_point(aes(x=estd4_ridge[1],y=estd4_ridge[2],colour=shrinkage1[5]),shape=16)+
  geom_contour(data=cand4_ridge$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[5]),
               breaks =cand4_ridge$break_lines_lam) +
  
  geom_contour(data=cand8_ridge$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour= shrinkage1[6]),
               breaks =cand8_ridge$break_lines_lam,linetype = "dashed") +
  
  scale_color_manual(values=c("black","orange","green","yellow","purple","gray")) +
  
  labs(x=quote(d[1]),y=quote(d[2]))+ggtitle("Ridge")+theme_bw()+ 
  xlim(c(-6,6)) + ylim(c(-6,6))+ 
  theme(plot.title = element_text(hjust = 0.5))


# Lasso ----
m1 <- 3
m2 <- 2
V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
Vinv <- solve(V)
m <- c(m1,m2)
alpha <- 0.05
n <- 1000


# Bayesian estimates ----
Bayes_est <- function(lam){
  
  f <- function(d, lam=1){( (t(d-m)%*% Vinv %*%(d-m)) + lam*sum(abs(d))  )}
  
  nlm1 <- nlm(function(d){f(d, lam=lam)}, c(m1,m2))
  nlm2 <- nlm(function(d){f(d, lam=lam)}, c(0,0))
  
  if(nlm1$min < nlm2$min){nlm.best <- nlm1} else { nlm.best <- nlm2 }
  
  c(nlm.best$estimate, nlm.best$min, lam=lam)
}

lam_cand <- c(0.5,10,16,30)

estd0 <- Bayes_est(0)
estd1 <- Bayes_est(lam_cand[1])
estd2 <- Bayes_est(lam_cand[2])
estd3 <- Bayes_est(lam_cand[3])
estd4 <- Bayes_est(lam_cand[4])


# Plots ------
plot.data <- function(phi=0,lam)
{
  rho <- 0
  s1 <- 0.5
  s2 <- 0.5
  m1 <- 3
  m2 <- 2
  V <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2,2)
  Vinv <- solve(V)
  m <- c(m1,m2)
  alpha <- 0.05
  n <- 1000
  
  # functions
  f_phi <- function(d, phi=1){log(1 + t(m-d) %*% Vinv %*% (m-d)) + phi*sum(abs(d))}
  f_lam <- function(d, lambda){
    ( (t(d-m)%*% Vinv %*%(d-m)) + lam*sum(abs(d)) )
  }  
  
  
  D1_pos <- vector()
  D2_pos <- vector()
  
  for (i in 1:n) {
    D1_pos[i] =rnorm(1, m1, 0.5);
    D2_pos[i] =rnorm(1, m2, 0.5);
  }
  
  
  # grid simulation
  d1=seq(from=-8,to=8,by=0.1)
  d2=seq(from=-8,to=8,by=0.1)
  f_new=vector()
  D1=vector()
  D2=vector()
  for (i in 1:length(d1)) {
    for (j in 1:length(d2)){
      D1=c(D1,d1[i]);
      D2=c(D2,d2[j]);
    }
  }
  
  
  lossvalue_phi <- vector()
  lossvalue_lam <- vector()
  
  for (i in 1:n) {
    lossvalue_phi[i] <- f_phi(c(D1_pos[i], D2_pos[i]),phi)
    lossvalue_lam[i] <- f_lam(c(D1_pos[i], D2_pos[i]),lam)
  }
  
  gamma_phi <- quantile(lossvalue_phi, probs = 1-alpha);
  gamma_lam <- quantile(lossvalue_lam, probs = 1-alpha);
  
  
  f_new_phi <- vector()
  f_new_lam <- vector()
  for (i in 1:length(D1)) {
    f_new_phi[i]=f_phi(c(D1[i],D2[i]),phi)
    f_new_lam[i]=f_lam(c(D1[i],D2[i]),lam)
  }
  
  F_new_phi=cbind(D1,D2,f_new_phi)
  F_new_phi <- data.frame(F_new_phi)
  colnames(F_new_phi) <- c("d1","d2","value")
  break_lines_phi <- c(min(F_new_phi$value),gamma_phi,max(F_new_phi$value))
  
  F_new_lam=cbind(D1,D2,f_new_lam)
  F_new_lam <- data.frame(F_new_lam)
  colnames(F_new_lam) <- c("d1","d2","value")
  break_lines_lam <- c(min(F_new_lam$value),gamma_lam,max(F_new_lam$value))
  
  return(list(F_new_lam=F_new_lam,break_lines_lam=break_lines_lam,
              F_new_phi=F_new_phi,break_lines_phi=break_lines_phi))
  
  
  
}


cand0 <- plot.data(lam = 0)
cand1 <- plot.data(lam = lam_cand[1])
cand2 <- plot.data(lam = lam_cand[2])
cand3 <- plot.data(lam = lam_cand[3])
cand4 <- plot.data(lam = lam_cand[4])
cand8 <- plot.data(lam = 1000)

shrinkage1 <- c("0.5","10","16","30","Inf")
shrinkage <- "0"


figure.lasso <-  ggplot() + 
  geom_contour(data=cand0$F_new_phi,mapping=aes(x=d1,y=d2,z=value,colour= shrinkage),
               breaks =cand0$break_lines_phi) + 
  geom_point(aes(x=estd0[1],y=estd0[2]), shape=16)+
  
  geom_point(aes(x=estd1[1],y=estd1[2],colour=shrinkage1[2]), shape=16)+
  geom_contour(data=cand1$F_new_lam, mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[2]),
               breaks =cand1$break_lines_lam) + 
  
  geom_point(aes(x=estd2[1],y=estd2[2],colour=shrinkage1[3]),shape=16)+
  geom_contour(data=cand2$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[3]),
               breaks =cand2$break_lines_lam)+
  
  geom_point(aes(x=estd3[1],y=estd3[2],colour=shrinkage1[4]),shape=16)+
  geom_contour(data=cand3$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[4]),
               breaks =cand3$break_lines_lam)+
  
  geom_point(aes(x=estd4[1],y=estd4[2],colour=shrinkage1[5]),shape=16)+
  geom_contour(data=cand4$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour=shrinkage1[5]),
               breaks =cand4$break_lines_lam) +
  
  geom_contour(data=cand8$F_new_lam,mapping=aes(x=d1,y=d2,z=value,colour= "Inf"),
               breaks =cand8$break_lines_lam,linetype = "dashed")  + 
  
  scale_color_manual(values=c("black","orange","green","yellow","purple","gray")) +
  
  labs(x=quote(d[1]),y=quote(d[2]))+ggtitle("Lasso")+theme_bw()+ 
  xlim(c(-6,6)) + ylim(c(-6,6))+ 
  theme(plot.title = element_text(hjust = 0.5))