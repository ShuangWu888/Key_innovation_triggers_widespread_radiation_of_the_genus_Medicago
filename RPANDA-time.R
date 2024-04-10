#library(phytools)
library(ape)
library(nlme)
library(permute)
library(lattice)
library(vegan)
#library(permute)
library(picante)
#library(lattice)
#library(nlme)
library(RPANDA)
sativa <- read.tree("r8s.84species-lowcop.newtime-smooth_100.out.timetree.rename.prune.reverse_order.Medicago.tre")
#data(sativa)
tot_time<-max(node.age(sativa)$ages)

# Fit the pure birth model (no extinction) with a constant speciation rate
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.36)
mu_par<-c()
result_bcst_d0 <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,cst.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
result_bcst_d0$model <- "pure birth with constant speciation rate"

# Fit a birth-death model with constant speciation rate and constant extinction
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.36)
mu_par<-c(0.02)
result_bcst_dcst <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,cst.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
result_bcst_dcst$model <- "birth-death model with constant speciation rate and constant extinction"

# Fit a birth-death model with constant speciation rate and exponential variation of the extinction rate with time
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){y[1] * exp(y[2] * t)}
lamb_par<-c(0.36)
mu_par<-c(-1.192e-8,5.377e-3)
result_bcst_dexp <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,cst.lamb=TRUE,expo.mu=TRUE,dt=1e-3)
result_bcst_dexp$model <- "birth-death model with constant speciation rate and exponential variation of the extinction rate with time"

# Fit a birth-death model with constant speciation rate and linear variation of the extinction rate with time
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){abs(y[1] + y[2] * t)}
lamb_par<-c(0.36)
mu_par<-c(0.025,-0.0016)
result_bcst_dlin <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,cst.lamb=TRUE,dt=1e-3)
result_bcst_dlin$model <- "birth-death model with constant speciation rate and linear variation of the extinction rate with time"

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.046,0.0072)
mu_par<-c()
result_bexp_d0 <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,expo.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
result_bexp_d0$model <- "birth model (no extinction) with exponential variation of the speciation rate with time"

# Fit a birth-death model with exponential variation of the speciation rate with time and constant extinction
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.046,0.0072)
mu_par<-c(0.02)
result_bexp_dcst <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
result_bexp_dcst$model <- "birth-death model with exponential variation of the speciation rate with time and constant extinction"

# Fit a birth-death model with exponential variation of the speciation rate with time and exponential variation of the extinction rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){y[1] * exp(y[2] * t)}
lamb_par<-c(0.046,0.0072)
mu_par<-c(1.182e-8,8.707e-3)
result_bexp_dexp <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,expo.lamb=TRUE,expo.mu=TRUE,dt=1e-3)
result_bexp_dexp$model <- "birth-death model with exponential variation of the speciation rate with time and exponential variation of the extinction rate with time"

# Fit a birth-death model with exponential variation of the speciation rate with time and linear variation of the extinction rate with time
#f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
#f.mu<-function(t,y){pmax(0,y[1] + y[2] * t)}
#lamb_par<-c(0.046,0.0072)
#mu_par<-c(0.025,-0.0016)
#result_bexp_dlin <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,expo.lamb=TRUE,dt=1e-3)
#result_bexp_dlin$model <- "birth-death model with exponential variation of the speciation rate with time and linear variation of the extinction rate with time"

# Fit the pure birth model (no extinction) with linear variation of the speciation rate with time 
f.lamb <-function(t,y){abs(y[1] + y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.0256,0.0017)
mu_par<-c()
result_blin_d0 <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,fix.mu=TRUE,dt=1e-3)
result_blin_d0$model <- "pure birth model (no extinction) with linear variation of the speciation rate with time"

# Fit a birth-death model with linear variation of the speciation rate with time and constant extinction
f.lamb <-function(t,y){abs(y[1] + y[2] * t)}
f.mu <-function(t,y){y[1]}
lamb_par<-c(0.0256,0.0017)
mu_par<-c(0.02)
result_blin_dcst <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,cst.mu=TRUE,dt=1e-3)
result_blin_dcst$model <- "birth-death model with linear variation of the speciation rate with time and constant extinction"

# Fit a birth-death model with linear variation of the speciation rate with time and exponential variation of the extinction rate with time
#f.lamb <-function(t,y){pmax(0,y[1] + y[2] * t)}
#f.mu<-function(t,y){y[1] * exp(y[2] * t)}
#lamb_par<-c(0.0256,-0.0017)
#mu_par<-c(1.182e-8,8.707e-3)
#result_blin_dexp <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,expo.mu=TRUE,dt=1e-3)
#result_blin_dexp$model <- "birth-death model with linear variation of the speciation rate with time and exponential variation of the extinction rate with time"

# Fit a birth-death model with linear variation of the speciation rate with time and linear variation of the extinction rate with time
#f.lamb <-function(t,y){pmax(0,y[1] + y[2] * t)}
#f.mu<-function(t,y){pmax(0,y[1] + y[2] * t)}
#lamb_par<-c(0.0256,-0.0017)
#mu_par<-c(0.025,-0.0016)
#result_blin_dlin <- fit_bd(sativa,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=55/90,dt=1e-3)
#result_blin_dlin$model <- "birth-death model with linear variation of the speciation rate with time and linear variation of the extinction rate with time"

# Find the best model
#index <- which.min(c(result_bcst_d0$aicc, result_bcst_dcst$aicc, result_bcst_dexp$aicc, result_bcst_dlin$aicc, result_bexp_d0$aicc, result_bexp_dcst$aicc, result_bexp_dexp$aicc, result_bexp_dlin$aicc, result_blin_d0$aicc, result_blin_dcst$aicc, result_blin_dexp$aicc, result_blin_dlin$aicc))
#rbind(result_bcst_d0, result_bcst_dcst, result_bcst_dexp, result_bcst_dlin, result_bexp_d0, result_bexp_dcst, result_bexp_dexp, result_bexp_dlin, result_blin_d0, result_blin_dcst, result_blin_dexp, result_blin_dlin)[index,]
index <- which.min(c(result_bcst_d0$aicc))
rbind(result_bcst_d0)[index,]
index <- which.min(c(result_bcst_dcst$aicc))
rbind(result_bcst_dcst)[index,]
index <- which.min(c(result_bcst_dexp$aicc))
rbind(result_bcst_dexp)[index,]
index <- which.min(c(result_bcst_dlin$aicc))
rbind(result_bcst_dlin)[index,]
index <- which.min(c(result_bexp_d0$aicc))
rbind(result_bexp_d0)[index,]
index <- which.min(c(result_bexp_dcst$aicc))
rbind(result_bexp_dcst)[index,]
index <- which.min(c(result_bexp_dexp$aicc))
rbind(result_bexp_dexp)[index,]
#index <- which.min(c(result_bexp_dlin$aicc))
#rbind(result_bexp_dlin)[index,]
index <- which.min(c(result_blin_d0$aicc))
rbind(result_blin_d0)[index,]
index <- which.min(c(result_blin_dcst$aicc))
rbind(result_blin_dcst)[index,]
#index <- which.min(c(result_blin_dexp$aicc))
#rbind(result_blin_dexp)[index,]
#index <- which.min(c(result_blin_dlin$aicc))
#rbind(result_blin_dlin)[index,]
