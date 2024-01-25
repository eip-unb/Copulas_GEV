library(extRemes)

#Julho
X = c(13.64, 39.32, 10.66, 224.07, 40.90, 22.22, 14.44, 23.59, 47.02, 37.01, 432.11,
      10.63, 28.51, 11.77, 25.35, 25.80, 39.73, 9.21, 22.36, 11.63, 33.35, 18.00, 18.62,
      17.71, 100.10, 23.32, 11.63, 10.20, 12.04, 11.63, 50.57, 11.63, 33.72, 14.69, 12.30,
      32.90, 179.75, 37.57, 7.95)
length(X)

#Setembro
Y = c(29.19, 8.49, 7.37, 82.93, 44.18, 13.82, 22.28, 28.06, 6.84, 12.14, 153.78,
      17.04, 13.47, 15.43, 30.36, 6.91, 22.12, 35.45, 44.66, 95.81, 6.18, 10.00, 58.39,
      24.05, 17.03, 38.65, 47.17, 27.99, 11.84, 9.60, 6.72, 13.74, 14.60, 9.65, 10.39, 60.14,
      15.51, 14.69, 16.44)
length(Y)



######################################################
##### Copula de Gumbel-Hougaard

#distribution function
C_GH=function(w,z,p1){
  exp(-(( (-log(w))^(1/p1) + (-log(z))^(1/p1) )^p1))
}

#density function
c_GH=function(u,v,p){
  num1=C_GH(w=u,z=v,p1=p)*(log(u)*log(v))^((1/p)-1)
  den=u*v*((-log(u))^(1/p) + (-log(v))^(1/p))^(2-p)
  num2=(((-log(u))^(1/p) + (-log(v))^(1/p))^(p)) +(1/p) -1
  
  return(num1*num2/den)
  
}


#
######################################################
##### Frank

#density function
c_F=function(u,v,p){
  
  num=-p*exp(-p*(u+v))*(exp(-p)-1)
  den=( exp(-p)-exp(-p*u)-exp(-p*v)+exp(-p*(u+v)) )^2
  
  return(num/den)
  
}


######################################################
##### Clayton

#density function
c_C=function(u,v,p){
  
  num1=(1+p)*(u*v)^(-1-p)
  num2=(-1+(u^(-p))+(v^(-p)))^(-2-1/p)
  
  return(num1*num2)
  
}


fit_copula=function(Fx, Fy, copula="Gumbel-Hougaard", N=10){
  
  if(copula=="Gumbel-Hougaard"){
    p_test=seq(0,1,by=10^(-3))[-1]
    
    vet=NULL
    for(i in 1:length(p_test)){
      vet[i]=sum(log(c_GH(Fx,Fy,p=p_test[i])))
    }
  }
  
  if(copula=="Frank"){
    p_test=seq(-N,N,by=10^(-2))
    p_test=p_test[p_test!=0]
    
    vet=NULL
    for(i in 1:length(p_test)){
      vet[i]=sum(log(c_F(Fx,Fy,p=p_test[i])))
    }
  }
  
  if(copula=="Clayton"){
    p_test=seq(0,N,by=10^(-2))
    p_test=p_test[p_test!=0]
    
    vet=NULL
    for(i in 1:length(p_test)){
      vet[i]=sum(log(c_C(Fx,Fy,p=p_test[i])))
    }
  }
  
  p_est=p_test[which.max(vet)]
  return(list(p_est=p_est, ll_max=vet[which.max(vet)]))
}



#######################################
#### Aplicacao

# GEV CDF function
gev_cdf <- function(x,  mu, sigma, k) {
  if (k == 0) {
    p <- exp(-exp(-(x - mu) / sigma))
  } else {
    p <- exp(-(1 + k * (x - mu) / sigma)^(-1 / k))
  }
  return(p)
}

#### Parameter estimation

(fit_x=fevd(X,type="GEV"))
fit_x$results$par

(fit_y=fevd(Y,type="GEV"))
fit_y$results$par

F_X=gev_cdf(X, mu=fit_x$results$par[1],
            sigma =fit_x$results$par[2],
            k=fit_x$results$par[3])

F_Y=gev_cdf(Y, mu=fit_y$results$par[1],
            sigma =fit_y$results$par[2],
            k=fit_y$results$par[3])


(copula_GH=fit_copula(F_X, F_Y, copula="Gumbel-Hougaard"))
(copula_F=fit_copula(F_X, F_Y, copula="Frank"))
(copula_C=fit_copula(F_X, F_Y, copula="Clayton"))


##################### fit_r function ############################
fit_R=function(fit_x, fit_y,theta,N=10^4, method="Gumbel-Hougaard"){
  
  mx=fit_x[1]
  sx =fit_x[2]
  gx=fit_x[3]
  
  my=fit_y[1]
  sy =fit_y[2]
  gy=fit_y[3]
  theta_est=theta
  
  v1=runif(N,0,1)
  F_Ginv=pevd(qevd(v1, loc = my, scale = sy, shape = gy), 
              loc = mx, scale = sx, shape = gx)
  u1=NULL
  for (i in 1:N) {
    u1[i]=runif(1, F_Ginv[i],1)
  }
  
  if(method=="Gumbel-Hougaard"){
    #theta_est=fit_copula(F_X, F_Y, copula="Gumbel-Hougaard")$p_est
    out=mean(c_GH(u1,v1, theta_est)*(1-0)*(1-F_Ginv), na.rm = TRUE)
  }
  if(method=="Frank"){
    #theta_est=fit_copula(F_X, F_Y, copula="Frank")$p_est
    out=mean(c_F(u1,v1, theta_est)*(1-0)*(1-F_Ginv), na.rm = TRUE)
  }
  if(method=="Clayton"){
    #theta_est=fit_copula(F_X, F_Y, copula="Clayton")$p_est
    out=mean(c_C(u1,v1, theta_est)*(1-0)*(1-F_Ginv), na.rm = TRUE)
  }
  
  return(out)
}

fit_R(fit_x$results$par, fit_y$results$par,
               theta=copula_GH$p_est, method="Gumbel-Hougaard")
fit_R(fit_x$results$par, fit_y$results$par,
      theta=copula_F$p_est, method="Frank")
fit_R(fit_x$results$par, fit_y$results$par,
      theta=copula_C$p_est, method="Clayton")
#####################non-parametric estimation of R #################
R_hat <- function(x, y) {
  sum(x < y) / length(x)
}
R_hat(x=X, y=Y)








