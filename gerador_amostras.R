

library(extRemes)

#derivada parcial da copula
c_1 = function(u1, v1, p1){
  exp(-
        (((-log(u1))^(1/p1) + (-log(v1))^(1/p1)))^(p1))*
    (((-log(u1))^(1/p1)) +(-log(v1))^(1/p1))^(p1-1)*((-log(v1))^(1/p1 -1))*(1/v1)
}
#c_1(0.4)

#copula 2 Frank
c_2 = function(u1,v1,p1){
  
  num1=(exp(-p1*u1) - 1)*(exp(-p1*v1) - 1)
  den1= (exp(-p1) - 1) 
  
  a1=(1 + num1/den1)^(-1)

  num2=(exp(-p1*u1) - 1)*(exp(-p1*v1))
  den2=(exp(-p1) - 1)
  
  out= a1*num2/den2
  return(out)
}
#c_2(0.4,0.2,0.8)

#copula 3 Clayton
c_3 = function(u1,v1,p1){
  (v1^(-p1-1))*((u1^(-p1) + v1^(-p1) - 1)^(-1/p1 -1))
}
#c_3(0.4,0.2,0.8)

#funcao inversa da derivada parcial
c1_inv=function(y,v1, p1){
  uniroot(function(x){c_1(x,v1, p1) - y}, c(0,1),extendInt="yes", trace=1)$root
}

c2_inv=function(y,v1, p1){
  uniroot(function(x){c_2(x,v1, p1) - y}, c(0,1),extendInt="yes", trace=1)$root
}

c3_inv=function(y,v1, p1){
  uniroot(function(x){c_3(x,v1, p1) - y}, c(0,1),extendInt="yes", trace=1)$root
}

# n=10#tamanho de amostra
# 
# #parametros
# lx=0;sx=1; shx=1
# ly=0.5;sy=1.5; shy=1.5
# p1=0.5#parametro de dependencia

#gerando amostra da copula Gumbel-Hougard
rcgh<-function(n=100,lx=0,sx=1, shx=1,
               ly=0,sy=1, shy=1,p1=0.5){
  p1_0=p1
  
  #passo 1
  U2=runif(n)
  # u1=NULL
  # for(i in 1:n){
  #   v1=u2[i]
  #   #passo 2: derivada parcial
  #   #passo 3: inversa da derivada parcial
  #   
  #   #Passo 4: gerar v
  #   v=runif(1)
  #   #Passo 5: inversa da c1 aplicado em v
  #   u1[i]=c1_inv(v)
  # }
  
  f_apply_i=function(i){
    v1_2=U2[i]
    #passo 2: derivada parcial
    #passo 3: inversa da derivada parcial
    
    #Passo 4: gerar v
    v=runif(1)
    #Passo 5: inversa da c1 aplicado em v
    U1=c1_inv(v, v1_2, p1=p1_0)
    return(U1)
  }
  U1=sapply(1:n, f_apply_i)
  
  
  
  
  
  
  #Passo 6: inversa da CDF da GEV
  out=data.frame(x=qevd(U1, loc=lx, scale=sx, shape = shx, type = "GEV"),
                 y=qevd(U2, loc=ly, scale=sy, shape = shy, type = "GEV"))
  
  return(out)
}

#Exemplo
rcgh(n=100,lx=0,sx=1, shx=1,
     ly=0.5,sy=1.5, shy=1.5,p1=0.5)



#gerando amostra da copula Frank
rcf<-function(n=100,lx=0,sx=1, shx=1,
               ly=0,sy=1, shy=1,p1=0.5){
  p1_0=p1
  
  #passo 1
  U2=runif(n)

  f_apply_i=function(i){
    v1_2=U2[i]
    #passo 2: derivada parcial
    #passo 3: inversa da derivada parcial
    
    #Passo 4: gerar v
    v=runif(1)
    #Passo 5: inversa da c1 aplicado em v
    U1=c2_inv(v, v1_2, p1=p1_0)
    return(U1)
  }
  U1=sapply(1:n, f_apply_i)
  
  
  
  
  
  
  #Passo 6: inversa da CDF da GEV
  out=data.frame(x=qevd(U1, loc=lx, scale=sx, shape = shx, type = "GEV"),
                 y=qevd(U2, loc=ly, scale=sy, shape = shy, type = "GEV"))
  
  return(out)
}
#Exemplo
rcf(n=100,lx=0,sx=1, shx=1,
     ly=0.5,sy=1.5, shy=1.5,p1=0.5)


#gerando amostra da copula Frank
rcc<-function(n=100,lx=0,sx=1, shx=1,
              ly=0,sy=1, shy=1,p1=0.5){
  p1_0=p1
  
  #passo 1
  U2=runif(n)
  
  f_apply_i=function(i){
    v1_2=U2[i]
    #passo 2: derivada parcial
    #passo 3: inversa da derivada parcial
    
    #Passo 4: gerar v
    v=runif(1)
    #Passo 5: inversa da c1 aplicado em v
    U1=c3_inv(v, v1_2, p1=p1_0)
    return(U1)
  }
  U1=sapply(1:n, f_apply_i)
  
  
  
  
  
  
  #Passo 6: inversa da CDF da GEV
  out=data.frame(x=qevd(U1, loc=lx, scale=sx, shape = shx, type = "GEV"),
                 y=qevd(U2, loc=ly, scale=sy, shape = shy, type = "GEV"))
  
  return(out)
}
#Exemplo
rcc(n=100,lx=0,sx=1, shx=1,
    ly=0.5,sy=1.5, shy=1.5,p1=0.5)
