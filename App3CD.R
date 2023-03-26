#Cook Distance for COVID-19 data in Porto Alegre

require(survival)

dados<-read.csv(file.choose(), header = TRUE,sep=";")

tempo=dados$tempo

n=length(tempo)
n
is.numeric(tempo)

censura=dados$censura

km<- survfit(Surv(tempo,censura) ~ 1) #Kaplan-Meier
plot(km,conf.int = F, xlab = "x", ylab=expression(hat(S)*"(x)"),mark.time = T)
summary(km)

x0=dados$x0
x1=dados$idade
x2=dados$uti

y<-log(tempo)

mlw23<-survreg(Surv(y,censura)~x1+factor(x2), dist="extreme") 
summary(mlw23) 


b0coef<- matrix(0,n,1)
b1coef<- matrix(0,n,1)
b2coef<- matrix(0,n,1)
sigcoef<- matrix(0,n,1)
alphacoef<- matrix(0,n,1)
pcoef<- matrix(0,n,1)

for(l in 1:n){
  dad<-dados[-l, ]
  tempo<-dad[ ,2]
  censura<-dad[ ,3]
  x0<-dad[ ,4]
  x1<-dad[ ,5]
  x2<-dad[ ,6]
  n<-length(tempo)
  
  ################
  #LEW-W
  ################
  log.veroewgu23 <- function(par){
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    sigma <- par[4]
    alpha <- par[5]
    p <- par[6]
    lambda <- 1
    
    mu <- beta0*x0 +beta1*x1 + beta2*x2
    w <- (y-mu)/sigma
    
    g=(1/sigma)*exp((-w)-exp(-w))
    G=exp(-exp(-w))
    
    df=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
    st= 1-(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha)
    
    lv <- censura*log(df)+(1-censura)*log(st)
    logv <- - sum(lv)
    if ((alpha>0)&&(sigma>0)&&(p>0))
      return(logv)
    else return (-Inf)
  }
  
  v23g <- optim(c(mlw23$coefficients[1], mlw23$coefficients[2], mlw23$coefficients[3], mlw23$scale, 1,1),log.veroewgu23, method = "SANN", hessian=T)
  
  b0coef[l]<- v23g$par[1]
  b1coef[l]<- v23g$par[2]
  b2coef[l]<- v23g$par[3]
  sigcoef[l]<- v23g$par[4]
  alphacoef[l]<- v23g$par[5]
  pcoef[l]<- v23g$par[6]
  
}

ck<-cbind(b0coef,b1coef,b2coef,sigcoef,alphacoef,pcoef)

inversag<-solve(v23g$hessian)
hes<-v23g$hessian
DC<-matrix(0,n,1)
DC1<-matrix(0,n,1)
for (u in 1:n){
  DC[u]<-t(ck[u,]-v23g$par)%*%inversag%*%(ck[u,]-v23g$par)
  DC1[u]<-t(ck[u,]-v23g$par)%*%hes%*%(ck[u,]-v23g$par)
}

#plots
#DC1
indice<-c(1:n)
plot(indice,DC1,type = "l",xlab="Index", ylab="Cook Distance")
identify(indice,DC1)
#DC
indice<-c(1:n)
plot(indice,DC,type = "l",xlab="Index", ylab="Cook Distance")
identify(indice,DC)

