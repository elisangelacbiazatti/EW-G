#App3 - Regression modeling applied to COVID-19 data in Porto Alegre, Brazil

require(survival)
require(moments)
require(AdequacyModel)

dados<-read.csv(file.choose(), header = TRUE,sep=";")

tempo=dados$tempo

n=length(tempo)
n
is.numeric(tempo)

TTT(tempo, col="red", lwd=2, grid=TRUE, lty = 2)

boxplot(tempo,ylab="Time")
summary(tempo)
hist(tempo,xlab="Time",col="white",main = " ")

truehist(tempo,nbins = 20)

tabelam=round(cbind(mean(tempo), median(tempo), sd(tempo), var(tempo),skewness(tempo), kurtosis(tempo), min(tempo), max(tempo)),4)
tabelam

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
mlw<-survreg(Surv(y,censura)~1,dist = "extreme")
summary(mlw)

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
set.seed(1729)
v23 <- optim(c(mlw23$coefficients[1], mlw23$coefficients[2], mlw23$coefficients[3], mlw23$scale, 1,1),log.veroewgu23, method = "SANN", hessian=T)

v23$par 

tabela1=rbind(v23$par[1],v23$par[2],v23$par[3],v23$par[4],v23$par[5],v23$par[6])
rownames(tabela1)=c("beta_0", "beta_1","beta_2" ,"sigma","alpha","p")
colnames(tabela1)=("Estimativas")
tabela1

logVeroew23<-(-1)*v23$value
logVeroew23
inversa23<-solve(v23$hessian)
varpar23<-diag(inversa23)
erropad23<-sqrt(varpar23)
erropad23

#p-value = 2*(1-pnorm(abs(MLE/SE)))
Z0=v23$par[1]/erropad23[1]
Z0
p0=2*(1-pnorm(abs(Z0))) 
p0
Z1=v23$par[2]/erropad23[2]
Z1
p1=2*(1-pnorm(abs(Z1)))
p1
Z2=abs(v23$par[3])/erropad23[3]
Z2
p2=2*(1-pnorm(abs(Z2)))
p2

tabela2=cbind(c(erropad23[1],erropad23[2],erropad23[3],erropad23[4],erropad23[5],erropad23[6]),c(p0,p1,p2,NA,NA,NA))
colnames(tabela2)=c("Erro Padrão", "p-valor")
tabela2
tabs<-cbind(tabela1,tabela2)
tabs

#AIC, CAIC, BIC
np23<-6
AIC23<-(-2*logVeroew23)+(2*np23)  

AICc23<-AIC23 + (((2*np23^2)+(2*np23))/(n-np23-1))

BIC23<-(-2*logVeroew23) + np23*log(n)


medidasew23<-rbind(c(AIC23, AICc23, BIC23),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))
colnames(medidasew23)=c("AIC","AICc","BIC")
medidasew23

tabelaew23=cbind(tabela1,tabela2,medidasew23)
tabelaew23

#calculo do IC

CI.matrix <- as.data.frame(matrix(NA, nrow = 3, ncol = 6))

CI.matrix[1,] <- v23$par
CI.matrix[2,] <- v23$par - 1.96 * erropad23
CI.matrix[3,] <- v23$par + 1.96 * erropad23
names(CI.matrix) <- c("beta_0", "beta_1","beta_2" ,"sigma","alpha","p")
rownames(CI.matrix) <- c("MLE", "95% Lower bound", "95% Upper bound")

CI.matrix


################
# KwGu or log Kw weibull 
################
log.verokw23 <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  sigma <- par[4]
  a <- par[5]
  b <- par[6]
  
  
  mu <- beta0*x0 +  beta1*x1 + beta2*x2 
  w <- (y-mu)/sigma
  g=(1/sigma)*exp((-w)-exp(-w))
  G=exp(-exp(-w))
  
  df=a*b*g*G^(a-1)*(1-G^a)^(b-1)
  
  st<- (1-(G^a))^b
  
  lv <- censura*log(df)+(1-censura)*log(st)
  logv <- - sum(lv)
  if ((sigma>0)&&(a>0)&&(b>0))
    return(logv)
  else return (-Inf)
}

set.seed(1729)
vkw23 <- optim(c(mlw23$coefficients[1], mlw23$coefficients[2], mlw23$coefficients[3], mlw23$scale, 1,1),log.verokw23, method = "SANN", hessian=T)
vkw23$par

tabelakw23_1=rbind(vkw23$par[1],vkw23$par[2],vkw23$par[3],vkw23$par[4],vkw23$par[5],vkw23$par[6])
rownames(tabelakw23_1)=c("beta_0", "beta_1", "beta_2","sigma","a","b")
colnames(tabelakw23_1)=("Estimativas")
tabelakw23_1

logVerokw23<-(-1)*vkw23$value
logVerokw23
inversakw23<-solve(vkw23$hessian)
varparkw23<-diag(inversakw23)
erropadkw23<-sqrt(varparkw23)
erropadkw23

#p-value = 2*(1-pnorm(abs(MLE/SE)))
Z0kw23=vkw23$par[1]/erropadkw23[1]
Z0kw23
p0kw23=2*(1-pnorm(abs(Z0kw23)))
p0kw23
Z1kw23=vkw23$par[2]/erropadkw23[2]
Z1kw23
p1kw23=2*(1-pnorm(abs(Z1kw23)))
p1kw23
Z2kw23=vkw23$par[3]/erropadkw23[3]
Z2kw23
p2kw23=2*(1-pnorm(abs(Z2kw23)))
p2kw23

tabelakw23_2=cbind(c(erropadkw23[1],erropadkw23[2],erropadkw23[3],erropadkw23[4],erropadkw23[5],erropadkw23[6]),c(p0kw23,p1kw23,p2kw23,NA,NA,NA))
colnames(tabelakw23_2)=c("Erro Padrão", "p-valor")
tabelakw23_2
tabskw<-cbind(tabelakw23_1,tabelakw23_2)
tabskw

#AIC, CAIC, BIC
npkw23<-6
AICkw23<-(-2*logVerokw23)+(2*npkw23)  

AICckw23<-AICkw23 + (((2*npkw23^2)+(2*npkw23))/(n-npkw23-1))

BICkw23<-(-2*logVerokw23) + npkw23*log(n)



medidaskw23<-rbind(c(AICkw23, AICckw23, BICkw23),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))
colnames(medidaskw23)=c("AIC","AICc","BIC")
medidaskw23

tabelakw23=cbind(tabelakw23_1,tabelakw23_2,medidaskw23)
tabelakw23

#calculo do IC

CI.matrixKw <- as.data.frame(matrix(NA, nrow = 3, ncol = 6))

CI.matrixKw[1,] <- vkw23$par
CI.matrixKw[2,] <- vkw23$par - 1.96 * erropadkw23
CI.matrixKw[3,] <- vkw23$par + 1.96 * erropadkw23
names(CI.matrixKw) <- c("beta_0", "beta_1", "beta_2","sigma","a","b")
rownames(CI.matrixKw) <- c("MLE", "95% Lower bound", "95% Upper bound")

CI.matrixKw


###############
#log beta Weibull 
###############
log.verob23 <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  sigma <- par[4]
  a<-par[5]
  b<-par[6]
  
  
  mu <- beta0*x0 + beta1*x1 + beta2*x2 
  w <- (y-mu)/sigma
  
  g=(1/sigma)*exp((-w)-exp(-w))
  G=exp(-exp(-w))
  df<-g/(beta(a,b))*G^(a-1)*(1-G)^(b-1)
  st<-1-pbeta(G, a, b)
  
  
  lv <- censura*log(df)+(1-censura)*log(st)
  
  logv <- - sum(lv)
  if ((sigma>0)&&(a>0)&&(b>0))
    return(logv)
  else return (-Inf)
}
set.seed(1729)
vb23 <- optim(c(3.7,0.00978876,0.473619831,1.804604522,0.09,1),log.verob23, method = "SANN", hessian=T)

vb23$par

tabela1b23=rbind(vb23$par[1],vb23$par[2],vb23$par[3],vb23$par[4],vb23$par[5],vb23$par[6])

rownames(tabela1b23)=c("beta_0", "beta_1", "beta_2","sigma","a","b")
colnames(tabela1b23)=("Estimativas")
tabela1b23

logVerowb23<-(-1)*vb23$value
logVerowb23
inversabw23<-solve(vb23$hessian)
varparbw23<-diag(inversabw23)
erropadbw23<-sqrt(varparbw23)
erropadbw23 #erros padrão

#p-valor = 2*(1-pnorm(abs(MLE/SE)))
Z0bw23=abs(vb23$par[1])/erropadbw23[1]
Z0bw23
p0bw23=2*(1-pnorm(abs(Z0bw23)))
p0bw23
Z1bw23=abs(vb23$par[2])/erropadbw23[2]
Z1bw23
p1bw23=2*(1-pnorm(abs(Z1bw23)))
p1bw23
Z2bw23=abs(vb23$par[3])/erropadbw23[3]
Z2bw23
p2bw23=2*(1-pnorm(abs(Z2bw23)))
p2bw23

tabela2b23=cbind(c(erropadbw23[1],erropadbw23[2],erropadbw23[3],erropadbw23[4],erropadbw23[5],erropadbw23[6]),c(p0bw23,p1bw23,p2bw23,NA,NA,NA))
colnames(tabela2b23)=c("Erro Padrão", "p-valor")
tabela2b23

cbind(tabela1b23,tabela2b23)

#AIC, CAIC, BIC
npbw23<-6
AICbw23<-(-2*logVerowb23)+(2*npbw23)  

AICcbw23<-AICbw23 + (((2*npbw23^2)+(2*npbw23))/(n-npbw23-1))
#
BICbw23<-(-2*logVerowb23) + npbw23*log(n)
#


medidasbwei23<-rbind(c(AICbw23, AICcbw23, BICbw23),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))
colnames(medidasbwei23)=c("AIC","AICc","BIC")
medidasbwei23

tabelabw23=cbind(tabela1b23,tabela2b23,medidasbwei23)
tabelabw23

#calculo do IC

CI.matrixlbw <- as.data.frame(matrix(NA, nrow = 3, ncol = 6))

CI.matrixlbw[1,] <- vb23$par
CI.matrixlbw[2,] <- vb23$par - 1.96 * erropadbw23
CI.matrixlbw[3,] <- vb23$par + 1.96 * erropadbw23
names(CI.matrixlbw) <- c("beta_0", "beta_1", "beta_2","sigma","a","b")
rownames(CI.matrixlbw) <- c("MLE", "95% Lower bound", "95% Upper bound")

CI.matrixlbw


#####
tabelaew23
tabelakw23
tabelabw23

###############
# Vuong teste
###############
pdf_ewgu23 <- function(par,x){
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
}

z1 = pdf_ewgu23(v23$par,y)


pdf_kw23 <- function(par,x){
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  sigma <- par[4]
  a <- par[5]
  b <- par[6]
  
  mu <- beta0*x0 +  beta1*x1 + beta2*x2 
  w <- (y-mu)/sigma
  g=(1/sigma)*exp((-w)-exp(-w))
  G=exp(-exp(-w))
  
  df=a*b*g*G^(a-1)*(1-G^a)^(b-1)
} 

z2 = pdf_kw23(vkw23$par,y)

#LEW-W x KwGu
dif_log8 = log(z1) - log(z2)
n = length(y)
est_teste8 = (sqrt(n)^(-1) * sum(dif_log8)) * 
  (n^(-1) * sum((dif_log8)^2) - ( n^(-1)* sum(dif_log8))^2)^(-1)
est_teste8



pdf_b23 <- function(par,x){
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  sigma <- par[4]
  a<-par[5]
  b<-par[6]
 
  mu <- beta0*x0 + beta1*x1 + beta2*x2 
  w <- (y-mu)/sigma
  
  g=(1/sigma)*exp((-w)-exp(-w))
  G=exp(-exp(-w))
  df<-g/(beta(a,b))*G^(a-1)*(1-G)^(b-1)
}

z3 = pdf_b23(vb23$par,y)

#LEW-W x  LBW
dif_log9 = log(z1) - log(z3)
n = length(y)
est_teste9 = (sqrt(n)^(-1) * sum(dif_log9)) * 
  (n^(-1) * sum((dif_log9)^2) - ( n^(-1)* sum(dif_log9))^2)^(-1)
est_teste9

