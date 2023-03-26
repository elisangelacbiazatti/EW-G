#App2 - COVID-19 Data in Natal - Brazi

library(AdequacyModel)
library(MASS)
require(GenSA)
require(survival)
require(moments)

dados<-read.csv(file.choose(), header = TRUE,sep=";")

tempo=dados$tempo

n=length(tempo)
n
is.numeric(tempo)

TTT(tempo, col="red", lwd=2, grid=TRUE, lty = 2)

boxplot(tempo,ylab="Time")
summary(tempo)
hist(tempo,xlab="Time",col="white",main = " ",breaks = 10)

truehist(tempo)

tabela=round(cbind(mean(tempo), median(tempo), sd(tempo), var(tempo),skewness(tempo), kurtosis(tempo), min(tempo), max(tempo)),4)



## Fit a 4 parameter distribution to the data
fit.sa4<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3],x[4]))) #minus the loglik  
  lower <- c(0,0,0,0) 
  upper <- c(100,100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}

#############################
#EW-W
#############################

pdf.eww<-function(x,alpha,p,k,beta){
  lambda<-1
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,pdf.eww)

pdf_eww = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
}
cdf_eww = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fa=(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha)
}
set.seed(1729)
resultsEWW  = goodness.fit(pdf = pdf_eww , cdf = cdf_eww ,
                           starts = c(1.064926, 2.236260, 0.532497, 9.734294),
                           data = tempo, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsEWW

chutes<-c(1.064926, 2.236260, 0.532497, 9.734294)
fit.gtnheww<- fitdistr(tempo,pdf.eww,start=list(alpha=chutes[1],p=chutes[2],k=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-10,1e-12),maxit=10000))
fit.gtnheww 

#############################
#EW-LL
#############################

pdf.ewll<-function(x,alpha,p,a,beta){
  lambda<-1
  g<-(beta/a^beta)*(x^(beta-1))*(1+(x/a)^beta)^(-2)
  G<-1-(1/(1+(x/a)^beta))
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,pdf.ewll)

pdf_ewll = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  beta = par[4]
  g<-(beta/a^beta)*(x^(beta-1))*(1+(x/a)^beta)^(-2)
  G<-1-(1/(1+(x/a)^beta))
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
}
cdf_ewll = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  beta = par[4]
  g<-(beta/a^beta)*(x^(beta-1))*(1+(x/a)^beta)^(-2)
  G<-1-(1/(1+(x/a)^beta))
  fa=(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha)
}
set.seed(1729)
resultsEWll  = goodness.fit(pdf = pdf_ewll , cdf = cdf_ewll ,
                            starts =c( 0.5927415, 5.0980851, 4.6372077, 0.4990583),
                            data = tempo, method = "SANN", domain = c(0, Inf),
                            mle = NULL);resultsEWll
chutes<-c(0.5927415, 5.0980851, 4.6372077, 0.4990583)
fit.gtnhewll<- fitdistr(tempo,pdf.ewll,start=list(alpha=chutes[1],p=chutes[2],a=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-12,1e-6),maxit=10000))
fit.gtnhewll

########################################################
#Kw-Wei 
#######################################################
kwei.pdf=function(x,a,b,k,beta){
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,kwei.pdf)

pdf_KwWei = function(par, x){
  a = par[1]
  b = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  a*b*g*G^(a-1)*(1-G^a)^(b-1)
}

cdf_KwWei = function(par, x){
  a = par[1]
  b = par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k) 
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwWei = goodness.fit(pdf = pdf_KwWei, cdf = cdf_KwWei,
                            starts = c(4.6641747, 20,   0.3312148, 20), 
                            data = tempo, method = "SANN", domain = c(0, Inf),
                            mle = NULL);resultsKwWei

chutes<-c(4.6641747, 20,   0.3312148, 20)
fit.gtnhkw<- fitdistr(tempo,kwei.pdf,start=list(a=chutes[1],b=chutes[2],k=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-9,1e-12,1e-6,1e-9),maxit=10000))
fit.gtnhkw


#################
# Beta Weibull
#################
bw.pdf=function(x,a,b,k,beta){
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}

set.seed(1729)
fit.sa4(tempo,bw.pdf)

cdfbw=function(par,x)
{
  a=par[1]
  b=par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  pbeta(G, a, b)
}
pdfbw=function(par,x)
{
  a=par[1]
  b=par[2]
  k = par[3]
  beta = par[4]
  g <- (k/beta)*(x/beta)^(k-1)*exp(-(x/beta)^k)
  G <- 1-exp(-(x/beta)^k)
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbw = goodness.fit(pdf = pdfbw, cdf = cdfbw,
                         starts = c(1.107557,  15,   1.159858,  20),
                         data = tempo, method = "SANN", domain = c(0, Inf),
                         mle = NULL);resultsbw

chutes<-c(1.107557,  15,   1.159858,  20)
fit.gtnhbw<- fitdistr(tempo,bw.pdf,start=list(a=chutes[1],b=chutes[2],k=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-8,1e-12),maxit=10000))
fit.gtnhbw



########################################################
#Kw-LL
#######################################################

pdfkwll<-function(x,a,b,alpha,beta){
  
  g<-(beta/alpha^beta)*(x^(beta-1))*(1+(x/alpha)^beta)^(-2)
  G<-1-(1/(1+(x/alpha)^beta))
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
}
set.seed(1729)
fit.sa4(tempo,pdfkwll)

pdf_KwLL = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g<-(beta/alpha^beta)*(x^(beta-1))*(1+(x/alpha)^beta)^(-2)
  G<-1-(1/(1+(x/alpha)^beta))
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
}
cdf_KwLL = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g<-(beta/alpha^beta)*(x^(beta-1))*(1+(x/alpha)^beta)^(-2)
  G<-1-(1/(1+(x/alpha)^beta))
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwLL = goodness.fit(pdf = pdf_KwLL, cdf = cdf_KwLL,
                           starts = c(0.2835631, 12.4882136, 9.3490248,  4.1553996), 
                           data = tempo, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsKwLL

chutes<-c(0.2835631, 12.4882136, 9.3490248,  4.1553996)
fit.gtnhkwll<- fitdistr(tempo,pdfkwll,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-12,1e-12,1e-6,1e-12),maxit=10000))
fit.gtnhkwll

#################
# Beta LL
#################
bll.pdf=function(x,a,b,alpha,beta){
  g<-(beta/alpha^beta)*(x^(beta-1))*(1+(x/alpha)^beta)^(-2)
  G<-1-(1/(1+(x/alpha)^beta))
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
set.seed(1729)
fit.sa4(tempo,bll.pdf)

cdfbll=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  G<-1-(1/(1+(x/alpha)^beta))
  pbeta(G, a, b)
}
pdfbll=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  g<-(beta/alpha^beta)*(x^(beta-1))*(1+(x/alpha)^beta)^(-2)
  G<-1-(1/(1+(x/alpha)^beta))
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbll = goodness.fit(pdf = pdfbll, cdf = cdfbll,
                          starts = c(0.92602,  19.64910, 20 ,  1.33090), 
                          data = tempo, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsbll

chutes<-c(0.92602,  19.64910, 20 ,  1.33090)
fit.gtnhbll<- fitdistr(tempo,bll.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-8,1e-12,1e-6),maxit=10000))
fit.gtnhbll

#densities
truehist(tempo,nbins = 10,
         ylim=c(0,.08),
         col = "white",ylab="f(x)",xlab = "x")
curve(pdf_ewll(fit.gtnhewll$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(pdf_eww(fit.gtnheww$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(pdf_KwWei(fit.gtnhkw$estimate,x),add=TRUE, lwd = 3, col="green",lty=3)
curve(pdfbw(fit.gtnhbw$estimate,x),add=TRUE, lwd = 3, col=7,lty=4)

legend(15,.08, legend = c("EW-LL", "EW-W", "Kw-W", "B-W"),
       col = c("red","blue","green",7), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")


#cdfs
km<- survfit(Surv(tempo) ~ 1) #Kaplan-Meier
plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l")#cdf emprirical
curve(cdf_ewll(fit.gtnhewll$estimate,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(cdf_eww(fit.gtnheww$estimate,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(cdf_KwWei(fit.gtnhkw$estimate,x),add=TRUE, lwd = 3, col="green",lty=3)
curve(cdfbw(fit.gtnhbw$estimate,x),add=TRUE, lwd = 3, col=7,lty=4)

legend(25,.8, legend = c( "EW-LL", "EW-W", "Kw-W", "B-W"),
       col = c("red","blue","green",7), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")


################
# Vuong test
################

##########Teste de Vuong
#EW-LL x EW-W
z1 = pdf_ewll(fit.gtnhewll$estimate,tempo)
z2 = pdf_eww(fit.gtnheww$estimate,tempo)
dif_log = log(z1) - log(z2)
n = length(tempo)
est_teste1 = (sqrt(n)^(-1) * sum(dif_log)) * 
  (n^(-1) * sum((dif_log)^2) - ( n^(-1)* sum(dif_log))^2)^(-1)
est_teste1

#EW-ll x Bw
z1 = pdf_ewll(fit.gtnhewll$estimate,tempo)
z3 = pdfbw(fit.gtnhbw$estimate,tempo)
dif_log2 = log(z1) - log(z3)
n = length(tempo)
est_teste2 = (sqrt(n)^(-1) * sum(dif_log2)) * 
  (n^(-1) * sum((dif_log2)^2) - ( n^(-1)* sum(dif_log2))^2)^(-1)
est_teste2

#EW-ll x KwW
z1 = pdf_ewll(fit.gtnhewll$estimate,tempo)
z4 = pdf_KwWei(fit.gtnhkw$estimate,tempo)
dif_log3 = log(z1) - log(z4)
n = length(tempo)
est_teste3 = (sqrt(n)^(-1) * sum(dif_log3)) * 
  (n^(-1) * sum((dif_log3)^2) - ( n^(-1)* sum(dif_log3))^2)^(-1)
est_teste3



