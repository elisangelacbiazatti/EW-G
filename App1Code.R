#App1 - COVID-19 Data in Goiania - Brazil

library(AdequacyModel)
library(MASS)
require(GenSA)
require(survival)
require(moments)

data<-read.csv(file.choose(), header = TRUE,sep=";")

tempo=data$time

n=length(tempo)
n
is.numeric(tempo)

TTT(tempo, col="red", lwd=2, grid=TRUE, lty = 2)

boxplot(tempo,ylab="Time")
summary(tempo)
hist(tempo,xlab="Time",col="white",main = " ",breaks = 20)

truehist(tempo)

tabela=round(cbind(mean(tempo), median(tempo), sd(tempo), var(tempo),skewness(tempo), kurtosis(tempo), min(tempo), max(tempo)),4)
tabela


## Fit a 4 parameter distribution to the data
fit.sa4<- function(data,density) {
  minusllike<-function(x) -sum(log(density(data,x[1],x[2],x[3],x[4]))) #minus the loglik  
  lower <- c(0,0,0,0) 
  upper <- c(100,100,100,100)
  out <- GenSA(lower = lower, upper = upper, fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}


###############################
#EW-BS
###############################
pdf.ewbs<-function(x,alpha,p,a,beta){
  lambda<-1
  G<-pnorm((1/a)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(a^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*a^2))/(2*a*sqrt(2*pi*beta))
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,pdf.ewbs)

pdf_ewbs = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  beta = par[4]
  G<-pnorm((1/a)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(a^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*a^2))/(2*a*sqrt(2*pi*beta))
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
}
cdf_ewbs = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  beta = par[4]
  G<-pnorm((1/a)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(a^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*a^2))/(2*a*sqrt(2*pi*beta))
  fa=(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha)
}
set.seed(1729)
resultsEWbs  = goodness.fit(pdf = pdf_ewbs , cdf = cdf_ewbs ,
                            starts = c(5, 0.1,  1,  1),
                            data = tempo, method = "SANN", domain = c(0, Inf),
                            mle = NULL);resultsEWbs

############
#EW-Gamma
###########
pdf.ewga<-function(x,alpha,p,a,b){
  lambda<-1
  g=dgamma(x,a,b)
  G=pgamma(x,a,b)
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,pdf.ewga)

pdf_ewga = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  b = par[4]
  g=dgamma(x,a,b)
  G=pgamma(x,a,b)
  fd=((alpha*p*g)/((1-G)*lambda^p))*(-log(1-G))^(p-1)*exp(-lambda^(-p)*(-log(1-G))^p)*(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha-1)
}
cdf_ewga = function(par, x){
  alpha = par[1]
  lambda = 1
  p = par[2]
  a = par[3]
  b = par[4]
  g=dgamma(x,a,b)
  G=pgamma(x,a,b)
  fa=(1-exp(-lambda^(-p)*(-log(1-G))^p))^(alpha)
}
set.seed(1729)
resultsEWga  = goodness.fit(pdf = pdf_ewga , cdf = cdf_ewga ,
                            starts = c(5.01348061, 0.09324391, 6.76007528, 0.43059929),
                            data = tempo, method = "SANN", domain = c(0, Inf),
                            mle = NULL);resultsEWga


#######################################################
#Kw-BS
#######################################################
kbs.pdf=function(x,a,b,alpha,beta){
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,kbs.pdf)

pdf_Kwbs = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  a*b*g*G^(a-1)*(1-G^a)^(b-1)
}

cdf_Kwbs = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwbs = goodness.fit(pdf = pdf_Kwbs, cdf = cdf_Kwbs,
                           starts = c(6.471824, 55.686099,  4.560596,  8.861920), 
                           data = tempo, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsKwbs
chutes<-c(6.471824, 55.686099,  4.560596,  8.861920)
fit.gtnhkbs<- fitdistr(tempo,kbs.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-6,1e-6),maxit=10000))
fit.gtnhkbs
#################
# Beta BS
#################
bbs.pdf=function(x,a,b,alpha,beta){
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}

set.seed(1729)
fit.sa4(tempo,bbs.pdf)

cdfbbs=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  pbeta(G, a, b)
}
pdfbbs=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  G<-pnorm((1/alpha)*((x/beta)^(0.5)-(x/beta)^(-0.5)))
  g<-exp(alpha^(-2))*x^(-1.5)*(x+beta)*exp(-((x/beta)+(beta/x))/(2*alpha^2))/(2*alpha*sqrt(2*pi*beta))
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbbs = goodness.fit(pdf = pdfbbs, cdf = cdfbbs,
                          starts = c(37.489633, 38.423280,  7.782113,  9.132127),
                          data = tempo, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsbbs

############
chutes<-c(37.489633, 38.423280,  7.782113,  9.132127)
fit.gtnhbbs<- fitdistr(tempo,bbs.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-6,1e-6,1e-6,1e-6),maxit=10000))
fit.gtnhbbs



#######################################################
#Kw-Gamma
#######################################################
kga.pdf=function(x,a,b,alpha,beta){
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fd<-a*b*g*G^(a-1)*(1-G^a)^(b-1)
  fd
}
set.seed(1729)
fit.sa4(tempo,kga.pdf)

pdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  a*b*g*G^(a-1)*(1-G^a)^(b-1)
  
}

cdf_Kwga = function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  beta = par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fa<-1-(1-G^a)^b
}
set.seed(1729)
resultsKwga = goodness.fit(pdf = pdf_Kwga, cdf = cdf_Kwga,
                           starts = c( 13,  8,  0.09,  0.008), 
                           data = tempo, method = "SANN", domain = c(0, Inf),
                           mle = NULL);resultsKwga

chutes<-c(13,  8,  0.09,  0.008)
fit.gtnhkga<- fitdistr(tempo,kga.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-8,1e-6,1e-6,1e-6),maxit=10000))
fit.gtnhkga

#################
# Beta Gamma
#################
bga.pdf=function(x,a,b,alpha,beta){
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  fd=(1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1) 
  fd
}
set.seed(1729)
fit.sa4(tempo,bga.pdf)

cdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  pbeta(G, a, b)
}
pdfbga=function(par,x)
{
  a=par[1]
  b=par[2]
  alpha=par[3]
  beta=par[4]
  g=dgamma(x,alpha,beta)
  G=pgamma(x,alpha,beta)
  (1/beta(a,b))*g*(G^(a-1))*(1-G)^(b-1)
}
set.seed(1729)
resultsbga = goodness.fit(pdf = pdfbga, cdf = cdfbga,
                          starts = c(0.1, 9,  0.19453,  0.0777195),
                          data = tempo, method = "SANN", domain = c(0, Inf),
                          mle = NULL);resultsbga

############
chutes<-c(0.1, 9,  0.19453,  0.0777195)
fit.gtnhbga<- fitdistr(tempo,bga.pdf,start=list(a=chutes[1],b=chutes[2],alpha=chutes[3], beta=chutes[4]),control=list(ndeps=c(1e-12,1e-12,1e-12,1e-10),maxit=10000))
fit.gtnhbga


#densities
truehist(tempo,
         ylim=c(0,.065),
         col = "white",ylab="f(x)",xlab = "x")
curve(pdf_ewbs(resultsEWbs$mle,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(pdf_ewga(resultsEWga$mle,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(pdfbbs(fit.gtnhbbs$estimate,x),add=TRUE, lwd = 3, col="green",lty=3)
curve(pdfbga(fit.gtnhbga$estimate,x),add=TRUE, lwd = 3, col=7,lty=4)

legend(65,.06, legend = c("EW-BS", "EW-Ga", "B-BS", "B-Ga"),
       col = c("red","blue","green",7), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")

#CDFs
km<- survfit(Surv(tempo) ~ 1) #Kaplan-Meier
plot(km$time, 1-km$surv, xlab = "x", ylab="F(x)",lwd=2,type = "l")#empirical
curve(cdf_ewbs(resultsEWbs$mle,x),add=TRUE, lwd = 3, col="red",lty=1)
curve(cdf_ewga(resultsEWga$mle,x),add=TRUE, lwd = 3, col="blue",lty=2)
curve(cdfbbs(fit.gtnhbbs$estimate,x),add=TRUE, lwd = 3, col="green",lty=3)
curve(cdfbga(fit.gtnhbga$estimate,x),add=TRUE, lwd = 3, col=7,lty=4)

legend(65,.65, legend = c( "EW-BS", "EW-Ga", "B-BS", "B-Ga"),
       col = c("red","blue","green",7), lwd = 3 ,lty = c(1,2,3,4) , bty ="n")

################
# Vuong test
################

#EW-Ga x EW-BS
z1 = pdf_ewga(resultsEWga$mle,tempo)
z2 = pdf_ewbs(resultsEWbs$mle,tempo)
dif_log = log(z1) - log(z2)
n = length(tempo)
est_teste1 = (sqrt(n)^(-1) * sum(dif_log)) * 
  (n^(-1) * sum((dif_log)^2) - ( n^(-1)* sum(dif_log))^2)^(-1)
est_teste1

#EW-BS x BGa
z1 = pdf_ewbs(resultsEWbs$mle,tempo)
z3 = pdfbga(fit.gtnhbga$estimate,tempo)
dif_log2 = log(z1) - log(z3)
n = length(tempo)
est_teste2 = (sqrt(n)^(-1) * sum(dif_log2)) * 
  (n^(-1) * sum((dif_log2)^2) - ( n^(-1)* sum(dif_log2))^2)^(-1)
est_teste2

#EW-BS x BBS
z1 = pdf_ewbs(resultsEWbs$mle,tempo)
z4 = pdfbbs(fit.gtnhbbs$estimate,tempo)
dif_log3 = log(z1) - log(z4)
n = length(tempo)
est_teste3 = (sqrt(n)^(-1) * sum(dif_log3)) * 
  (n^(-1) * sum((dif_log3)^2) - ( n^(-1)* sum(dif_log3))^2)^(-1)
est_teste3

#EW-BS x Kw-BS
z1 = pdf_ewbs(resultsEWbs$mle,tempo)
z31 = pdf_Kwbs(fit.gtnhkbs$estimate,tempo)
dif_log31 = log(z1) - log(z31)
n = length(tempo)
est_teste31 = (sqrt(n)^(-1) * sum(dif_log31)) * 
  (n^(-1) * sum((dif_log31)^2) - ( n^(-1)* sum(dif_log31))^2)^(-1)
est_teste31

#EW-BS x Kw-Ga
z1 = pdf_ewbs(resultsEWbs$mle,tempo)
z41 = pdf_Kwga(fit.gtnhkga$estimate,tempo)
dif_log41 = log(z1) - log(z41)
n = length(tempo)
est_teste41 = (sqrt(n)^(-1) * sum(dif_log41)) * 
  (n^(-1) * sum((dif_log41)^2) - ( n^(-1)* sum(dif_log41))^2)^(-1)
est_teste41
