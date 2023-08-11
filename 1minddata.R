setwd(" ")
BMO<-as.matrix(read.csv("  .csv"))
d5price<-BMO[,c(2,4,6)]
install.packages("astsa")
require(astsa)
install.packages("cubature")
install.packages("quantmod")
install.packages("fArma")
install.packages("fGarch")

library( quantmod )
library( fArma )
require(cubature)
require(fGarch)

plot.ts(BMO[,2])
plot.ts(BMO[,4])
plot.ts(BMO[,6])
#acf2(BMO[,2])

price<-as.numeric(BMO[1:1890,2])
bid<-as.numeric(BMO[1:1890,4])
ask<-as.numeric(BMO[1:1890,6])

retp<-function(t){((price[t]-price[t-1])/price[t-1])*100}
retup<-as.ts(retp(2:1890))
bir<-function(t){((bid[t]-bid[t-1])/bid[t-1])*100}
bidr<-bir(2:1890)
asr<-function(t){((ask[t]-ask[t-1])/ask[t-1])*100}
askr<-asr(2:1890)
acf2(askr)
acf2(bidr)
#####################################################################
armaSearch = function(
xx,
minOrder=c(0,0),
maxOrder=c(5,5),
trace=FALSE )
{
bestAic = 1e9
len = NROW( xx )
for( p in minOrder[1]:maxOrder[1] ) for( q in minOrder[2]:maxOrder[2] )
{
if( p == 0 && q == 0 )
{
next
}
 
formula = as.formula( paste( sep="", "xx ~ arma(", p, ",", q, ")" ) )
 
fit = tryCatch( armaFit( formula, data=xx ),
error=function( err ) FALSE,
warning=function( warn ) FALSE )
if( !is.logical( fit ) )
{
fitAic = fit@fit$aic
if( fitAic < bestAic )
{
bestAic = fitAic
bestFit = fit
bestModel = c( p, q )
}
 
if( trace )
{
ss = paste( sep="", "(", p, ",", q, "): AIC = ", fitAic )
print( ss )
}
}
else
{
if( trace )
{
ss = paste( sep="", "(", p, ",", q, "): None" )
print( ss )
}
}
}
 
if( bestAic < 1e9 )
{
return( list( aic=bestAic, fit=bestFit, model=bestModel ) )
}
 
return( FALSE )
} 
####################GARCH############################################
pp= price
require("astsa")
acf2(pp)
pt = as.ts( tail( pp, 500 ) )
armaSearch(pp)
ppArma = armaFit( pp ~ arma(1,2), data=pp )
ppArma@fit$aic
b<-as.numeric( predict( ppArma, n.ahead=1, doplot=F )$pred )
#as.numeric(BMO[1871,2])-price[1870]
ppGarch = garchFit(~arma(1,2) + garch(1,2), data=as.ts(tail(pp, 1878)))

predict(ppGarch, n.ahead=1, doplot=F)

garchFit(~garch(1,1),data=pp)
varp<-(pp-mean(pp))^2
res<-residuals(ppArma)
omega<-6.304e-04 
alpha1<-8.357e-01
beta1<- 1.701e-01 
#beta2<- 1.701e-01
ht1<-function(t){as.numeric(omega+alpha1*res[t]^2+beta1*varp[t])}    
+beta2*varp[t-1])}
#####################################################################
armaSearch(retup)
ppArma = armaFit( retup ~ arma(3,3), data=retup )
ppArma@fit$aic
b<-as.numeric( predict( ppArma, n.ahead=1, doplot=F )$pred )

ppGarch = garchFit(~arma(3,3) + garch(1,1), data=retup)
predict(ppGarch, n.ahead=1, doplot=F)

pba<-as.matrix(cbind(retup,bid[-1],ask[-1]))
pba<-ts(as.matrix(cbind(retup,bidr,askr)))
plot(pba)
plot.ts(pba[1:200,], plot.type = "single", lty=c(1:3),col = 1:3,ylab="Returns")

ts.plot(as.ts(price),as.ts(bid),as.ts(ask),gpars=list(xlab="time", ylab="money", lty=c(1:3),col=c(1:3)))
ts.plot(as.ts(retup),as.ts(bidr),as.ts(askr),gpars=list(xlab="time", ylab="money", lty=c(1:3),col=c(1:3)))

install.packages("MTS")
library(MTS)
archTest(pba)
ccm(pba,level=T,output=T)

m1=comVol(pba,p=2)
names(m1)

dccPre(pba,p=2)

dpba<-diffM(pba,d=1)

Eccm(pba, maxp = 5, maxq = 6, include.mean = FALSE, rev = TRUE) #(1,3)


####Error-Correction VAR Models start here page 14 of MTS pdf 
n1<-ECMvar(pba,2,ibeta=alpha)
names(n1)
######################################################################
acf2(price)
acf2(bid)
acf2(ask)
y<-price[4:1890]
w<-price[3:1889]
#v<-price[2:1888]
#u<-price[1:1887]
x<-bid[4:1890]
z<-ask[4:1890]
acf2(x)
acf2(y)
acf2(z)
#eh<-c(varp[2],ht1(2:1889))
eh<-ht1(4:1890)
n<-1887
l1<-c(1:1887)
kh<-function(x,h){(1/h)*exp(-0.5*(x^2/h^2))}
khAR2<-function(l,j,h){kh(y[l]-y[j],h)}
nkhAR2<-function(l,j,h){(1/h)*exp(-0.5*((x[l]-x[j])^2+(z[l]-z[j])^2+(w[l]-w[j])^2+(eh[l]-eh[j]^2))/h^2)}
                                        +(u[l]-u[j]^2)+(v[l]-v[j]^2))/h^2)} 
AR2Newr<-function(l,j,h){nkhAR2(l,j,h[2])*khAR2(l,j,h[1])*(1/(n-1))}
AR2Newr1<-function(l,h){log((sum(AR2Newr(l,l1,h)))-AR2Newr(l,l,h))}
AR2Newr2<-function(h){(1/n)*Reduce("+",lapply(l1,AR2Newr1,h=h))}
AR2Newr3<-function(h){-AR2Newr2(h)}
optim(c(1,1),AR2Newr3,control=list(trace=1,maxit=300))
#h<-c(0.01596412,0.06181732) #AR(1)
#h<-c(0.07746348,0.02516263) #AR(3)
h<-c(0.01038926,0.23234213) #non-differenced data 
mconst<-1/sqrt(2*pi)
fdeltan<-function(s){mconst*(sum(kh(s-y,h[1])*kh(sqrt((z[length(l1)]-z)^2+(x[length(l1)]-x)^2+(eh[length(l1)]-eh)^2),h[2])))}
                                                      


+(w[length(l1)]-w)^2+(v[length(l1)]-v)^2+(u[length(l1)]-u)^2),h[2])))}
fdeltad<-sum(kh(sqrt((z[length(l1)]-z)^2+(x[length(l1)]-x)^2+(eh[length(l1)]-eh)^2),h[2]))
                     +(w[length(l1)]-w)^2+(v[length(l1)]-v)^2+(u[length(l1)]-u)^2),h[2]))
fdelta<-function(s){fdeltan(s)/fdeltad}
fd<-Vectorize(function(s){fdelta(s)})
zyx<-curve(fd(x),80,83,col=2)
plot(density(y))
lines(zyx,type="l",col=2)
adaptIntegrate(fd,-1,1)
##as.numeric(BMO[1891,2])
