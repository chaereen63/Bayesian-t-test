SD=function(group1,		# data vector of group 1
		group2=0,		# data vector of group 2
		iters=1000000,	# the amount of iterations
		burns=800001,	# the amount of burnin
		chains=1, 		# the amount of chains
		thins=1, 		# the thinning
		sample,		# 1=one sample 2=two sample
		sig=2,		# both groups have (1) sigma or (2) sigma's
		dcheck=1,	# do an extra check using density and spline (2)
		prior='cauchy',	# prior can be cauchy or normal(0,1)
		plot=T,		# plot or no plot, that's the question
		wod=getwd(),	# working directory 
		bugsdir="c:/program files/winbugs14"		# WinBUGS directory: default is c:/program files/winbugs14
		)		
{
library(R2WinBUGS)
setwd(wod)

tvalue=c()
BF=c()

if (sample==1)
{

n=length(group1)
if (length(group2)==1)
{
X=group1
} else
{
X=group1-group2
}

X=X/sd(X)

data=list('X','n')

if (prior=='cauchy')
{
modelString=
"
model{delta~dnorm(0,lambdaDelta)
	lambdaDelta ~ dchisqr(1)
	sigma.~dnorm(0,sigmaL)
	sigmaL~dchisqr(1)
	sigma<-abs(sigma.)
	lambdaData<-pow(sigma,-2)
	mu<-delta*sigma
	for (i in 1:n){X[i] ~ dnorm(mu,lambdaData)}}
"
}

if (prior=='normal')
{
modelString=
"
model{delta~dnorm(0,1)
	sigma.~dnorm(0,sigmaL)
	sigmaL~dchisqr(1)
	sigma<-abs(sigma.)
	lambdaData<-pow(sigma,-2)
	mu<-delta*sigma
	for (i in 1:n){X[i] ~ dnorm(mu,lambdaData)}}
"
}

write(modelString, file='modelfile.bug')
inits=function()
	{
	list(delta=rnorm(1,0,1),sigma.=runif(1,0,5))
	}
parameters=c('delta','mu','sigma')
}

if ((sample==2)&&(sig==1))
{

n1=length(group1)
n2=length(group2)

group2=group2-mean(group1)
group2=group2/sd(group1)
group1=as.vector(scale(group1))

data=list('group1','group2','n1','n2')

if (prior=='cauchy')
{

modelString=
"
model{delta~dnorm(0,lambdaDelta)
	lambdaDelta ~ dchisqr(1)
	sigma.~dnorm(0,sigmaL)
	sigmaL~dchisqr(1)
	sigma<-abs(sigma.)
	var<-pow(sigma,2)
	mu~dnorm(0,muL)
	muL~dchisqr(1)
	lambdaData<- 1/var
	alpha<-delta*sigma	
	muData1<-mu+alpha*0.5	
	muData2<-mu-alpha*0.5	
	for (i in 1:n1){group1[i] ~ dnorm(muData1,lambdaData)}
	for (i in 1:n2){group2[i] ~ dnorm(muData2,lambdaData)}}
"

}

if (prior=='normal')
{

modelString=
"
model{delta~dnorm(0,1)
	sigma.~dnorm(0,sigmaL)
	sigmaL~dchisqr(1)
	sigma<-abs(sigma.)
	var<-pow(sigma,2)
	mu~dnorm(0,muL)
	muL~dchisqr(1)
	lambdaData<- 1/var
	alpha<-delta*sigma	
	muData1<-mu+alpha*0.5	
	muData2<-mu-alpha*0.5	
	for (i in 1:n1){group1[i] ~ dnorm(muData1,lambdaData)}
	for (i in 1:n2){group2[i] ~ dnorm(muData2,lambdaData)}}
"

}

write(modelString, file='modelfile.bug')

inits=function()
	{
	list(delta=rnorm(1,0,1),mu=rnorm(1,0,1),sigma.=runif(1,0,5))
	}
parameters=c('delta','mu','sigma','alpha')
}


if ((sample==2)&&(sig==2))
{

n1=length(group1)
n2=length(group2)

group2=group2-mean(group1)
group2=group2/sd(group1)
group1=as.vector(scale(group1))

data=list('group1','group2','n1','n2')
if (prior=='cauchy')
{
modelString=
"
model{delta~dnorm(0,lambdaDelta)
	lambdaDelta ~ dchisqr(1)
	sigma.1~dnorm(0,sigmaL1)
	sigmaL1~dchisqr(1)
	sigma1<-abs(sigma.1)
	sigma.2~dnorm(0,sigmaL2)
	sigmaL2~dchisqr(1)
	sigma2<-abs(sigma.2)
	var1<-pow(sigma1,2)
	var2<-pow(sigma2,2)
	mu~dnorm(0,muL)
	muL~dchisqr(1)
	lambdaData1<- 1/var1
	lambdaData2<- 1/var2
	alpha<-delta*( sqrt( ( var1*(n1-1) + var2*(n2-1) )/(n1+n2-2) ) )	
	muData1<-mu+alpha*0.5	
	muData2<-mu-alpha*0.5	
	for (i in 1:n1){group1[i] ~ dnorm(muData1,lambdaData1)}
	for (i in 1:n2){group2[i] ~ dnorm(muData2,lambdaData2)}}
"
}

if (prior=='normal')
{
modelString=
"
model{delta~dnorm(0,1)
	sigma.1~dnorm(0,sigmaL1)
	sigmaL1~dchisqr(1)
	sigma1<-abs(sigma.1)
	sigma.2~dnorm(0,sigmaL2)
	sigmaL2~dchisqr(1)
	sigma2<-abs(sigma.2)
	var1<-pow(sigma1,2)
	var2<-pow(sigma2,2)
	mu~dnorm(0,muL)
	muL~dchisqr(1)
	lambdaData1<- 1/var1
	lambdaData2<- 1/var2
	alpha<-delta*( sqrt( ( var1*(n1-1) + var2*(n2-1) )/(n1+n2-2) ) )	
	muData1<-mu+alpha*0.5	
	muData2<-mu-alpha*0.5	
	for (i in 1:n1){group1[i] ~ dnorm(muData1,lambdaData1)}
	for (i in 1:n2){group2[i] ~ dnorm(muData2,lambdaData2)}}
"
}

write(modelString, file='modelfile.bug')

inits=function()
	{
	list(delta=rnorm(1,0,1),mu=rnorm(1,0,1),sigma.1=runif(1,0,5),sigma.2=runif(1,0,5))
	}
parameters=c('delta','mu','sigma1','sigma2','alpha')
}

OST = bugs(data, inits, parameters,	model.file ="modelfile.bug",
	 			n.chains=chains, n.iter=iters, n.burnin=burns, n.thin=thins,
	 			bugs.directory=bugsdir,	codaPkg=F,debug=F)

dat1=matrix(,chains,3)
rn=c()
for (i in 1:chains)
{
ChC=paste('chain',i)
rn=append(rn,ChC)
}
rownames(dat1)=rn
colnames(dat1)=c('no OR','OR1','OR2')

if (prior=='cauchy')
{
	d0_prior1=dcauchy(0)
	d0_prior2=dcauchy(0)*2
}

if (prior=='normal')
{
	d0_prior1=dnorm(0,0,1)
	d0_prior2=dnorm(0,0,1)*2
}
for (i in 1:chains)
	{
	post.dat=OST$sims.array[,i,1]
	m=mean(post.dat)
	sd=sd(post.dat)
	d0_post=dnorm(0,m,sd)
	dat1[i,1]=d0_post/d0_prior1
	d0_post1=d0_post/pnorm(0,m,sd)
	dat1[i,2]=d0_post1/d0_prior2
	d0_post2=d0_post/(1-pnorm(0,m,sd))
	dat1[i,3]=d0_post2/d0_prior2
	}
if (dcheck!=1)
{

dat2=matrix(,chains,3)
rn=c()
for (i in 1:chains)
{
ChC=paste('chain',i)
rn=append(rn,ChC)
}
rownames(dat2)=rn
colnames(dat2)=c('no OR','OR1','OR2')

for (i in 1:chains)
	{
	post.dat=OST$sims.array[,i,1]
	d=density(post.dat)
	f=splinefun(d)
	d0_post=f(0)
	dat2[i,1]=d0_post/d0_prior1
	d0_post1=d0_post/integrate(f,min(d$x),0)$value
	dat2[i,2]=d0_post1/d0_prior2
	d0_post2=d0_post/(1-integrate(f,min(d$x),0)$value)
	dat2[i,3]=d0_post2/d0_prior2
	}
}

dat.SD=new.env()
dat.SD$Post.delta=OST$sims.array[,,1]
dat.SD$Post.mu=OST$sims.array[,,2]
if (length(parameters)==3){
	(dat.SD$Post.sigma=OST$sims.array[,,3])
	meanBF=mean(dat1[,1])
	PH0D=meanBF/(meanBF+1)
	PH1D=1-PH0D
	med=mean(X)
	dcount=0
	if (((iters-burns)/thins)>1000) (sa=1000)
	if (((iters-burns)/thins)<1000) (sa=(iters-burns)/thins)
	ind=sample(1:((iters-burns)/thins),sa)
	for (i in 1:sa)
	{
	mu.=OST$sims.array[ind[i],1,2]
	si.=OST$sims.array[ind[i],1,3]
	dd=rnorm(n,mu.,si.)
	delta=mean(dd)/sd(dd)
	if (delta*med>0) (dcount=dcount+1)
	}
	prep=PH0D*0.5 + PH1D*(dcount/sa)
}
if (length(parameters)==4){
	(dat.SD$Post.sigma=OST$sims.array[,,3])
	(dat.SD$Post.alpha=OST$sims.array[,,4])

	meanBF=mean(dat1[,1])
	PH0D=meanBF/(meanBF+1)
	PH1D=1-PH0D
	med=mean(group1)-mean(group2)
	dcount=0
	if (((iters-burns)/thins)>10000) (sa=10000)
	if (((iters-burns)/thins)<10000) (sa=(iters-burns)/thins)
	ind=sample(1:((iters-burns)/thins),sa)
	for (i in 1:sa)
	{
	mu.=OST$sims.array[ind[i],1,2]
	si.=OST$sims.array[ind[i],1,3]
	al.=OST$sims.array[ind[i],1,4]
	dd1=rnorm(n1,mu.+al./2,si.)
	dd2=rnorm(n2,mu.-al./2,si.)
	delta=(mean(dd1)-mean(dd2))/sqrt( ( (sd(dd1))^2*(n1-1) + (sd(dd2))^2*(n2-1) )/(n1+n2-2) )
	if (delta*med>0) (dcount=dcount+1)
	}
	prep=PH0D*0.5 + PH1D*(dcount/sa)
}
if (length(parameters)==5)
	{
	(dat.SD$Post.sigma1=OST$sims.array[,,3])
	(dat.SD$Post.sigma2=OST$sims.array[,,4])
	(dat.SD$Post.alpha=OST$sims.array[,,5]) 
	meanBF=mean(dat1[,1])
	PH0D=meanBF/(meanBF+1)
	PH1D=1-PH0D
	med=mean(group1)-mean(group2)
	dcount=0
	if (((iters-burns)/thins)>10000) (sa=10000)
	if (((iters-burns)/thins)<10000) (sa=(iters-burns)/thins)
	ind=sample(1:((iters-burns)/thins),sa)
	for (i in 1:sa)
	{
	mu.=OST$sims.array[ind[i],1,2]
	si1.=OST$sims.array[ind[i],1,3]
	si2.=OST$sims.array[ind[i],1,3]
	al.=OST$sims.array[ind[i],1,5]
	dd1=rnorm(n1,mu.+al./2,si1.)
	dd2=rnorm(n2,mu.-al./2,si2.)
	delta=(mean(dd1)-mean(dd2))/sqrt( ( (sd(dd1))^2*(n1-1) + (sd(dd2))^2*(n2-1) )/(n1+n2-2) )
	if (delta*med>0) (dcount=dcount+1)
	}
	prep=PH0D*0.5 + PH1D*(dcount/sa)
}

dat.SD$BF=dat1
dat.SD$summary=OST$summary
if (dcheck==2)
{
	dat.SD$BF_check=dat2
}


dat.SD=as.list(dat.SD)


cat('Results of a', sample, 'sample SD Bayesian t-test.','\n')
if (chains>1)
{ 
cat('Rhat, calculated from',chains,'chains is', round(OST$summary[1,8],3),'(should be <1.1 ).','\n')
cat('The mean Bayes factor, calculated from',chains,'chains=',apply(dat1,2,mean),'.','\n')
cat('Bayesian p-rep =',prep,'.','\n')
}
cat('All Bayes Factors are displayed below:','\n','\n')
print(dat1)
if (dcheck==2)
{
cat('\n')
cat('All Bayes Factors are displayed below (using splinefunction):','\n','\n')
print(dat2)
}


if (plot==T)
{

d=density(dat.SD$Post.delta)
f=splinefun(d)
d0_post=f(0)

dev.new(height=1,width=16)
par(cex.main = 2, mar = c(5, 5, 4, 1) + 0.1, mgp = c(3.5, 1, 0),
cex.lab = 2, font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
layout(matrix(1:3,1,3))

plot(d,lwd=2,lty=2,xlim=c(-5,5),ylim=c(0,3),ylab='Density',xlab=expression(paste('effect size ',delta)),axes=F,main=substitute(paste(BF["01"]," = ", k),list(k=round(mean(dat1[,1]),2))))
if (prior=='cauchy')
{
plot(function(c) dcauchy(c,0,1),-5,5,add=T,lwd=2)
}
if (prior=='normal')
{
plot(function(c) dnorm(c,0,1),-5,5,add=T,lwd=2)
}
points(0,d0_prior1,pch=20,cex=3)
points(0,d0_post,pch=20,cex=3)
text(-5,(4/10)*3,expression(paste(H[0],':',delta,'=0')),cex=2,pos=4)
text(-5,(3/10)*3,expression(paste(H[1],':',delta!=0)),cex=2,pos=4)
axis(1,labels=T,at=c(-5,0,5),cex.axis=2,cex.lab=2)
axis(2,labels=F,tcl=0)

plot(function(c) f(c)*(1/integrate(f,min(d$x),0)$value), d$x[1],0,lwd=2,lty=2,xlim=c(-5,0),ylim=c(0,max(spline(d)$y)*(1/integrate(f,min(d$x),0)$value)),ylab='',xlab=expression(paste('effect size ',delta)),main=substitute(paste(BF["01"]," = ", k),list(k=round(mean(dat1[,2]),2))),axes=F)
lines(c(0,0),c(0,f(0)*(1/integrate(f,min(d$x),0)$value)),lwd=2,lty=2)

if (prior=='cauchy')
{
plot(function(c) dcauchy(c,0,1)*2,-5,0,add=T,lwd=2)
lines(c(0,0),c(0,dcauchy(0)*2),lwd=2)
}
if (prior=='normal')
{
plot(function(c) dnorm(c,0,1)*2,-5,0,add=T,lwd=2)
lines(c(0,0),c(0,dnorm(0)*2),lwd=2)
}
points(0,d0_prior2,pch=20,cex=3)
points(0,d0_post*(1/integrate(f,min(d$x),0)$value),pch=20,cex=3)
text(-5,(4/10)*max(spline(d)$y)*(1/integrate(f,min(d$x),0)$value),expression(paste(H[0],':',delta,'=0')),cex=2,pos=4)
text(-5,(3/10)*max(spline(d)$y)*(1/integrate(f,min(d$x),0)$value),expression(paste(H[1],':',delta<0)),cex=2,pos=4)
axis(1,labels=T,at=c(-5,-2.5,0),cex.axis=2,cex.lab=2)
axis(2,labels=F,tcl=0)
legend('topleft',c('posterior','prior'),lty=2:1,lwd=2,cex=2,bty='n')

plot(function(c) f(c)*(1/integrate(f,0,max(d$x))$value), d$x[length(d$x)],0,lwd=2,lty=2,xlim=c(0,5),ylim=c(0,max(spline(d)$y)*(1/integrate(f,0,max(d$x))$value)),ylab='',xlab=expression(paste('effect size ',delta)),main=substitute(paste(BF["01"]," = ", k),list(k=round(mean(dat1[,3]),2))),axes=F)
lines(c(0,0),c(0,f(0)*(1/integrate(f,0,max(d$x))$value)),lwd=2,lty=2)
if (prior=='cauchy')
{
plot(function(c) dcauchy(c,0,1)*2,0,5,add=T,lwd=2)
lines(c(0,0),c(0,dcauchy(0)*2),lwd=2)
}
if (prior=='normal')
{
plot(function(c) dnorm(c,0,1)*2,0,5,add=T,lwd=2)
lines(c(0,0),c(0,dnorm(0)*2),lwd=2)

}
points(0,d0_prior2,pch=20,cex=3)
points(0,d0_post*(1/integrate(f,0,max(d$x))$value),pch=20,cex=3)
text(3,(4/10)*max(spline(d)$y)*(1/integrate(f,0,max(d$x))$value),expression(paste(H[0],':',delta,'=0')),cex=2,pos=4)
text(3,(3/10)*max(spline(d)$y)*(1/integrate(f,0,max(d$x))$value),expression(paste(H[1],':',delta>0)),cex=2,pos=4)
axis(1,labels=T,at=c(0,2.5,5),cex.axis=2,cex.lab=2)
axis(2,labels=F,tcl=0)

}
return(dat.SD)
}



