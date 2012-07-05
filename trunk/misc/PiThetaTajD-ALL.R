t=read.table("PI")
t$V5[t$V5=="na"]=NA
t$V5=as.numeric(as.character(t$V5))

u=read.table("THETA")
u$V5[u$V5=="na"]=NA
u$V5=as.numeric(as.character(u$V5))

w=read.table("TAJD")
w$V5[w$V5=="na"]=NA
w$V5=as.numeric(as.character(w$V5))

chroms<-chromosomes

graphs=sum(length(sapply(chromosomes,length)))

pdf("PLOTS-PI.pdf",width=15)
par(mfrow=c(3,(graphs+3)/3))
maxval=max(t$V5,na.rm=T)
minval=min(t$V5,na.rm=T)
xmaxtotal=max(t$V2,na.rm=T)
for(chro in chroms)
{
	xmax=0
	xmaxindiv=max(t[t$V1==chro,2])
	ifelse(scale==1,xmax<-xmaxtotal,xmax<-xmaxindiv)
	plot(t[t$V1==chro,2],t[t$V1==chro,5],t="l",xlab="position",ylab="pi",main=chro,xlim=c(0,xmax),ylim=c(minval,maxval),col="black")
}
dev.off()


pdf("PLOTS-THETA.pdf",width=15)
par(mfrow=c(3,(graphs+3)/3))
maxval=max(u$V5,na.rm=T)
minval=min(u$V5,na.rm=T)
xmaxtotal=max(u$V2,na.rm=T)
for(chro in chroms)
{
	xmax=0
	xmaxindiv=max(u[u$V1==chro,2])
	ifelse(scale==1,xmax<-xmaxtotal,xmax<-xmaxindiv)
	plot(u[u$V1==chro,2],u[u$V1==chro,5],t="l",xlab="position",ylab="theta",main=chro,xlim=c(0,xmax),ylim=c(minval,maxval),col="black")
}
dev.off()

pdf("PLOTS-TAJD.pdf",width=15)
par(mfrow=c(3,(graphs+3)/3))
maxval=max(w$V5,na.rm=T)
minval=min(w$V5,na.rm=T)
xmaxtotal=max(w$V2,na.rm=T)
for(chro in chroms)
{
	xmax=0
	xmaxindiv=max(w[w$V1==chro,2])
	ifelse(scale==1,xmax<-xmaxtotal,xmax<-xmaxindiv)
	plot(w[w$V1==chro,2],w[w$V1==chro,5],t="l",xlab="position",ylab="taj d",main=chro,xlim=c(0,xmax),ylim=c(minval,maxval),col="black")
}
dev.off()

pdf("PLOTS-PI_&_THETA.pdf",width=15)
par(mfrow=c(3,(graphs+3)/3))
tmax=max(t$V5,na.rm=T)
umax=max(u$V5,na.rm=T)
maxvaltu=max(tmax,umax)
tmin=min(t$V5,na.rm=T)
umin=min(u$V5,na.rm=T)
minvaltu=min(tmin,umin)
txmaxtotal=max(t$V2,na.rm=T)
uxmaxtotal=max(u$V2,na.rm=T)
error<-("error")
ifelse(txmaxtotal==uxmaxtotal,xmaxtotal<-max(txmaxtotal),error)
for(chro in chroms)
{
	xmax=0
	xmaxindiv=max(t[t$V1==chro,2])
	ifelse(scale==1,xmax<-xmaxtotal,xmax<-xmaxindiv)
	plot(t[t$V1==chro,2],t[t$V1==chro,5],t="l",xlab="position",ylab="pi and theta",main=chro,xlim=c(0,xmax),ylim=c(minvaltu,maxvaltu),col="black")
	lines(u[u$V1==chro,2],u[u$V1==chro,5],t="l",xlab="position",ylab="pi and theta",main=chro,xlim=c(0,xmax),ylim=c(minvaltu,maxvaltu),col="red")
}
dev.off()
