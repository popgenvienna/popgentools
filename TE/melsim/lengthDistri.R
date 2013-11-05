d<-read.table("/Volumes/Volume_4/analysis/DsimTE/refgen/TEannotation/stat/lengthdistri.R")
ids=c("ssr","g2","g1","g05")

par(mfrow=c(2,2))

avleng<-function(df){
  c<-sum(df$V4)
  ls<-sum(df$V4*df$V3)
  al<-ls/c
  return(al)
}

histable<-function(df){
  v<-c()
  for(i in 1:nrow(df))
  {
    cur<-df[i,]
    ele<-rep(cur$V3,cur$V4)
    v<-c(v,ele)
  }
  return(v)
}

for(id in ids)
{
  ss<-d[d$V2==id,]
  le=sum(ss$V4)
  al<-avleng(ss)
  al<-as.integer(al*100)
  al<-as.numeric(al)/100
  print(al)
  header=paste(id," count=",le," av.leng.=",al,sep="")
  

  plot(ss$V3,ss$V4,type="l",main=header,log="x",xlim=c(1,15000),ylim=c(1,350),xlab="TE length",ylab="count")
  lines(c(50,50),c(0,350),col="red")
  
  #histdat<-histable(ss)
  #print(histdat)
  #bre<-seq(0,30000,100)
  #bre<-seq(2,5,0.1)
  #bre<-c(0,10^bre)
  #print(bre)
  #hist.data<-hist(histdat, plot=F,breaks=bre)
  #hist.data$counts = log10(hist.data$counts)
  #plot(hist.data,xlim=c(0,3000),ylim=c(0,5))
  #hist(histable(ss),breaks=bre,xlim=c(0,3000))
  
}

