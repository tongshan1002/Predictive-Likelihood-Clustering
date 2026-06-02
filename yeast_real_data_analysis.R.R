rm(list = ls())

genedata <- read.table("C:/research/GENE/genemult.txt", header=F)

lgenemult<-log(genedata)

library(cluster)
library(splines)
library(fossil)




meanlgenemult<-apply(lgenemult,1,mean)
meanlgenemat<-matrix(meanlgenemult,nrow=78,ncol=18,byrow=F)
lgenemultc<-lgenemult-meanlgenemat


endtime<-7*18/60
starttime<- 7/60

timevec<-seq(starttime,endtime, 7/60)

tranobscurves<-t(lgenemultc)

VDPdata <- cbind(timevec,tranobscurves)



btgroup<-c(rep(1,13),rep(2,39),rep(3,8),rep(4,7),rep(5,11))
valuelikeli<-c()

VDPmat <- VDPdata
VDPmat.resp <- VDPmat[,-1]
size=5
clust.vec=size
sam<-nrow(VDPmat)
col<-ncol(VDPmat.resp)
k<-max(clust.vec)
group<-sample(1:k,col,replace=TRUE)
count=0
trend<-c()
value<-c()
bestgroup<-c()
rand<-c()
grouprand<-c()

beta0<-c()
beta1<-c()
beta2<-c()
beta3<-c()
beta4<-c()
beta5<-c()
beta6<-c()
beta7<-c()



count.odd<-1
matrix50<-matrix(rep(0,ncol(VDPmat.resp)*50),nrow=50)
value50<-c(rep(0,50))
valuetop20<-c(rep(0,20))


#Generate fourier T

x1=VDPmat[,1]
T1=matrix(c(rep(1,sam)),ncol=1)
T2=matrix(x1,ncol=1)
T3=matrix(x1^2,ncol=1)
T4=matrix(x1^3,ncol=1)
T5=matrix(cos(2*pi*x1),ncol=1)
T6=matrix(sin(2*pi*x1),ncol=1)
T7=matrix(cos(4*pi*x1),ncol=1)
T8=matrix(sin(4*pi*x1),ncol=1)
T=cbind(T1,T2,T3,T4,T5,T6,T7,T8)




predictlikelihood<-function(gr){
  
  n<-vector("list",1)
  u<-vector("list",1)
  m<-c()
  
  #number of each group
  for (i in 1:length(unique(gr))){
    m[i]<-as.data.frame(table(gr))[i,2]
    n[[i]]<-vector()
  }
  
  
  #put data into k groups
  
  for (i in 1:length(unique(gr))){
    for (j in 1:ncol(VDPmat.resp)) {
      if (gr[j]==i) {
        n[[i]]=c(n[[i]],j)
      }
    }}
  
  for (i in 1:length(unique(gr))){
    u[[i]]<-VDPmat.resp[,n[[i]]]
  }
  
  for (i in 1:length(unique(gr))){
    if (m[i]==1){u[[i]]=matrix(u[[i]])}
  }
  
  
  # get sigma_ij for k groups
  sigma<-c()
  sigma_n<-c()
  
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      linear.model<-lm(u[[i]][,j]~x1+I(x1^2)+I(x1^3)+I(cos(2*pi*x1))+I(sin(2*pi*x1))+I(cos(4*pi*x1))+I(sin(4*pi*x1)))
      sig<-summary(linear.model)$sigma^2
      sigma<-c(sigma,sig)
      
    }
    sigma_n[i]<-mean(sigma)
    sigma<-c()
    
  }
  
  
  
  A_n<-vector("list",1)
  S_n<-vector("list",1)
  I=diag(sam)
  
  #Generate S_n
  for (i in 1:length(unique(gr))){
    A_n[[i]]=(-1/(m[i]-1))*T%*%(solve((t(T)%*%T)))%*%t(T)
    
    #S_n
    S_n[[i]]=matrix(rep(0,(m[i]*sam)*(m[i]*sam)),nrow=(m[i]*sam))
  }
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      for(k in 1:m[i]){
        S_n[[i]][((j-1)*sam+1):(j*sam),((k-1)*sam+1):(k*sam)]<-A_n[[i]]  
      }
    }
  }
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      S_n[[i]][((j-1)*sam+1):(j*sam),((j-1)*sam+1):(j*sam)]<-I
    }
  }
  
  #get M_n
  mvector<-vector("list",1)
  M<-vector("list",1)
  
  for (i in 1:length(unique(gr))){
    mvector[[i]]<-vector()
    M[[i]]<-vector()
  }
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      mvector[[i]]<-c(mvector[[i]],u[[i]][,j])
    }
    mvector[[i]]<-matrix(mvector[[i]])
    M[[i]]=t(mvector[[i]])%*%S_n[[i]]%*%mvector[[i]]
  }
  
  
  #get sigmasse and old c
  sigmasse<-vector("list",1)
  Nu<-c()
  logGpart<-c()
  
  for (i in 1:length(unique(gr))){
    if (m[i]==1){sigmasse[[i]]=sigma_n[i]}
    else{
      sigmasse[[i]]<- 1/m[i]*sigma_n[i]+1/(m[i]*sam-length(unique(gr)))*((m[i]-1)/m[i]*M[[i]])
    }
    Nu[i]=(m[i]*sam-length(unique(gr)))/2
    logGpart[i]= lfactorial(m[i])-0.5* log(m[i]) - (Nu[i]-1)*log(sigmasse[[i]])+lgamma(Nu[i])-Nu[i]*log(Nu[i])
  }
  
  oldc<-0
  for (i in 1:length(unique(gr))){
    oldc<-sum(oldc,logGpart[i])
  }
  
  return(oldc)
}



predictlikelihoodsilhouette<-function(gr){
  
  n<-vector("list",1)
  u<-vector("list",1)
  m<-c()
  
  #number of each group
  for (i in 1:length(unique(gr))){
    m[i]<-as.data.frame(table(gr))[i,2]
    n[[i]]<-vector()
  }
  
  
  #put data into k groups
  
  for (i in 1:length(unique(gr))){
    for (j in 1:ncol(VDPmat.resp)) {
      if (gr[j]==i) {
        n[[i]]=c(n[[i]],j)
      }
    }}
  
  for (i in 1:length(unique(gr))){
    u[[i]]<-VDPmat.resp[,n[[i]]]
  }
  
  for (i in 1:length(unique(gr))){
    if (m[i]==1){u[[i]]=matrix(u[[i]])}
  }
  
  
  # get sigma_ij for k groups
  sigma<-c()
  sigma_n<-c()
  
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      linear.model<-lm(u[[i]][,j]~x1+I(x1^2)+I(x1^3)+I(cos(2*pi*x1))+I(sin(2*pi*x1))+I(cos(4*pi*x1))+I(sin(4*pi*x1)))
      sig<-summary(linear.model)$sigma^2
      b0<-summary(linear.model)$coefficients[1]
      b1<-summary(linear.model)$coefficients[2]
      b2<-summary(linear.model)$coefficients[3]
      b3<-summary(linear.model)$coefficients[4]
      b4<-summary(linear.model)$coefficients[5]
      b5<-summary(linear.model)$coefficients[6]
      b6<-summary(linear.model)$coefficients[7]
      b7<-summary(linear.model)$coefficients[8]
      
      sigma<-c(sigma,sig)
      beta0<-c(beta0,b0)
      beta1<-c(beta1,b1)
      beta2<-c(beta2,b2)
      beta3<-c(beta3,b3)
      beta4<-c(beta4,b4)
      beta5<-c(beta5,b5)
      beta6<-c(beta6,b6)
      beta7<-c(beta7,b7)
      
      
      
    }
    sigma_n[i]<-mean(sigma)
    sigma<-c()
    
  }
  
  
  A_n<-vector("list",1)
  S_n<-vector("list",1)
  I=diag(sam)
  
  #Generate S_n
  for (i in 1:length(unique(gr))){
    A_n[[i]]=(-1/(m[i]-1))*T%*%(solve((t(T)%*%T)))%*%t(T)
    
    #S_n
    S_n[[i]]=matrix(rep(0,(m[i]*sam)*(m[i]*sam)),nrow=(m[i]*sam))
  }
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      for(k in 1:m[i]){
        S_n[[i]][((j-1)*sam+1):(j*sam),((k-1)*sam+1):(k*sam)]<-A_n[[i]]  
      }
    }
  }
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      S_n[[i]][((j-1)*sam+1):(j*sam),((j-1)*sam+1):(j*sam)]<-I
    }
  }
  
  #get M_n
  mvector<-vector("list",1)
  M<-vector("list",1)
  
  for (i in 1:length(unique(gr))){
    mvector[[i]]<-vector()
    M[[i]]<-vector()
  }
  
  for (i in 1:length(unique(gr))){
    for (j in 1:m[i]){
      mvector[[i]]<-c(mvector[[i]],u[[i]][,j])
    }
    mvector[[i]]<-matrix(mvector[[i]])
    M[[i]]=t(mvector[[i]])%*%S_n[[i]]%*%mvector[[i]]
  }
  
  
  #get sigmasse and old c
  sigmasse<-vector("list",1)
  Nu<-c()
  logGpart<-c()
  
  for (i in 1:length(unique(gr))){
    if (m[i]==1){sigmasse[[i]]=sigma_n[i]}
    else{
      sigmasse[[i]]<- 1/m[i]*sigma_n[i]+1/(m[i]*sam-length(unique(gr)))*((m[i]-1)/m[i]*M[[i]])
    }
    Nu[i]=(m[i]*sam-length(unique(gr)))/2
    logGpart[i]= lfactorial(m[i])-0.5* log(m[i]) - (Nu[i]-1)*log(sigmasse[[i]])+lgamma(Nu[i])-Nu[i]*log(Nu[i])
  }
  
  oldc<-0
  for (i in 1:length(unique(gr))){
    oldc<-sum(oldc,logGpart[i])
  }
  
  newlist <- list(oldc,beta0,beta1,beta2,beta3,beta4,beta5,beta6,beta7)
  
  return(newlist)
  
}




for (Tm in 1:100000){
  
  if ((Tm/3)==round(Tm/3)){count.odd<-abs(1-count.odd)}
  
  if (count.odd==1){
    
    #get old c
    oldvalue<- predictlikelihood(group)
    
    #find new partition c+1
    delta<-0.4^(log(col))
    groupnew<-c()
    if (k>2){alpha<-(k-1)/(k*(k-1)-1)} else {alpha=0.5}
    p=0
    
    for (i in 1:ncol(VDPmat.resp)){
      random<-runif(1)
      if (random<=1 & (1-(1-(k-1)*alpha)*delta)<random){
        groupnew[i]<-k+1
        p=1
      }
      if (p==0){    
        if (random < (1-delta)) {groupnew[i]=group[i]}
        for (j in (1:(k-1))){
          if ( (random >= (1-delta)+(j-1)*alpha*delta) & (random < (1-delta)+j*alpha*delta))
          { groupnew[i]=group[i]+1}
        }
        if (groupnew[i]==k+1){groupnew[i]=1}
      }
      p=0    
    }
    
    if (max(groupnew)!=length(unique(groupnew))){
      for (i in 1: length(unique(groupnew))){
        groupnew[groupnew %in% unique(groupnew)[i]]<-seq_along(unique(groupnew))[i]
      }  
    }
    
    #get groupnew c
    
    Tp=100/log(1+Tm)
    newvalue<- predictlikelihood(groupnew)
    
    probab<-min(exp((newvalue-oldvalue)/Tp),1)
    
    
  } else{
    beta0<-c()
    beta1<-c()
    beta2<-c()
    beta3<-c()
    beta4<-c()
    beta5<-c()
    beta6<-c()
    beta7<-c()
    
    #get old c
    
    oldvalue<- predictlikelihoodsilhouette(group)[[1]]
    
    if (length(unique(groupnew))==1 | length(unique(group))==1){groupnew=group
    count.odd<-abs(1-count.odd)}
    else{
      
      
      
      #find new partition c+1 using silhouette
      beta0<-predictlikelihoodsilhouette(group)[[2]]
      beta1<-predictlikelihoodsilhouette(group)[[3]]
      beta2<-predictlikelihoodsilhouette(group)[[4]]
      beta3<-predictlikelihoodsilhouette(group)[[5]]
      beta4<-predictlikelihoodsilhouette(group)[[6]]
      beta5<-predictlikelihoodsilhouette(group)[[7]]
      beta6<-predictlikelihoodsilhouette(group)[[8]]
      beta7<-predictlikelihoodsilhouette(group)[[9]]
      
      
      
      clustering.vec<-group
      
      matrix.of.betahats<-cbind(beta0,beta1,beta2,beta3,beta4,beta5,beta6,beta7)
      
      my.D<-dist(matrix.of.betahats)
      
      a<-silhouette(clustering.vec,dist=my.D)
      
      g1<-c()
      g2<-c()
      g3<-c()
      
      for (i in 1:ncol(VDPmat.resp)){
        g1<-c(g1,a[i,1])
        g2<-c(g2,a[i,2])
        g3<-c(g3,a[i,3])
      }
      
      g1<-matrix(g1)
      g2<-matrix(g2)
      g3<-matrix(g3)
      
      
      g1<-t(g1)
      g2<-t(g2)
      g3<-t(g3)
      
      ming<-min(g3)
      maxg<-max(g3)
      
      if (maxg==0){maxg=10000}
      
      
      probsil<-c()
      groupnew<-c()
      for (i in 1:col){
        probsil[i]<- max((g3[i]-ming)/(maxg-ming),g3[i]/maxg)
        uc<-runif(1,0,1)
        if (probsil[i]>=uc){groupnew[i]<-g1[i]}
        else {groupnew[i]<-g2[i]}
        
      }
      
      
      
      if (max(groupnew)!=length(unique(groupnew))){
        for (i in 1: length(unique(groupnew))){
          groupnew[groupnew %in% unique(groupnew)[i]]<-seq_along(unique(groupnew))[i]
        }  
      }
    }
    
    # get groupnew value
    Tp=100/log(1+Tm)
    
    newvalue<- predictlikelihoodsilhouette(groupnew)[[1]]
    probab<-min(exp((newvalue-oldvalue)/Tp),1)
  }
  
  if (all(group == groupnew)){count=count-1
  Tm=Tm-1}
  
  randin<-rand.index(btgroup,groupnew)
  rand<-c(rand,randin)
  
  
  u1<-runif(1,0,1)
  if((u1<probab)){print("accept")
    group=groupnew
    count=count+1
    value<-c(value,newvalue)
  }else{
    print("reject")
    group=group
    value<-c(value,value[length(value)])
  }
  
  
  
  if (Tm == 0) {
    rom <- 1
  } else {
    rom <- Tm %% 50
  }
  
  matrix50[rom,]<-groupnew
  value50[rom]<-newvalue
  
  if (rom==0) {
    valuetop20<-rev(tail(sort(value50),20))
    sum20<-sum(valuetop20)
    provalue20<- valuetop20/sum20
    provalue20<-c(provalue20)
    choose<-sample(1:20,size=1,prob=provalue20)
    p<-which(value50==valuetop20[choose])
    if (length(p)>1){p<-p[1]}
    group=matrix50[p,]}
  
  if (newvalue==max(value)){bestgroup<-groupnew}
  
  
  
  
  
  
  
  
  
  valuelikeli<-c(valuelikeli,value)
  
  acceptrate=count/Tm
  print(acceptrate)
  print(Tm)
  
}

write.csv(rand,"C:/research/GENE/for1130rand50tho.csv",row.names = FALSE)
write.csv(bestgroup,"C:/research/GENE/for1130bestgroup50tho.csv",row.names = FALSE)
write.csv(value,"C:/research/GENE/for1130value50tho.csv",row.names = FALSE)

