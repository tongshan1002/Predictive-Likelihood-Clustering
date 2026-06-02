library(cluster)
library(fossil)
library(splines)

mylambda<-0.4
mydf<-2
simul<-100  # Increase "simul" if you want to create many simulated functional data sets
randtot<-c()

# initializing counters:
totpropor.ward<-0; totpropor.single<-0; totpropor.complete<-0; totpropor.average<-0;
totpropor.median<-0; totpropor.centroid<-0;

vecnorm<-function(x){
  mynorm<-sum(x^2)^(1/2)
  return(mynorm)
}

for (iter in 1:simul)
{
  times<-50
  endtime<-5
  dt<-endtime/times
  sigma<-0.5
  drift<-1
  
  myk<-5            # CHANGED: now 5 clusters
  objects<-40       # still 40 total curves
  
  normalvec<-matrix(rnorm(objects*times, sd=sigma), nrow=objects, ncol=times)
  ouvalue<-matrix(0, nrow=objects, ncol=times)
  
  # Defining the cluster sizes and the true clustering structure:
  # CHANGED: 5 clusters, each has 8 curves => cumulative sizes:
  clussizes<-c(8, 16, 24, 32, objects)  # objects=40
  
  trueclust<-c(
    rep(1, clussizes[1]),
    rep(2, clussizes[2]-clussizes[1]),
    rep(3, clussizes[3]-clussizes[2]),
    rep(4, clussizes[4]-clussizes[3]),
    rep(5, clussizes[5]-clussizes[4])
  )
  
  timevec<-seq(0, endtime, length=times)   # tgrid
  
  #################################################################
  # loop to create approximation of Ornstein-Uhlenbeck process
  # (same OU structure, now for 5 groups)
  
  # cluster 1
  for (timeint in 2:times)
  {
    ouvalue[1:clussizes[1], timeint] <-
      -drift*ouvalue[1:clussizes[1], timeint-1]*dt +
      sigma*normalvec[1:clussizes[1], timeint]*sqrt(dt)
  }
  
  # cluster 2
  for (timeint in 2:times)
  {
    ouvalue[(clussizes[1]+1):clussizes[2], timeint] <-
      -drift*ouvalue[(clussizes[1]+1):clussizes[2], timeint-1]*dt +
      sigma*normalvec[(clussizes[1]+1):clussizes[2], timeint]*sqrt(dt)
  }
  
  # cluster 3
  for (timeint in 2:times)
  {
    ouvalue[(clussizes[2]+1):clussizes[3], timeint] <-
      -drift*ouvalue[(clussizes[2]+1):clussizes[3], timeint-1]*dt +
      sigma*normalvec[(clussizes[2]+1):clussizes[3], timeint]*sqrt(dt)
  }
  
  # cluster 4
  for (timeint in 2:times)
  {
    ouvalue[(clussizes[3]+1):clussizes[4], timeint] <-
      -drift*ouvalue[(clussizes[3]+1):clussizes[4], timeint-1]*dt +
      sigma*normalvec[(clussizes[3]+1):clussizes[4], timeint]*sqrt(dt)
  }
  
  # cluster 5 (NEW)
  for (timeint in 2:times)
  {
    ouvalue[(clussizes[4]+1):clussizes[5], timeint] <-
      -drift*ouvalue[(clussizes[4]+1):clussizes[5], timeint-1]*dt +
      sigma*normalvec[(clussizes[4]+1):clussizes[5], timeint]*sqrt(dt)
  }
  
  #################################################################
  # Simulations of Signal curves (true underlying curves)
  
  curves<-matrix(0, nrow=objects, ncol=times)
  
  # cluster 1: 1-8
  for (rowi in 1:clussizes[1])
  {
    curves[rowi,] <- -sin(timevec-1)*log(timevec+0.5)
  }
  
  # cluster 2: 9-16
  for (rowi in (clussizes[1]+1):clussizes[2])
  {
    curves[rowi,] <- cos(timevec)*log(timevec+0.5)
  }
  
  # cluster 3: 17-24
  for (rowi in (clussizes[2]+1):clussizes[3])
  {
    curves[rowi,] <- -0.25 - 0.1*cos(0.5*(timevec-1))*sqrt(5*sqrt(timevec)+0.5)*(timevec^1.5)
  }
  
  # cluster 4: 25-32
  for (rowi in (clussizes[3]+1):clussizes[4])
  {
    curves[rowi,] <- -1 + 0.3*timevec
  }
  
  # cluster 5: 33-40 (NEW signal)
  mu5 <- 0.2*(timevec - 2.5)^2 - 1
  for (rowi in (clussizes[4]+1):clussizes[5])
  {
    curves[rowi,] <- mu5
  }
  
  # observed curves
  obscurves <- curves + ouvalue
  tranobscurves<-t(obscurves)
  
  VDPdata<-cbind(timevec,tranobscurves)
  VDPmat <- VDPdata
  VDPmat.resp <- VDPmat[,-1]
  size=3
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
  bestrand<-c()
  btgroup<-c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8))
  
  
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
  #B-splines
  #Knots at equally spaced quantiles of time point and get B-spline matrix
  
  B.matrix<-function(TP, K, r){
    require(splines)
    u <- quantile(x=TP, probs=(0:K)/K)
    u.all <- c(rep(min(u),r),u,rep(max(u),r))
    B <- splineDesign(knots=u.all, x=TP, ord = r+1)
    
    return(B)
  }
  
  
  #tp time point
  #knots
  #r =degree
  
  #knots+degree=number of column=knots+order-1
  #order=degree+1
  TP<-timevec
  T=B.matrix(TP,5,3)
  
  
  
  
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
        linear.model<-lm(u[[i]][,j]~x1+I(x1^2)+I(x1^3))
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
        linear.model<-lm(u[[i]][,j]~x1+I(x1^2)+I(x1^3))
        sig<-summary(linear.model)$sigma^2
        b0<-summary(linear.model)$coefficients[1]
        b1<-summary(linear.model)$coefficients[2]
        b2<-summary(linear.model)$coefficients[3]
        b3<-summary(linear.model)$coefficients[4]
        
        sigma<-c(sigma,sig)
        beta0<-c(beta0,b0)
        beta1<-c(beta1,b1)
        beta2<-c(beta2,b2)
        beta3<-c(beta3,b3)
        
        
        
        
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
    
    newlist <- list(oldc,beta0,beta1,beta2,beta3)
    
    return(newlist)
    
  }
  
  
  
  
  for (Tm in 1:10000){
    
    if ((Tm/3)==round(Tm/3)){count.odd<-abs(1-count.odd)}
    
    if (count.odd==1){
      
      #get old c
      oldvalue<- predictlikelihood(group)
      
      #find new partition c+1
      delta<-0.4^(log(col))
      groupnew<-c()
      if (k>2){alpha<-(k-1)/(k*(k-1)-1)} else {alpha=0.5}# no run
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
      
      Tp=100/log(log(1+Tm))
      newvalue<- predictlikelihood(groupnew)
      
      probab<-min(exp((newvalue-oldvalue)/Tp),1)
      
      
    } else{
      beta0<-c()
      beta1<-c()
      beta2<-c()
      beta3<-c()
      
      
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
        
        
        
        clustering.vec<-group
        
        matrix.of.betahats<-cbind(beta0,beta1,beta2,beta3)
        
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
      Tp=100/log(log(1+Tm))
      
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
    
    
    
    rom<-Tm%%50
    matrix50[rom,]<-groupnew
    value50[rom]<-newvalue
    
    if (rom==49) {
      valuetop20<-rev(tail(sort(value50),20))
      sum20<-sum(valuetop20)
      provalue20<- valuetop20/sum20
      provalue20<-c(provalue20)
      choose<-sample(1:20,size=1,prob=provalue20)
      p<-which(value50==valuetop20[choose])
      if (length(p)>1){p<-p[1]}
      group=matrix50[p,]}
    
    if (newvalue==max(value)){bestgroup<-groupnew}
    
    acceptrate=count/Tm
    print(iter)
    print(acceptrate)
    
  }
  
  randtot<-c(randtot,max(rand))
  
}


write.csv(randtot, "C:/research/2026/randtotbsplinebasis40curvesigma5k5.csv", row.names = FALSE)

