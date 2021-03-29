

###########################################
## sequential advantage selection
###########################################
I_SS<-function(X,Y,trt,g_true,step,rho){
  #some constants
  n<-length(Y);p<-dim(X)[2]
  
  #store the treatment regime
  #critical<-rep(0,step);direction<-rep(0,step);err.rate<-rep(0,step)
  #V_outcome: use the true model; V_est: use the A-learning method to estimate
  V_outcome<-rep(0,step+1);V_est<-matrix(0,nrow=n,ncol=step+1)
  a_opt<-matrix(0,nrow=n,ncol=(step+1));S<-rep(0,step+1)
  regime_coef<-matrix(0,nrow=(step+1),ncol=(p+1))
  err.rate<-rep(0,step+1)
  
  reg0<-lm(Y~trt)
  a_opt[,1]<-as.numeric(reg0$coefficients[2]>=0)
  regime_coef[1,]<-c(as.numeric(reg0$coefficients[2]>=0),rep(0,p))
  V_est[,1]<-reg0$coef[1]+reg0$coef[2]*a_opt[,1]
  err.rate[1]<-mean(a_opt[,1]!=g_true)
  S[1]<-sum(V_est[,1]-Y)
  
  #select the most important variable
  #S score
  S_temp<-rep(0,p);coef_all<-matrix(0,nrow=p,ncol=2)
  #critical_temp<-rep(0,p);direction_temp<-rep(0,p)
  V_temp<-matrix(0,nrow=n,ncol=p);err.rate_temp<-rep(0,p)
  for(j in 1:p){
    reg<-lm(Y~X[,j]*trt)
    coef<-reg$coef
    coef_all[j,]<-coef[3:4]
    diff<-coef[3]+coef[4]*X[,j]
    g_opt<-as.numeric(diff>=0)
    S_temp[j]<-sum(diff*(g_opt-a_opt[,1]))
    #critical_temp[j]<--coef[3]/coef[4]
    #direction_temp[j]<-as.numeric(coef[4]>=0)*2-1 #1 for >, -1 for <
    V_temp[,j]<-Y+diff*(g_opt-trt)
    err.rate_temp[j]<-mean(g_opt!=g_true)
  }
  #the most important variable
  imp<-which.max(S_temp);coef_temp<-coef_all[imp,]
  regime_coef[2,c(1,imp+1)]<-coef_temp
  a_opt[,2]<-as.numeric(cbind(1,X[,imp])%*%coef_temp>0)
  V_est[,2]<-V_temp[,imp];S[2]<-S_temp[imp];err.rate[2]<-err.rate_temp[imp]
  
  #time.start<-proc.time()[3]
  #select the rest important variables sequentially
  for (j in 2:step){
    S_temp<-rep(0,p) #in the jth search
    coef_all<-matrix(0,nrow=p,ncol=j+1)
    V_temp<-matrix(0,nrow=n,ncol=p)
    err.rate_temp<-rep(0,p)
    for (k in (1:p)[-imp]){
      #search the maximum, a_opt is the current trt regime
      X_temp<-X[,c(imp,k)]
      reg<-lm(Y~X[,c(imp,k)]*trt)
      coef<-reg$coef
      coef_all[k,]<-reg$coef[(j+2):(2*(j+1))]
      diff<-cbind(1,X[,c(imp,k)])%*%coef_all[k,]
      g_opt<-as.numeric(diff>=0)
      S_temp[k]<-sum(diff*(g_opt-a_opt[,j]))
      V_temp[,k]<-Y+diff*(g_opt-trt)
      err.rate_temp[k]<-mean(g_opt!=g_true)
    }
    #important variables
    j_opt<-which.max(S_temp)
    imp<-c(imp,j_opt);coef_temp<-coef_all[j_opt,]
    regime_coef[j+1,c(1,imp+1)]<-coef_temp
    a_opt[,j+1]<-as.numeric(cbind(1,X[,imp])%*%coef_temp>0)
    V_est[,j+1]<-V_temp[,j_opt];S[j+1]<-S_temp[j_opt]
    err.rate[j+1]<-err.rate_temp[j_opt]
    #time.stop<-proc.time()[3]
    #cat(j,"th execution time: ",time.stop-time.start,"s\n")
  }
  #time.stop<-proc.time()[3]
  #cat("execution time: ",time.stop-time.start,"s\n")
  
  out<-list(imp=imp,V_est=V_est,err.rate=err.rate,
            S=S,a_opt=a_opt,regime_coef=regime_coef)
  return(out)
}


#####################################################################









##########################################lose####################################
##########################################lose####################################
loss<-function(X,coef)
  
{
  
  
  g<-as.numeric(I(as.matrix(X)%*%coef>0))
  
  
  loss<-mean(wt*(b-g)^2)
  
  return(loss)
  
}

#######################################lose #######################################

split<-function(x)
{
  
  loss<-function(cutoff)
    
  {
    g<-as.numeric(I(x>cutoff))
    
    loss<-mean(wt*(b-g)^2)
    
    return(loss)
    
  }
  
  
  impurity<- min(do.call(rbind,lapply(x, loss)))
  idx<-which.min(lapply(x,loss))
  cutoff<--x[idx]
  
  return(data.frame(impurity,cutoff))
  
}
#######################################################################


rank.impurity<-function(X)
{
  


split<-function(x)
{


loss1<-function(cutoff)

{
g<-as.numeric(I(x>cutoff))

loss<-mean(wt*(b-g)^2)

return(loss)

}




loss0<-function(cutoff)

{
g<-as.numeric(I(x<cutoff))

loss<-mean(wt*(b-g)^2)

return(loss)

}



impurity1<- min(do.call(rbind,lapply(x, loss1)))
idx<-which.min(lapply(x,loss1))
cutoff1<-x[idx]


impurity0<- min(do.call(rbind,lapply(x, loss0)))
idx<-which.min(lapply(x,loss0))
cutoff0<-x[idx]


if(impurity1<impurity0)
{
impurity<-impurity1
cutoff<--cutoff1
coef<-1
}else

{
impurity<-impurity0
cutoff<-cutoff0
coef<--1
}




return(data.frame(impurity,cutoff,coef))



}

  
  cons<-data.frame(rep(1,n))
  names(cons)<-"cons"
  
  
  if(loss(cons,1)<loss(cons,-1))
  {impurity<-loss(cons,1)
  
  coef<-1
  } else
    
  {
    impurity<-loss(cons,-1)
    
    coef<--1
    
  }
  
  
  Xn<-dim(X)[2]
  
  Xname<-names(X)
  
  
  
  
  Xgini<-lapply(1:Xn,function(i) split(X[,i]))
  im.table<-cbind(Xname, do.call(rbind,Xgini))

  im.table$reduction<-impurity-im.table$imp
  im.table$percent<-im.table$reduc/impurity
  im.table<-im.table[order(im.table$imp),]
  
  
  return(im.table)
  
}

#######################################lose #######################################


##########################################add#######################################

add<-function(X0,x,coef)
{
  
  X<-cbind(X0,x)
  
  Xn<-dim(X)[2]
  
  
  
  Xs=X[,-1]
  
  glm.fit=glm(b~.,data=Xs,family="binomial",weight=wt)
  
  
  coef<-glm.fit$coef[-1]
  
  comx<-as.matrix(Xs)%*%as.numeric(coef)
  split.result<-split(comx)
  
  cons<-split.result$cutoff
  impurity<-split.result$impurity
  
  coef=c(cons,coef)
  
  
  
  result<-data.frame(matrix(c(impurity,coef),1,Xn+1))
  names(result)<-c("impurity",names(X0),"coef")
  return(result)
  
}
##########################################add######################################



#########################################forward###################################

forward<-function(X=X)
{
  
  cons<-data.frame(rep(1,n))
  names(cons)<-"cons"
  
  
  Xnames<-names(X)
  Xn<-dim(X)[2]
  
  X0.name<-Xnames[1]
  X0<-cbind(cons,X[X0.name])
  
  
  impurity<-rank[1,2]
  impurity0<-rank[1,2]+rank[1,5]
  
  coef<-rank[1,3:4]
  
  composite.var<-data.frame(impurity,coef)
  names(composite.var)<-c("impurity", "cons", names(X)[1])
  
  
  
  
  if(Xn>1)
  {
    
    X1.name<-Xnames[2:Xn]
    X1<-X[X1.name]
    
    
    
    Xn<-dim(X1)[2]
    
    finished<-0
    
    
    while(!finished)
    {
      
      im.ta<-lapply(1:Xn,function(i) add(X0=X0,x=X1[,i],coef=as.numeric(coef)))
      
      im.ta<-cbind(X1.name,do.call(rbind,im.ta))
      im.ta$reduction<-impurity-im.ta$imp
      im.ta$percent<-im.ta$reduc/impurity
      im.ta<-im.ta[order(im.ta$impur),]
      
      
      
      print(im.ta)
      
      
      
      if(im.ta$percent[1]<alpha | im.ta$imp[1] <0.01)
      {
        
        
        break
      }else
      {
        
        
        impurity<-im.ta$impurity[1]
        
        Xname<-as.vector(im.ta[,1])[1]
        xn<-dim(im.ta)[2]
        composite.var<-im.ta[1,(2:(xn-2))]
        names(composite.var)[length(composite.var)]<-c(Xname)
        
        
        coef<-as.numeric(composite.var[-1])
        X0<-cbind(X0,X1[Xname])
        
        
        if(dim(X0)[2]-1==dim(X)[2]) break
        
        if(impurity==0) break
        
        
        X1<-X1[-which(names(X1)%in%c(Xname))]
        X1.name<-names(X1)
        Xn<-dim(X1)[2]
        
      }
      
      
    }
  }
  
  
  return(composite.var[-1])
}




