###########################################################################

loss<-function(X,coef)
  
{
  
  
  g<-as.numeric(I(as.matrix(X)%*%coef>0))
  
  
  loss<-mean(wt*(b-g)^2)
  
  return(loss)
  
}


########################################################################

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


#########################################################################

backward<-function(X,coef)
{


Xnames<-names(X)


finished<-0
Xn<-dim(X)[2]
if (Xn<2) finished<-1


while(!finished)
{

Xn<-dim(X)[2]


comx<-as.matrix(X)%*%as.numeric(coef)
split.result<-split(comx)

cons<-split.result$cutoff
impurity<-split.result$impurity

result<-c(impurity,cons,coef)
names(result)<-c("err","cons",Xnames)


composite.var<-result
print(composite.var)



coef<-backward.delete(X=X,composite.var=composite.var)
Xnames<-names(coef)
X<-X[Xnames]
Xn.new<-length(Xnames)
if(Xn.new==Xn | Xn.new==1)

{ finished=1}


}







Xn<-dim(X)[2]

if(Xn == 1){

Xs=X

glm.fit=glm(b~.,data=Xs,family="binomial",weight=wt)

coef<-glm.fit$coef[-1]

comx<-as.matrix(X)%*%as.numeric(coef)
split.result<-split(comx)

cons<-split.result$cutoff
impurity1<-split.result$impurity

result<-c(impurity1,cons,as.numeric(coef))
names(result)<-c("err","cons",Xnames)

composite.var<-result
print(composite.var)


}


if(Xn == 0){


ones<-data.frame(rep(1,n))
names(ones)<-"cons"


if(loss(ones,1)<loss(ones,-1))
{
impurity0<-loss(ones,1)
cons=1
} else{

impurity0<-loss(ones,-1)
cons=-1
}



result<-c(impurity0,cons)
names(result)<-c("err","cons")
composite.var<-result
print(composite.var)


}



return(composite.var[-1])



}


############################################################################


del.impurity<-function(i,X,coef,cons)
{
  
  
  coef<-as.matrix(coef)
  
  del.imp<-loss(cbind(1,X[,-i]),c(cons,coef[-i]))
  return(data.frame(del.imp))
}




########################################################

backward.delete<-function(X,composite.var)
{
  
  
  Xnames<-names(composite.var)[-c(1,2)]
  coef<-composite.var[-c(1,2)]
  impurity<-as.numeric(composite.var[1])
  cons<-as.numeric(composite.var[2])
  
  Xn<-length(Xnames)
  del.imp.table<-lapply(1:Xn,function(i) del.impurity(i,X=X,coef=coef,cons=cons))
  
  im.table<-cbind(Xnames,do.call(rbind,del.imp.table))
  im.table$deterioration<-im.table[,2]-as.numeric(impurity)
  im.table<-im.table[order(im.table$deter,decreasing = TRUE),]
  im.table$percent<-im.table$deterioration/im.table$deterioration[1]
  im.table$status<-I(im.table$percent>=alpha1)
  
  
  
  print(im.table)
  Xnames<-as.vector(im.table$Xnames[im.table$percent>=alpha1])
  coef<-coef[Xnames]
  
  
  return(coef)
  
  
  
}



