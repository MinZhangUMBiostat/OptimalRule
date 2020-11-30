


rm(list=ls()) 

source("Rcode_CALBfunctions.r")

library(glmnet)

set.seed(123456)







####################################################################
####################   Scenario I   ###############################
####################################################################
n=600


rho=0.2
alpha1=0.10


p1<-50



beta = c(0.1,0,0,rep(0,6),1,1,rep(0,p1-10))
gamma = c(1,-1,rep(0,p1-2))

##############################



X = matrix(0,nrow=n,ncol=p1)                                           


X[,1] = rnorm(n)
  
for(j in 2:p1)
{
X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
}
  

c=0.1+0.25*X[,1]+0.25*X[,5]

p=exp(c)/(1+exp(c))


a = rbinom(n,1,p)
epsilon = rnorm(n)
C_X=1*as.numeric(X[,7]^2+1.5*X[,8]^2+2*X[,9]+1.5*X[,10])
g_true = as.numeric(C_X>=0)


y =exp( 2+ as.numeric(X%*%gamma)-abs(1+1.5*X[,1]-2*X[,2])*(a-g_true)^2 ) + epsilon




#####################################################################

X<-data.frame(X)
m<-dim(X)[2]
names(X)[(1):(m)]<-paste("X",1:(m),sep='')

XX=data.frame(X^2)
m<-dim(XX)[2]
names(XX)[(1):(m)]<-paste("XX",1:(m),sep='')

##############################################################################


lasso_cv <- cv.glmnet(as.matrix(X), a, family = "binomial")
coef1=coef(lasso_cv, s = "lambda.min")


c=as.matrix(cbind(1,X))%*%coef1

ph=as.numeric(exp(c)/(1+exp(c)))




#####################################################################

m<-dim(X)[2]

AX=a*cbind(X)
names(AX)[(1):(m)]<-paste("AX",1:(m),sep='')


AXX=a*cbind(XX)
names(AXX)[(1):(m)]<-paste("AXX",1:(m),sep='')

data.trn=data.frame(X,XX,AX,AXX,a)


cvfit = cv.glmnet(as.matrix(data.trn), y, type.measure = "mse", nfolds=20)

coef1=coef(cvfit, s = "lambda.min")

Xn=dim(data.trn)[2]

select=c(1:Xn)[coef1[-1]!=0]


coef=coef1[select+1]




m1=as.matrix(cbind(X,XX,1*X,XX,1)[,select])%*%coef+coef1[1]

m0=as.matrix(cbind(X,XX,0*X,0*XX,0)[,select])%*%coef+coef1[1]




################################################################




ym<-a*y/ph-(a-ph)/ph*m1-(1-a)*y/(1-ph)-(a-ph)/(1-ph)*m0
b<-ifelse(ym>=0,1,0)
wt<-abs(ym)






#################################################################################

data.X=cbind(X,XX)


cvfit = cv.glmnet(as.matrix(data.X), b,weights=c(wt), family = "binomial", type.measure = "class", nfolds=20)



coef1=coef(cvfit, s = "lambda.min")

Xn=dim(data.X)[2]

select=c(1:Xn)[coef1[-c(1)]!=0]


coef=coef1[select+1]
Xnames=names(data.X)[select]
names(coef)=Xnames

#############################################################

Xs=data.X[Xnames]

result<-backward(Xs,coef)

coef<-result
Xnames<-names(result)[-1]


g<-as.numeric(I(as.matrix(cbind(1,data.X[Xnames]))%*%as.numeric(coef)>0))






####################################################################
####################   Scenario II   ###############################
####################################################################

n=600


rho=0.2
alpha1=0.10


p1<-1000



beta = c(0.1,0,0,rep(0,6),1,1,rep(0,p1-10))
gamma = c(1,-1,rep(0,p1-2))

##############################



X = matrix(0,nrow=n,ncol=p1)                                           


X[,1] = rnorm(n)
  
for(j in 2:p1)
{
X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
}
  

a = rbinom(n,1,0.5)
epsilon = rnorm(n)
C_X=as.numeric(cbind(1,X)%*%beta)
g_true = as.numeric(C_X>=0)


y =exp( 2+ as.numeric(X%*%gamma)-abs(1+1.5*X[,1]-2*X[,2])*(a-g_true)^2 ) + epsilon




#####################################################################

X<-data.frame(X)
m<-dim(X)[2]
names(X)[(1):(m)]<-paste("X",1:(m),sep='')



##############################################################################


ph=mean(a)




#####################################################################
m<-dim(X)[2]

AX=a*cbind(X)
names(AX)[(1):(m)]<-paste("AX",1:(m),sep='')


data.trn=data.frame(X,AX,a)


cvfit = cv.glmnet(as.matrix(data.trn), y, type.measure = "mse", nfolds=20)

coef1=coef(cvfit, s = "lambda.min")

Xn=dim(data.trn)[2]

select=c(1:Xn)[coef1[-1]!=0]


coef=coef1[select+1]



m1=as.matrix(cbind(X,1*X,1)[,select])%*%coef+coef1[1]

m0=as.matrix(cbind(X,0*X,0)[,select])%*%coef+coef1[1]



################################################################




ym<-a*y/ph-(a-ph)/ph*m1-(1-a)*y/(1-ph)-(a-ph)/(1-ph)*m0
b<-ifelse(ym>=0,1,0)
wt<-abs(ym)






#################################################################################

data.X=cbind(X)



cvfit = cv.glmnet(as.matrix(data.X), b,weights=c(wt), family = "binomial", type.measure = "class", nfolds=20)



coef1=coef(cvfit, s = "lambda.min")

Xn=dim(data.X)[2]

select=c(1:Xn)[coef1[-c(1)]!=0]


coef=coef1[select+1]
Xnames=names(data.X)[select]
names(coef)=Xnames

#############################################################

Xs=data.X[Xnames]

result<-backward(Xs,coef)

coef<-result
Xnames<-names(result)[-1]


g<-as.numeric(I(as.matrix(cbind(1,X[Xnames]))%*%as.numeric(coef)>0))

