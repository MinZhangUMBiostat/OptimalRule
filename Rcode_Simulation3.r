


rm(list=ls())
 
library(glmnet)


set.seed(12345)

source("othermethods_functions.r")
source("Rcode_CALBfunctions.r")

alpha=0.05
rho=0.2



################################
n=600
Monte<-500
s=500



p1<-50


beta = c(0.1,0,0,rep(0,6),1,1,rep(0,p1-10))
gamma = c(1,-0.8,1,0.9,0.8,1,0.9,0.8,rep(0,p1-8))


##############################




result1<-matrix(0,Monte,12)
result2<-matrix(0,Monte,12)
result3<-matrix(0,Monte,12)


generate<-function()
{

X = matrix(0,nrow=n,ncol=p1)                                           


X[,1] = rnorm(n)
  
for(j in 2:p1)
{
X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
}
  


p=0.5
A = rbinom(n,1,p)
epsilon = rnorm(n)
C_X=1*as.numeric(X[,7]^2-1*X[,8]^2+1*X[,9]+1*X[,10])
g_true = as.numeric(C_X>=0)

Y = 1+ as.numeric(X%*%gamma)+ A*C_X + epsilon
return(list(X=X,Y=Y,A=A,g_true=g_true))

}



data<-lapply(1:Monte,function(i) generate())




####################################################################################################################
mt=100000



testX = matrix(0,nrow=mt,ncol=p1)

testX[,1]=rnorm(mt)

for(j in 2:p1)
{
testX[,j] = rho*testX[,j-1]+sqrt(1-rho^2)*rnorm(mt)
}



C_test=1*as.numeric(1*testX[,7]^2-1*testX[,8]^2+1*testX[,9]+1*testX[,10])


g_test=as.numeric(I(C_test>0))



V_true=mean(1+ as.numeric(testX%*%gamma)+ g_test*C_test)




calV<-function(regime_coef)
{
X=testX
XX=X^2    
A_opt<-as.numeric(cbind(1,X,XX)%*%regime_coef>=0)
 
V_true=mean(1+ as.numeric(testX%*%gamma)+ A_opt*C_test)

err.rate= mean(A_opt!=g_test)
   
FP=sum(g_test[A_opt==1]==0)/sum(g_test==0)
FN=sum(g_test[A_opt==0]==1)/sum(g_test==1)
   
return(c(V_true,err.rate,FP,FN))
}








#########################################################
for(i in 1:s)
{


print("i=")
print(i)


a<-data[[i]]$A
y<-data[[i]]$Y

XX1<-data[[i]]$X
g_true<-data[[i]]$g_true




X<-data.frame(XX1)
m<-dim(X)[2]
names(X)[(1):(m)]<-paste("X",1:(m),sep='')




XX=data.frame(XX1^2)
m<-dim(XX)[2]
names(XX)[(1):(m)]<-paste("XX",1:(m),sep='')



#####################################################################
#                                SAS
#####################################################################

data.X=cbind(XX1,XX1^2)

step=30
result<-I_SS(X=data.X,Y=y,trt=a,g_true=g_true,step=30,rho=rho)




BIC<--log(cumsum(result$S))+(0:step)*log(n)/n
BIC_number1 <- which.min(BIC) - 1


Vresult<-calV(result$regime_coef[BIC_number1+1,])



size<-BIC_number1
TP<-sum(result$imp[1:BIC_number1]%in%c(9,10))

VR<-Vresult[1]/V_true
ER<-Vresult[2]


FP=Vresult[3]
FN=Vresult[4]


result1[i,1]<-size
result1[i,2]<-TP

result2[i,1]<-FP
result2[i,2]<-FN

result3[i,1]<-VR
result3[i,2]<-ER


#####################################################################
#                                Tian 
#####################################################################


t=2*a-1

W=X*(t/2)




data.X=cbind(X,XX)

cvfit = cv.glmnet(as.matrix(data.X), y, type.measure = "mse", nfolds=20)

coef1=coef(cvfit, s = "lambda.min")

Xn=dim(data.X)[2]

select=c(1:Xn)[coef1[-1]!=0]


coef=coef1[select+1]
Xnames=names(data.X)[select]
names(coef)=Xnames


mm=as.matrix(cbind(data.X)[,select])%*%coef+coef1[1]




data.X=cbind(1,X,XX)

data.X=data.X*(t/2)

cvfit = cv.glmnet(as.matrix(data.X), y-mm,intercept=FALSE, type.measure = "mse", nfolds=20)

coef1=coef(cvfit, s = "lambda.min")


Xn=dim(data.X)[2]-1

select=c(1:Xn)[coef1[-c(1,2)]!=0]



coef=coef1[select+2]
Xnames=names(data.X)[-1][select]
names(coef)=Xnames

beta0=coef1[2]

size<-length(Xnames)
TP<-sum(c(Xnames)%in%c("X9","X10","XX7","XX8"))

coef<-c(beta0,as.numeric(coef))
names(coef)=c("cons",Xnames)


Xn<-dim(data.X)[2]-1
count<-data.frame(t(rep(0,Xn+1)))
names(count)<-c("cons",names(data.X)[-1])
count[names(coef)]<-coef



Vresult<-calV(as.numeric(count))



VR<-Vresult[1]/V_true
ER<-Vresult[2]


FP=Vresult[3]
FN=Vresult[4]


result1[i,3]<-size
result1[i,4]<-TP

result2[i,3]<-FP
result2[i,4]<-FN

result3[i,3]<-VR
result3[i,4]<-ER

#####################################################################
#                                RegL
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

Xn=dim(AX)[2]
coef=coef1[select+1]
Xnames=names(data.trn)[select]
names(coef)=Xnames


m1=as.matrix(cbind(X,XX,1*X,XX,1)[,select])%*%coef+coef1[1]

m0=as.matrix(cbind(X,XX,0*X,0*XX,0)[,select])%*%coef+coef1[1]




AX.s=Xnames[Xnames%in%c(names(AX))]
AX.ss=substr(AX.s, 2, 6)


AX.s.coef=coef[AX.s]


AXX.s=Xnames[Xnames%in%c(names(AXX))]
AXX.ss=substr(AXX.s, 2, 6)



AXX.s.coef=coef[AXX.s]

Inter=coef1[1]
a.coef=0
if("a"%in%names(coef))
{a.coef=coef["a"]}



data.X=cbind(X,XX)



Xn<-dim(data.X)[2]
count<-data.frame(t(rep(0,Xn+1)))
names(count)<-c("cons",names(data.X))
count[AX.ss]<-AX.s.coef
count[AXX.ss]<-AXX.s.coef
count[1]=a.coef



Xnames=c(AX.ss,AXX.ss)


size<-length(AX.ss)+length(AXX.ss)
TP<-sum(c(Xnames)%in%c("X9","X10","XX7","XX8"))


Vresult<-calV(as.numeric(count))





VR<-Vresult[1]/V_true
ER<-Vresult[2]

FP=Vresult[3]
FN=Vresult[4]

result1[i,5]<-size
result1[i,6]<-TP

result2[i,5]<-FP
result2[i,6]<-FN


result3[i,5]<-VR
result3[i,6]<-ER


################################################################
#               fit the propensity score model
#################################################################


ph=mean(a)





#####################################################################
#                                ForMMER
#####################################################################

X.1<-X[a==1,]
X.0<-X[a==0,]

y.1<-y[a==1]
y.0<-y[a==0]



fit <- lm(y.0~ 1, data = X.0)


xnam <- names(X.0)
formula <- as.formula(paste(" ~ ", paste(xnam, collapse = "+")))
fit0=step(fit,scope=list(upper=formula,lower=~1),direction="forward",step=10,trace=0)




Xnames0<-rownames(summary(fit0)$coef)[-1]


m0_coef=summary(fit0)$coef[,1]
m0.lm<-as.matrix(cbind(1,X[Xnames0]))%*%summary(fit0)$coef[,1]



fit <- lm(y.1~ 1, data = X.1)


xnam <- names(X.1)
formula <- as.formula(paste(" ~ ", paste(xnam, collapse = "+")))
fit0=step(fit,scope=list(upper=formula,lower=~1),direction="forward",step=10,trace=0)



Xnames1<-rownames(summary(fit0)$coef)[-1]


m1_coef=summary(fit0)$coef[,1]
m1.lm<-as.matrix(cbind(1,X[Xnames1]))%*%summary(fit0)$coef[,1]




################################################################






ym<-a*y/ph-(a-ph)/ph*m1.lm-(1-a)*y/(1-ph)-(a-ph)/(1-ph)*m0.lm
b<-ifelse(ym>=0,1,0)
wt<-abs(ym)




n=length(b)

#################################################################################

data.X=cbind(X,XX)

rank<-rank.impurity(data.X)


Xnames<-as.vector((rank$Xname)[1:10])


XX2<-data.X[Xnames]




result<-forward(XX2)



Xnames<-names(result)[-1]



size<-length(Xnames)
TP<-sum(c(Xnames)%in%c("X9","X10","XX7","XX8"))




Xn<-dim(data.X)[2]
count<-data.frame(t(rep(0,Xn+1)))
names(count)<-c("cons",names(data.X))
count[names(result)]<-result


Vresult<-calV(as.numeric(count))




VR<-Vresult[1]/V_true
ER<-Vresult[2]

FP=Vresult[3]
FN=Vresult[4]

result1[i,7]<-size
result1[i,8]<-TP

result2[i,7]<-FP
result2[i,8]<-FN

result3[i,7]<-VR
result3[i,8]<-ER


 




#####################################################################
#                                CAL
#####################################################################




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



beta0=coef1[1]
size<-length(Xnames)
TP<-sum(c(Xnames)%in%c("X9","X10","XX7","XX8"))



coef<-c(beta0,as.numeric(coef))
names(coef)=c("cons",Xnames)

Xn<-dim(data.X)[2]
count<-data.frame(t(rep(0,Xn+1)))
names(count)<-c("cons",names(data.X))
count[names(coef)]<-coef




Vresult<-calV(as.numeric(count))




VR<-Vresult[1]/V_true
ER<-Vresult[2]

FP=Vresult[3]
FN=Vresult[4]

result1[i,9]<-size
result1[i,10]<-TP

result2[i,9]<-FP
result2[i,10]<-FN


result3[i,9]<-VR
result3[i,10]<-ER




#####################################################################
#                                CALB
#####################################################################

if(length(coef)>2)
{

Xs=data.X[Xnames]


b.w=b
wt.w=wt

Xs.w=Xs

X.w=X

coef.w=coef

n.w=n





alphas=c(0.05,0.1,0.15,0.2)

value_alphas=c(0,0,0,0)

 for (k in 1:length(alphas)){

   alpha1 = alphas[k];
   value=0;


samplen = sample(n.w);

folds=5



for (j in 1:folds){

print(j)

	tst_idx = samplen[seq(j,n,by=folds)]

	trn_idx = setdiff(1:n,tst_idx)


n=length(trn_idx)

wt=wt.w[trn_idx]
b=b.w[trn_idx]

Xs=Xs.w[trn_idx,]



result<-backward(Xs,coef.w[-1])



Xn<-length(result)
coef<-result
Xnames<-names(result)[2:(Xn)]





g<-as.numeric(I(as.matrix(cbind(1,Xs.w[Xnames]))%*%as.numeric(coef)>0))



value=value+sum(wt.w[tst_idx]*(b.w[tst_idx]-g[tst_idx])^2)

}

value_alphas[k]=value

}



alpha1=alphas[which.min(value_alphas)];


n=n.w

wt=wt.w
b=b.w

Xs=Xs.w
coef=coef.w


result<-backward(Xs,coef[-1])



Xnames<-names(result)[-1]

} else

{

result=coef
}



size<-length(Xnames)
TP<-sum(c(Xnames)%in%c("X9","X10","XX7","XX8"))



Xn<-dim(data.X)[2]
count<-data.frame(t(rep(0,Xn+1)))
names(count)<-c("cons",names(data.X))
count[names(result)]<-result


Vresult<-calV(as.numeric(count))





VR<-Vresult[1]/V_true
ER<-Vresult[2]

FP=Vresult[3]
FN=Vresult[4]

result1[i,11]<-size
result1[i,12]<-TP

result2[i,11]<-FP
result2[i,12]<-FN

result3[i,11]<-VR
result3[i,12]<-ER
###########################################################################


###########################################################################



############################################################################

}





mt3<-result3
mt2<-result2
mt1<-result1

i=1



for(i in c(5,7,1,3,9,11))
{

print(round(c(mean(mt1[,i]),mean(mt1[,i+1]),mean(mt2[,i])*100,mean(mt2[,i+1])*100,mean(mt3[,i+1])*100,sd(mt3[,i+1])*100,mean(mt3[,i])*100,sd(mt3[,i])*100    ),1))

}

