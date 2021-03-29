#==============================================#
#         README.txt                           #
#==============================================#

####################################################### enclosed in this folder ###########
   Rcode_Simulation1, Rcode_Simulation2, Rcode_Simulation3: code used for simulation studies (including cross-validation for choosing alpha)
   Rcode_CALB: illustrating code for the proposed method CALB (for a fixed alpha). Please refer to simulation code for choosing alpha using cross-validation
   Rcode_CALBfunctions: Defined functions used for implementing the proposed method
   othermethods_functions: Defined functions used for other methods.
###########################################################################################






This README file will show you the R script Rcode_CALB.r  and Rcode_CALBfunctions.r details.
The goal of this R script is to guide readers to apply the method proposed in the paper, 
"SUBGROUP IDENTIFICATION AND VARIABLE SELECTION FOR TREATMENT DECISION MAKING", written by Baqun zhang, and Min zhang, submitted to Annals of Applied Statistics in 2020. 
Notice that this R script  is generated 
for our  simulated data as described in the paper and below. Therefore, one should edit properly some part of the scripts  according to user's data setting.  
The authors do our best to avoid potential errors and to help readers to use the provided R scripts, but please note that the authors can not take any responsibilities for any problems resulted from the R scripts.


1. General data setting
   This code is for the data described as in Section of simulation study in the paper. Briefly, the followings are assumed:

   a. Randomized clinical trial data  or observational data
   b. One decision points
   c. two treatments A=1 or 0 (e.g. treatment v.s. no treatment) are available. 
   d. A set of markers X (continuous) is given ; 
   e. A outcome Y is given as response variable; Large value of outcome is more preferable.


2. Simulated data setting

   a. Treatment a generated from Bernoulli
distributions with probability pi(X) as specified below or with probability 0.5.
   b. Markers X were generated from multivarate normal distribution .
   c. Scenario I:Response was generated from normal distribution with E(Y|Treatment, X)=exp(2+X1-X2-|1+1.5*X1-2*X2|*(a-g)^2) where g=I(X7^2+ 1.5X8^2+ 2X9+ 1.5X10>0),  logit{¦Ð(X)} =0.1 + 0.25X1+ 0.25X5
      Scenario II:Response was generated from normal distribution with E(Y|Treatment, X)=exp(2+X1-X2-|1+1.5*X1-2*X2|*(a-g)^2) where g=I(0.1+X9+X10>0), pi(X)=0.5
      Scenario III:Response was generated from normal distribution with E(Y|Treatment, X)=1+X1-0.8X2+X3+0.9X4+0.8X5+X6+0.9X7+0.8X8+a(X7^2-X8^2+ X9+ X10), pi(X)=0.5



3. R library to be installed: 
   To run this code, you must install following R library:

   a.library(glmnet) 


4. Run the example code "Rcode_CALB.r".

   a.Include r file "Rcode_CALBfunctions.r" by 
source("your_folder:/Rcode_CALBfunctions.r")






