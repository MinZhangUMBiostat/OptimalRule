#==============================================#
#         README.txt                           #
#==============================================#

Code for the proposed methods in `` Subgroup Identification and Variable Selection  for Treatment Decision Making'' by Baqun Zhang and Min Zhang

1. General data setting
   This code is for the data described as in Section of simulation study in the paper. Briefly, the followings are assumed:

   a. Randomized clinical trial data  or observational data
   b. One decision point
   c. two treatments A=1 or 0 (e.g. treatment v.s. control) are available. 
   d. A set of covariate X is given ; 
   e. Y: response or outcome; Large values are more preferable.


2. Simulated data setting

   a. Treatment a generated from Bernoulli
distributions with probability model as specified below or with probability 0.5.
   b. Markers X were generated from multivarate normal distribution .
   c. Scenario I:Response was generated from normal distribution with E(Y|Treatment, X)=exp(2+X1-X2-|1+1.5*X1-2*X2|*(a-g)^2) where g=I(X7^2+ 1.5X8^2+ 2X9+ 1.5X10>0),  logit{¦Ð(X)} =0.1 + 0.25X1+ 0.25X5
     Scenario II:Response was generated from normal distribution with E(Y|Treatment, X)=exp(2+X1-X2-|1+1.5*X1-2*X2|*(a-g)^2) where g=I(0.1+X9+X10>0), ¦Ð(X)=0.5

 


3. R library to be installed: 
   To run this code, you must install following R library:

   a.library(glmnet) 


4. Run the example code "Rcode_CALB.r".

   a.Include r file "Rcode_CALBfunctions.r" by 
source("your_folder:/Rcode_CALBfunctions.r")



