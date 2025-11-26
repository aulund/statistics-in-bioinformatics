rm(list=ls())
PSA=c(3.8, 3.4, 2.9, 2.8, 2.7, 2.1, 1.6,2.5, 2.0, 1.7, 1.4, 1.2, 0.9, 0.8) 
Group=factor(rep(c("Cancer","Healthy"),c(7,7)),levels=c("Healthy","Cancer"))
df=data.frame(Group,PSA)
#install.packages("neuralnet")
library(neuralnet)
nn=neuralnet(Group~PSA,data=df, hidden=0, act.fct = "logistic",
             linear.output = FALSE, err.fct="ce", threshold = 1e-6)
print(nn)
plot(nn)
df2=data.frame(PSA=2)
predict(nn,df2)

# Logistic regression
glm(Group~PSA, data=df, family=binomial)

# Plot
state=ifelse(Group=="Cancer",1,0) # Recode to Cancer -> 1, Healthy -> 0
plot(PSA,state,col=Group,xlab="PSA µg/L",ylab="Probability",xlim=c(0,5))
x=seq(0,5,0.01)
pred=predict(nn,data.frame(PSA=x))[,1]
lines(x,pred,col="red")
