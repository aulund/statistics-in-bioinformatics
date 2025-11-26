rm(list=ls())
Protein<-c(3.8, 3.4, 2.9, 2.8, 2.7,    2.1, 2.0, 1.7, 1.6, 1.4,    0.9, 0.8, 0.2,0.2,0.1) # PSA
Group<-factor(rep(c("Cancer","Healthy","Cancer"),c(5,5,5)),levels=c("Healthy","Cancer"))
df<-data.frame(Group,Protein)
state=ifelse(Group=="Cancer",1,0) # Recode to Cancer -> 1, Healthy -> 0


set.seed(500)
library(neuralnet)
# fit neural network
nn=neuralnet(Group~Protein,data=df, hidden=c(2,2),act.fct = "logistic",
             linear.output = FALSE,err.fct="ce",threshold = 1e-8, rep=1)
plot(nn)


plot(Protein,state,col=Group,xlab="",ylab="Probability",xlim=c(0,5))
legend("bottomright",levels(Group),col=1:2,pch=1)
x=seq(0,5,0.01)
pred=predict(nn,data.frame(Protein=x))
lines(x,pred[,1],col="red")

