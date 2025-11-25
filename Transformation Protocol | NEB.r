install.packages(c("reshape","pROC","MASS", "pvclust","tibble","ggplot2","dplyr","tidyr"))
load(c("reshape","pROC","MASS", "pvclust","tibble","ggplot2","dplyr","tidyr"))


df <- read.table(file.choose(),header=TRUE,sep=";",row.names=1)

mode(df)
print(df)

df["STPG1", 8]


library(tibble)
tib = as_tibble(df)
tib

tdf = t(df)

tdf=data.frame(tdf)
tdf

Group = rep(c("CVID", "HD"),c(9,9))
Group

Group=factor(Group)

Group=="CVID"


tdf$STPG1[Group=="CVID"]


library("dplyr")
select(tdf, STPG1) %>% filter(Group=="CVID")

mean(tdf$STPG1)
median(tdf$STPG1)
sd(tdf$STPG1)
IQR(tdf$STPG1)


median(tdf$FUCA2)
#174.5

mean(tdf$STPG1[Group=="CVID"])

tdf %>% filter(Group=="CVID") %>% pull(STPG1) %>% mean()


IQR(tdf$ANKIB1[Group=="HD"])
#387


apply(tdf,2,median)
apply(tdf,1,median)

tdf %>% filter(Group == "CVID") %>% summarise(across(where(is.numeric), median))


mean_CVID=apply(tdf[Group=="CVID",],2,mean)
mean_HD=apply(tdf[Group=="HD",],2,mean)
FC=mean_CVID/mean_HD
FC

log2FC=log2(FC)

sumdf=data.frame(mean_CVID, mean_HD, log2FC)
round(sumdf,2)

sort(log2FC, decreasing = TRUE)
# upregulated = HS3ST1         CFH  
#downregulated = (CFTR) SEMA3F DPM1

boxplot(tdf, las=3, ylab="Counts")

library(tidyr)
library(ggplot2)

tdf_long <- tdf%>%
  pivot_longer(cols =where(is.numeric),
               names_to = "variable",
               values_to="value")
ggplot(tdf_long, aes(x=variable, y=value)) +
  geom_boxplot()+
  labs(y = "Counts", x = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        
        

boxplot(df, las=3, ylab="Counts")

boxplot(tdf$LAS1L~Group, las=2,ylab="Read count", names=c("CVID", "HD"), main="LAS1L")

tdf_long <- tdf %>%
  mutate(Group = Group) %>% pivot_longer(
    cols = where(is.numeric),
    names_to = "variable",
    values_to = "value"
    # add group vector temporarily
  )
ggplot(tdf_long, aes(x = variable, y = value, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(y = "Counts", x = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


stripchart(tdf$LAS1L~Group,vertical =T, method="jitter", pch=1, col=c(1:2), ylab="Read count", main="LAS1L")

ggplot(tdf_long, aes(x = variable, y = value, color = Group)) +
  geom_jitter(size=2,width=0.1,height=0,alpha = 0.5) +
  labs(y = "Counts", x = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

par(mfrow=c(2,1))
  hist(tdf$LAS1L[Group=="CVID"],main="LAS1L CVID", xlab=c("Read Count"))
  hist(tdf$LAS1L[Group=="HD"],main="LAS1L HD",xlab=c("Read count"))     

  
t.test(tdf$TSPAN6~Group, var.equal = FALSE)

t.test(tdf$TSPAN6~Group, var.equal = FALSE)$p.value

pvals=apply(tdf,2,function(x){t.test(x~Group, var.equal = FALSE)$p.value})     

sumdf$pvals = pvals

sumdf_order_p=sumdf[order(sumdf$pvals),]
round(sumdf_order_p,5)

rownames(sumdf[sumdf$pvals <= 0.05, ])
#[1] "DPM1"    "SCYL3"   "FGR"     "FUCA2"   "NFYA"   
#[6] "NIPAL3"  "LAS1L"   "ENPP4"   "SEMA3F"  "ANKIB1" 
#[11] "CYP51A1" "LAP3"    "HS3ST1" 

sumdf$pvals_BF=pvals*25
sumdf_order_p=sumdf[order(sumdf$pvals),]
round(sumdf_order_p,5)

rownames(sumdf_order_p[sumdf_order_p$pvals_BF <= 0.05,])

> rownames(sumdf_order_p[sumdf_order_p$pvals_BF < 0.05,])
#[1] "ANKIB1"  "NFYA"    "HS3ST1"  "NIPAL3"  "FGR"    
#[6] "LAS1L"   "DPM1"    "CYP51A1" "LAP3"   

sumdf$pvals_BH=p.adjust(pvals,method="BH")
sumdf_order_p=sumdf[order(sumdf$pvals),]
round(sumdf_order_p, 5)

rownames(sumdf_order_p[sumdf_order_p$pvals_BH <= 0.05,])
#[1] "ANKIB1"  "NFYA"    "HS3ST1"  "NIPAL3"  "FGR"    
#[6] "LAS1L"   "DPM1"    "CYP51A1" "LAP3"    "ENPP4"  
#[11] "FUCA2"   "SEMA3F"  "SCYL3"  

plot(tdf$NFYA, tdf$ANKIB1, xlab="counts of NFYA", ylab="counts of ANKIB1", col=Group)
legend("topleft", legend=levels(Group), pch=1, col=1:2)

cor.test(tdf$LAS1L,tdf$DPM1, method = "pearson")

cor(tdf, method = "spearman")
