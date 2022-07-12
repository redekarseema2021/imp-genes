install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)

setwd("F:\\PHD (1)\\Cancer")
#cdf=fread('Significant genes for network construction1.csv')
key_gene<-read.csv("F:\\PHD (1)\\Cancer\\final_key.csv", sep=",")
View(key_gene)
summary(key_gene)
attach(key_gene)
library(survival)
install.packages("survminer")
library(survminer)
list_new<-colnames(key_gene)[4:13]
list_new
x1<-NULL
x2<-NULL
for (i in 1:length(list_new)){
  print(list_new[i])
  res.cox <- coxph(Surv(OS.time, OS) ~ get(list_new[i]), data=key_gene)
  print(summary(res.cox))
  #res.cox.sum <- data.frame(summary(res.cox)$coefficients)
 # print(res.cox.sum[,5])
  #if (res.cox.sum[, 5] <= 0.01) {
  # print(list_new[i])
  #print(res.cox.sum[,5])
  # x1 <- c(x1, list_new[i])
  #x2 <- c(x2, res.cox.sum[, 5])
  #}
}
x1_name <- "Gene"
x2_name <- "P-value"   
df <- data.frame(x1,x2)
names(df) <- c(x1_name,x2_name)
print(df)
write.csv(df,"F:\\PHD (1)\\dataset\\TCGA_myth\\TCGA-GBM.methylation27_final1_new_30_pvalu0.05_subset_final_april.csv", row.names = FALSE)



key_gene<-read.csv("F:\\PHD (1)\\Cancer\\final_key.csv", sep=",")
View(key_gene)
res.cox1 <- coxph(Surv(OS.time, OS) ~ SLC44A2.+KRTAP4.2, data=key_gene)
print(summary(res.cox1))


res.cox2 <- coxph(Surv(OS.time, OS) ~ SEMA3B+APS+MARK2+PITPNM2.+SFRP1+PRLH+DIP2C+CTSZ, data=key_gene)
print(summary(res.cox2))
