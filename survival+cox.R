for (i in 1:31){
  print(i)
  res.cut<-surv_cutpoint(summarylist_rfs[[i]], time = "rfs.time", event = "rfs",
                         variables = "MSIscore")
  res_cat<-surv_categorize(res.cut)
  summarylist_rfs[[i]]<-cbind(summarylist_rfs[[i]],res_cat$score)
  names(summarylist_rfs[[i]])[8]<-"MSI_group"
}

rfsname<-names(summarylist_rfs)

library(regplot)
library(survival)
library(survminer)
data1<-read.csv()
cox.1 <- coxph(Surv(time,status)~gender+age+stage+score,data= LUADdata)
obs<-data1[5,]
nomo<-regplot(cox.1,observation=obs, plots = c("density","boxes"),failtime =c(12,36,60), prfail = TRUE,
              boxcol="#ADD8E6",cexvars=1.2,cexscales=1.2,cexcats=1.0,droplines=TRUE,
              title = "Nomogram to predict CSS in LUAD patients")