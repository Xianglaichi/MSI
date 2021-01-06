DEGsummary<-data.frame(DEG$ACC)
for (i in 1:29){
  print(i)
  DEGsummary<-rbind(DEGsummary,DEG[[i+1]])
}
DEGsummary$ID_group<-paste(DEGsummary$ID,DEGsummary$group,sep = "_")
class(DEGsummary$ID_group)
DEGsummary$ID_group<-as.factor(DEGsummary$ID_group)
DEGcommon<-as.data.frame(summary(DEGsummary$ID_group,maxsum = 100000))
DEGcommon$ID_group<-rownames(DEGcommon)
DEGcommon<-DEGcommon%>%filter(DEGcommon$`summary(DEGsummary$ID_group, maxsum = 1e+05)`>=10)
DEGcommon<-read.csv("/DEGcommon.csv")
commongene<-data.frame(DEGcommon$ID_group)
commongene_final<-data.frame(gsub("_.*$","",commongene$DEGcommon.ID_group))
names(commongene_final)<-"ID"

GDSCexp<-as.data.frame(GDSCexp)
GDSCexp_filter<-GDSCexp[rownames(GDSCexp) %in% commongene_final$ID,]
GDSCexp_filter<-as.data.frame(t(GDSCexp_filter))

GDSC_gene_drug_cor<-as.data.frame(matrix(1:114395,ncol=167,nrow=685))
colnames(GDSC_gene_drug_cor)<-colnames(GDSCsummary[,-c(168:169)])
rownames(GDSC_gene_drug_cor)<-colnames(GDSCexp_filter)

GDSC_gene_drug_cor[!is.na(GDSC_gene_drug_cor)]<-NA 
for (i in 1:167){
  print(i)
  for (j in 1:685){
    print(j)
    GDSC_gene_drug_cor[j,i] = cor(GDSCsummary[i],GDSCexp_filter[j],use = "complete.obs")
  }
}
test<-as.data.frame(apply(abs(GDSC_gene_drug_cor),2,max))
test<-test[order(test$`apply(abs(GDSC_gene_drug_cor), 2, max)`),]
test<-as.data.frame(test)


GDSC_gene_drug_cor_test<-GDSC_gene_drug_cor
GDSC_gene_drug_cor<-GDSC_gene_drug_cor_test
GDSC_gene_drug_cor$gene<-NULL
GDSC_gene_drug_cor<-apply(GDSC_gene_drug_cor,2,function(x){
  x<-ifelse(abs(x)<= 0.3,NA,x) 
})
GDSC_gene_drug_cor<-as.data.frame(t(GDSC_gene_drug_cor))

hx<-apply(GDSC_gene_drug_cor,2,function(x){
  if(summary(is.na(x))[[2]]>=3) print("x") else x<-NA
}) 
hx<-as.data.frame(hx)
hx$ID<-rownames(hx)
hx=hx%>%filter(is.na(hx)==FALSE)
GDSC_gene_drug_cor$drug<-rownames(GDSC_gene_drug_cor)