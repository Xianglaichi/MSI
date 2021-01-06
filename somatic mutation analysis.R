DEGbar<-read.csv("/DEGbar.csv")
DEGbar1<-DEGbar
library(maftools)
MAFdata<-list()
MAFdatahigh<-list()
MAFdatalow<-list()
for(i in 29:length(cancername1)){
  print(i)
  MAFdata[[cancername1[i]]] = read.maf(paste("/",cancername1[i],'.maf',sep = ''))
  datas<-MAFdata[[cancername1[i]]]@data
  datas<-as.data.frame(datas)
  datas$Tumor_Sample_Barcode<-gsub('-','.',datas$Tumor_Sample_Barcode)
  datas$Tumor_Sample_Barcode<-gsub('............$','',datas$Tumor_Sample_Barcode)
  MAFdatahigh[[cancername1[i]]] = datas[datas$Tumor_Sample_Barcode %in% rownames(pan_group[[cancername1[i]]])[which(pan_group[[cancername1[i]]]$High == 1)],]
  MAFdatahigh[[cancername1[i]]] = read.maf(MAFdatahigh[[cancername1[i]]])
  MAFdatalow[[cancername1[i]]] = datas[datas$Tumor_Sample_Barcode %in% rownames(pan_group[[cancername1[i]]])[which(pan_group[[cancername1[i]]]$High == 0)],]
  MAFdatalow[[cancername1[i]]] = read.maf(MAFdatalow[[cancername1[i]]])
}

MAFdata_compare<-list()
for (i in 17:length(cancername)){
  print(i)
  MAFdata_compare[[cancername[i]]] <- mafCompare(m1 = MAFdatahigh[[cancername[i]]], m2 = MAFdatalow[[cancername[i]]], 
                                                 m1Name = 'High', m2Name = 'Low', minMut = 5)
}
MAFdata_compare_test<-MAFdata_compare

MAFdata_compare_summary<-lapply(MAFdata_compare,function(df){
  df<-as.data.frame(df$results)
  return(df)
})

MAFdata_compare_summary<-lapply(MAFdata_compare_summary,function(df){
  df<-df%>%filter(df$adjPval<0.05)
})

for (i in 1:30){
  print(i)
  MAFdata_compare_summary[[i]]<-within(MAFdata_compare_summary[[i]],{
    group<-NA
    group[or > 1] = 'high'
    group[or < 1] = 'low'
  })
}
MAFdata_compare_summary<-lapply(MAFdata_compare_summary,function(df){
  df$group<-as.factor(df$group)
  return(df)
})

MAFdata_compare_summarybar<-data.frame()
for (i in 1:30){
  print(i)
  MAFdata_compare_summarybar[2*i-1,1]<-cancername[i]
  MAFdata_compare_summarybar[2*i-1,2]<-summary(MAFdata_compare_summary[[i]]$group)[1]
  MAFdata_compare_summarybar[2*i-1,3]<-'High'
  MAFdata_compare_summarybar[,2][is.na(MAFdata_compare_summarybar[,2])] <- 0
  MAFdata_compare_summarybar[2*i,1]<-cancername[i]
  MAFdata_compare_summarybar[2*i,2]<-summary(MAFdata_compare_summary[[i]]$group)[2]
  MAFdata_compare_summarybar[2*i,3]<-'Low'
  MAFdata_compare_summarybar[,2][is.na(MAFdata_compare_summarybar[,2])] <- 0
  MAFdata_compare_summarybar[2*i-1,4]<-MAFdata_compare_summarybar[2*i-1,2]+MAFdata_compare_summarybar[2*i,2]
  MAFdata_compare_summarybar[2*i,4]<-MAFdata_compare_summarybar[2*i-1,2]+MAFdata_compare_summarybar[2*i,2]
}
colnames(MAFdata_compare_summarybar)<-c('cancer','value','Group','sum')
MAFdata_compare_summarybar$rank<-DEGbar$sum

MAFdata_compare_summarybar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% 
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ 
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))