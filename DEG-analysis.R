for (i in 1:30){
  print(i)
  datalist_new1[[i]]<-datalist_new1[[i]][,match(summaryrfs[[i]]$id,colnames(datalist_new1[[i]]))]}

pan_group<-list()
for (i in 1:30){
  print(i)
  pan_group[[cancername1[i]]]<-summaryrfs[[i]]
  pan_group[[cancername1[i]]]<-within(pan_group[[cancername1[i]]],{
    High<-NA
    Low<-NA
    High[MSI == 'MSI-H'] = '1'
    High[MSI == 'MSI-L'] = '0'
    Low[MSI == 'MSI-H'] = '0'
    Low[MSI == 'MSI-L'] = '1'
  })
  rownames(pan_group[[cancername1[i]]])<-pan_group[[cancername1[i]]]$id
  pan_group[[cancername1[i]]]$id<-NULL
  pan_group[[cancername1[i]]]<-pan_group[[cancername1[i]]][,-c(1:8)]
}

pan_group<-lapply(pan_group,function(df){
  df[] <- lapply(df, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else as.numeric(as.character(x))
  })
  return(df)
})
pan_group<-lapply(pan_group,as.matrix)


DEG<-list()
for (i in 27:30){
  print(i)
  fit <- lmFit(datalist_new1[[i]],pan_group[[i]])
  contrast.matrix <- makeContrasts(High - Low,levels=pan_group[[i]])
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  all_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 60000,sort.by = 'logFC')
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  DEG[[cancername1[i]]]<-all_diff
  DEG[[cancername1[i]]]$ID<-rownames(DEG[[cancername1[i]]])
}
cancernameDEG<-names(DEG)

DEGtest<-DEG
DEG<-DEGtest
for (i in 1:26){
  print(i)
  DEG[[i]]$cancer<-cancernameDEG[i]
  DEG[[i]]$gene<-rownames(DEG[[i]])
}

DEG<-lapply(DEG,function(df){
  df<-df%>%filter(df$adj.P.Val<0.05)
})

DEG<-lapply(DEG,function(df){
  df<-df%>%filter(abs(df$logFC)>0.5)
})

for (i in 1:26){
  print(i)
  DEG[[i]]<-within(DEG[[i]],{
    group<-NA
    group[logFC > 0 ] = 'up-regulated'
    group[logFC < 0 ] = '-down-regulated'
  })
}
DEG<-lapply(DEG,function(df){
  df$group<-as.factor(df$group)
  return(df)
})

DEGbar<-data.frame()
for (i in 1:26){
  print(i)
  DEGbar[2*i-1,1]<-cancernameDEG[i]
  DEGbar[2*i-1,2]<-summary(DEG[[i]]$group)[1]
  DEGbar[2*i-1,3]<-'High'
  DEGbar[2*i-1,4]<-summary(DEG[[i]]$group)[1]+summary(DEG[[i]]$group)[2]
  DEGbar[2*i,1]<-cancernameDEG[i]
  DEGbar[2*i,2]<-summary(DEG[[i]]$group)[2]
  DEGbar[2*i,3]<-'Low'
  DEGbar[2*i,4]<-summary(DEG[[i]]$group)[1]+summary(DEG[[i]]$group)[2]
}
colnames(DEGbar)<-c('cancer','value','Group','sum')

DEGbar%>%mutate(cancer = fct_reorder(cancer, desc(sum)))%>% #此处表示按反向排列！管道连接???
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ #表示手动调节颜色???
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))
write.csv(DEGbar,file = "/DEGbar.csv")
