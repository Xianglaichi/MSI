miRids<-read.csv("/miRids.csv")
miRids<-miRids[!apply(miRids == "", 1, all),]

miRdata<-list()
for(i in 1:length(cancername1)){
  print(i)
  miRdata[[cancername1[i]]] = read.table(paste("/",cancername1[i],'_miR.tsv',sep = ''),
                                         sep = '\t',
                                         header = T,
                                         quote = '',
                                         fill = T, 
                                         comment.char = "!",
                                         stringsAsFactors = FALSE)
}
BIrownames<-function(df){
  rownames(df)<-df[,1]
  df[,1]<-NULL
  return(df)
}
miRdata<-lapply(miRdata,BIrownames)
miR_group<-lapply(pan_group,as.data.frame)

miRdata<-lapply(miRdata,function(df){
  colnames(df)<-paste0(colnames(df),"A")
  return(df)
})

for (i in 1:30){
  print(i)
  miRdata[[i]]<-miRdata[[i]][,colnames(miRdata[[i]]) %in% rownames(miR_group[[i]])]
}
for (i in 1:30){
  print(i)
  miR_group[[i]]<-miR_group[[i]][rownames(miR_group[[i]]) %in% colnames(miRdata[[i]]),]
}
for (i in 1:30){
  print(i)
  miRdata[[i]]<-miRdata[[i]][,match(rownames(miR_group[[i]]),colnames(miRdata[[i]]))]
}
miRdata<-lapply(miRdata,as.matrix)
miR_group<-lapply(miR_group,as.matrix)

miRDEG<-list()
miR_diff<-data.frame()
for (i in 27:30){
  print(i)
  fit <- lmFit(miRdata[[i]],miR_group[[i]])
  contrast.matrix <- makeContrasts(High - Low,levels=miR_group[[i]])
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  miR_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 100000,sort.by = 'logFC')
  miR_diff$logP<- -log10(miR_diff$adj.P.Val)
  miRDEG[[cancername1[i]]]<-miR_diff
  miRDEG[[cancername1[i]]]$ID<-rownames(miRDEG[[cancername1[i]]])
}
cancernamemir<-names(miRDEG)

DEGtest<-DEG
for (i in 1:25){
  print(i)
  miRDEG[[i]]$cancer<-cancernamemir[i]
}

miRDEG<-lapply(miRDEG,function(df){
  df<-df%>%filter(df$adj.P.Val<0.05)
})

miRDEG<-lapply(miRDEG,function(df){
  df<-df%>%filter(abs(df$logFC)>0.5)
})

for (i in 1:25){
  print(i)
  miRDEG[[i]]<-within(miRDEG[[i]],{
    group<-NA
    group[logFC > 0 ] = 'high'
    group[logFC < 0 ] = 'low'
  })
}
miRDEG<-lapply(miRDEG,function(df){
  df$group<-as.factor(df$group)
  return(df)
})

miRDEGbar<-data.frame()
for (i in 1:25){
  print(i)
  miRDEGbar[2*i-1,1]<-cancernamemir[i]
  miRDEGbar[2*i-1,2]<-summary(miRDEG[[i]]$group)[1]
  miRDEGbar[2*i-1,3]<-'High'
  miRDEGbar[,2][is.na(miRDEGbar[,2])] <- 0
  miRDEGbar[2*i,1]<-cancernamemir[i]
  miRDEGbar[2*i,2]<-summary(miRDEG[[i]]$group)[2]
  miRDEGbar[2*i,3]<-'Low'
  miRDEGbar[,2][is.na(miRDEGbar[,2])] <- 0
  miRDEGbar[2*i-1,4]<-miRDEGbar[2*i-1,2]+miRDEGbar[2*i,2]
  miRDEGbar[2*i,4]<-miRDEGbar[2*i-1,2]+miRDEGbar[2*i,2]
}
colnames(miRDEGbar)<-c('cancer','value','Group','sum')

DEGbar1<-DEGbar
DEGbar1<-as.data.frame(t(DEGbar1))
DEGbar1$`21`<-NULL
DEGbar1$`22`<-NULL
miRDEGbar$rank<-DEGbar1$sum

miRDEGbar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% 
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))
