TCGAprotein<-read.table("/TCGAprotein.txt",
                        sep = '\t',
                        header = T,
                        quote = '',
                        fill = T, 
                        comment.char = "!",
                        stringsAsFactors = FALSE)
rownames(TCGAprotein)<-TCGAprotein$SampleID
TCGAprotein[,1]<-NULL
colnames(TCGAprotein)<-paste(colnames(TCGAprotein),"A",sep = '')

pan_group<-lapply(pan_group,as.data.frame)

Proteindata<-list()
for (i in 1:30){
  print(i)
  Proteindata[[cancername1[i]]]<-TCGAprotein[colnames(TCGAprotein) %in% rownames(pan_group[[cancername1[i]]])]
}

Protein_group<-lapply(pan_group,as.data.frame)

for (i in 1:30){
  print(i)
  Proteindata[[i]]<-Proteindata[[i]][,colnames(Proteindata[[i]]) %in% rownames(Protein_group[[i]])]
}
for (i in 1:30){
  print(i)
  Protein_group[[i]]<-Protein_group[[i]][rownames(Protein_group[[i]]) %in% colnames(Proteindata[[i]]),]
}
for (i in 1:30){
  print(i)
  Proteindata[[i]]<-Proteindata[[i]][,match(rownames(Protein_group[[i]]),colnames(Proteindata[[i]]))]
}
Proteindata<-lapply(Proteindata,as.matrix)
Protein_group<-lapply(Protein_group,as.matrix)


ProteinDEG<-list()
Protein_diff<-data.frame()
for (i in 1:24){
  print(i)
  fit <- lmFit(Proteindata[[cancernamepr[i]]],Protein_group[[cancernamepr[i]]])
  contrast.matrix <- makeContrasts(High - Low,levels=Protein_group[[cancernamepr[i]]])
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  Protein_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 250,sort.by = 'logFC')
  Protein_diff$logP<- -log10(Protein_diff$adj.P.Val)
  ProteinDEG[[cancernamepr[i]]]<-Protein_diff
  ProteinDEG[[cancernamepr[i]]]$ID<-rownames(ProteinDEG[[cancernamepr[i]]])
}

cancernamepr<-names(ProteinDEG)
DEGtest<-DEG
for (i in 1:24){
  print(i)
  ProteinDEG[[i]]$cancer<-cancernamepr[i]
}

ProteinDEG<-lapply(ProteinDEG,function(df){
  df<-df%>%filter(df$adj.P.Val<0.05)
})

ProteinDEG<-lapply(ProteinDEG,function(df){
  df<-df%>%filter(abs(df$logFC)>0.5)
})

for (i in 1:24){
  print(i)
  ProteinDEG[[i]]<-within(ProteinDEG[[i]],{
    group<-NA
    group[logFC > 0 ] = 'high'
    group[logFC < 0 ] = 'low'
  })
}
ProteinDEG<-lapply(ProteinDEG,function(df){
  df$group<-as.factor(df$group)
  return(df)
})

ProteinDEGbar<-data.frame()
for (i in 1:24){
  print(i)
  ProteinDEGbar[2*i-1,1]<-cancernamepr[i]
  ProteinDEGbar[2*i-1,2]<-summary(ProteinDEG[[i]]$group)[1]
  ProteinDEGbar[2*i-1,3]<-'High'
  ProteinDEGbar[,2][is.na(ProteinDEGbar[,2])] <- 0
  ProteinDEGbar[2*i,1]<-cancernamepr[i]
  ProteinDEGbar[2*i,2]<-summary(ProteinDEG[[i]]$group)[2]
  ProteinDEGbar[2*i,3]<-'Low'
  ProteinDEGbar[,2][is.na(ProteinDEGbar[,2])] <- 0
  ProteinDEGbar[2*i-1,4]<-ProteinDEGbar[2*i-1,2]+ProteinDEGbar[2*i,2]
  ProteinDEGbar[2*i,4]<-ProteinDEGbar[2*i-1,2]+ProteinDEGbar[2*i,2]
}
colnames(ProteinDEGbar)<-c('cancer','value','Group','sum')
DEGbar1<-as.data.frame(t(DEGbar1))
DEGbar1<-DEGbar
DEGbar1<-DEGbar1[match(ProteinDEGbar$cancer,DEGbar1$cancer),]
ProteinDEGbar$rank<-DEGbar1$sum

ProteinDEGbar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% 
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ 
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))