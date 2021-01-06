Methydata<-data.frame()
Methy_group<-lapply(summarylist_rfs,function(df){
  df<-df[c(1,8)]
})
myLoad<-list()
diff_Methy<-list()

for(i in 1:length(cancername)){
  print(i)
  print(cancername[i])
  Methydata = read.table(paste("/TCGA-",cancername[i],'.methylation450.tsv',sep = ''),
                         sep = '\t',
                         header = T,
                         quote = '',
                         fill = T, 
                         comment.char = "!",
                         stringsAsFactors = FALSE)
  rownames(Methydata)<-Methydata[,1]
  Methydata<-Methydata[,-1]
  Methydata<-Methydata[,colnames(Methydata) %in% rownames(Methy_group[[cancername[i]]])]
  Methy_group[[cancername[i]]]<-Methy_group[[cancername[i]]][rownames(Methy_group[[cancername[i]]]) %in% colnames(Methydata),]
  Methydata<-Methydata[,match(rownames(Methy_group[[cancername[i]]]),colnames(Methydata))]
  Methydata<-as.matrix(Methydata)
  myLoad$beta<-Methydata
  myLoad$pd<-Methy_group[[cancername[i]]]
  print("successful") 
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=4)
  myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$score_group,adjPVal = 0.05,adjust.method = "BH",arraytype = "450K")
  diff_Methy[[cancername[i]]]<-myDMP
}

Methydata<-Methydata1
Methy_group<-Methy_group1
for (i in 1:17) {
  print(i) 
  Methydata[[i]]<-lapply(Methydata[[i]],BIrownames)
}  


for (i in 1:30){
  print(i)
  Methydata[[i]]<-Methydata[[i]][,colnames(Methydata[[i]]) %in% rownames(Methy_group[[i]])]
}
for (i in 1:30){
  print(i)
  Methy_group[[i]]<-Methy_group[[i]][rownames(Methy_group[[i]]) %in% colnames(Methydata[[i]]),]
}
for (i in 1:30){
  print(i)
  Methydata[[i]]<-Methydata[[i]][,match(rownames(Methy_group[[i]]),colnames(Methydata[[i]]))]
}

Methydata<-lapply(Methydata,as.matrix)
Methy_group<-lapply(Methy_group,as.matrix)
Methy_group<-lapply(Methy_group,function(df){
  df[,1]<-as.numeric(df[,1])
  df[,2]<-as.numeric(df[,2])
  return(df)
})

Methydata<-lapply(Methydata,function(df){
  rownames(df)<-Methyids$gene
})


MethyDEG<-list()
Methy_diff<-data.frame()
for (i in 1:17){
  print(i)
  fit <- lmFit(Methydata[[i]],Methy_group[[i]])
  contrast.matrix <- makeContrasts(High - Low,levels=Methy_group[[i]])
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  Methy_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 1881,sort.by = 'logFC')
  Methy_diff$logP<- -log10(Methy_diff$adj.P.Val)
  MethyDEG[[cancername[i]]]<-Methy_diff
}
DEGhahaha<-MethyDEG

for (i in 1:17){
  print(i)
  MethyDEG[[i]]$cancer<-cancername[i]
}

MethyDEG_test<-MethyDEG
MethyDEG<-MethyDEG_test
MethyDEG<-lapply(MethyDEG,function(df){
  df<-df%>%filter(df$adj.P.Val<0.05)
})

MethyDEG<-lapply(MethyDEG,function(df){
  df<-df%>%filter(abs(df$logFC)>0.2)
})

for (i in 1:17){
  print(i)
  MethyDEG[[i]]<-within(MethyDEG[[i]],{
    group<-NA
    group[logFC > 0 ] = 'high'
    group[logFC < 0 ] = 'low'
  })
}
MethyDEG<-lapply(MethyDEG,function(df){
  df$group<-as.factor(df$group)
  return(df)
})


MethyDEGbar<-data.frame()
for (i in 1:17){
  print(i)
  MethyDEGbar[2*i-1,1]<-cancername[i]
  MethyDEGbar[2*i-1,2]<-summary(MethyDEG[[i]]$group)[1]
  MethyDEGbar[2*i-1,3]<-'High'
  MethyDEGbar[,2][is.na(MethyDEGbar[,2])] <- 0
  MethyDEGbar[2*i,1]<-cancername[i]
  MethyDEGbar[2*i,2]<-summary(MethyDEG[[i]]$group)[2]
  MethyDEGbar[2*i,3]<-'Low'
  MethyDEGbar[,2][is.na(MethyDEGbar[,2])] <- 0
  MethyDEGbar[2*i-1,4]<-MethyDEGbar[2*i-1,2]+MethyDEGbar[2*i,2]
  MethyDEGbar[2*i,4]<-MethyDEGbar[2*i-1,2]+MethyDEGbar[2*i,2]
}
colnames(MethyDEGbar)<-c('cancer','value','Group','sum')
MethyDEGbar$rank<-DEGbar$sum

MethyDEGbar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% 
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ 
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))
