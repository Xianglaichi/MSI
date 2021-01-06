DEG_filter<-lapply(DEG,function(df){
  df<-df[df$adj.P.Val<0.05,]
  return(df)
})

Go_result<-list()
Go_result$ACC<-data.frame()
Go_result$MESO<-data.frame()
for (i in 17:length(cancername)){
  print(i)
  genelist_test<-as.character(DEG_filter[[i]]$ID)
  genelist_test <- bitr(genelist_test, fromType="SYMBOL", toType=c("ENSEMBL"), OrgDb="org.Hs.eg.db")
  genelist_test<-genelist_test[,2]
  print("老巴???")
  Go_resultBP <- enrichGO(genelist_test, 'org.Hs.eg.db', 
                          keyType = "ENSEMBL", 
                          ont="BP", 
                          pvalueCutoff=0.05) 
  result<-Go_resultBP@result
  result<-result[result$p.adjust < 0.05,]
  Go_result[[cancername[i]]]<-result
}

for (i in 1:30){
  print(i)
  Go_result[[i]]$cancer<-cancername[i]
}

Go_summary<-data.frame(Go_result$ACC)
for (i in 1:29){
  print(i)
  Go_summary<-rbind(Go_summary,Go_result[[i+1]])
}
Go_summary$ID<-as.factor(Go_summary$ID)
Go_common<-as.data.frame(summary(Go_summary$ID,maxsum = 100000))
Go_common$ID<-rownames(Go_common)
Go_summary<-Go_summary[match(Go_common$ID,Go_summary$ID),]
Go_common$Description<-Go_summary$Description
Go_common<-Go_common%>%filter(Go_common$`summary(Go_summary$ID, maxsum = 1e+05)`>=10)

radardata <- as.data.frame(matrix(c(23,20,10,16,13,12,25), ncol=7))
colnames(radardata) <- c("Response to hypoxia" , "T cell activation" , "T cell cytotoxicity" , "Antigen processing and presentation" , "Lipid metabolism","Glycoprotein metabolism","EMT transition")
rownames(radardata) <- "GO pathways"
radardata <- rbind(rep(30,7) , rep(0,7) , radardata)
colors_border=c(rgb(0.2,0.4,0.6,0.9))
colors_in=c(rgb(0.2,0.4,0.6,0.3))
radarchart(radardata, axistype=1 , 
           pcol=colors_border, 
           pfcol=colors_in, 
           plwd=4, 
           centerzero = F,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,28,7), cglwd=0.8, 
           vlcex=1.5,
           paxislabels = 4)
