library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
CNVdata<-list()
for (i in 1:length(cancername)){
  print(i)
  CNVdata[[cancername[i]]]<-read.table(paste('/TCGA-',cancername[i],'.gistic.tsv',sep = ''),
                                       sep = '\t',
                                       header = T,
                                       quote = '',
                                       fill = T, 
                                       stringsAsFactors = FALSE)
  CNVdata[[cancername[i]]]$Gene.Symbol<-gsub("\\..*","",CNVdata[[cancername[i]]]$Gene.Symbol)
  genename <- as.character(CNVdata[[cancername[i]]]$Gene.Symbol) 
  genename <- bitr(genename, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
  genename<-genename[match(CNVdata[[cancername[i]]]$Gene.Symbol,genename$ENSEMBL),]
  CNVdata[[cancername[i]]]$Gene.Symbol<-genename$SYMBOL
  CNVdata[[cancername[i]]]%>%filter(is.na(Gene.Symbol)==FALSE)
  genename=genename%>%filter(is.na(SYMBOL)==FALSE)
  CNVdata[[cancername[i]]]<-as.matrix(CNVdata[[cancername[i]]])
  rownames(CNVdata[[cancername[i]]])<-CNVdata[[cancername[i]]][,1]
  CNVdata[[cancername[i]]]<-CNVdata[[cancername[i]]][,-1]
  CNVdata[[cancername[i]]]<-t(CNVdata[[cancername[i]]])
  CNVdata[[cancername[i]]]<-as.data.frame(CNVdata[[cancername[i]]])
  CNVdata[[cancername[i]]]$id<-rownames(CNVdata[[cancername[i]]])
  CNVdata[[cancername[i]]]<-melt(CNVdata[[cancername[i]]],id.vars=c("id"),variable.name="gene",value.name="value")
  CNVdata[[cancername[i]]]$value<-as.character(CNVdata[[cancername[i]]]$value)
  CNVdata[[cancername[i]]]$value<-as.numeric(CNVdata[[cancername[i]]]$value)
  CNVdata[[cancername[i]]] = CNVdata[[cancername[i]]]%>%filter(value != 0)
  CNVdata[[cancername[i]]]$CN<-NA
  CNVdata[[cancername[i]]]$CN[which(CNVdata[[cancername[i]]]$value ==1)] = "Amp"
  CNVdata[[cancername[i]]]$CN[which(CNVdata[[cancername[i]]]$value ==-1)] = "Del"
  CNVdata[[cancername[i]]]<-CNVdata[[cancername[i]]][,c("gene","id","CN")]
  colnames(CNVdata[[cancername[i]]])<-c("Gene","Sample_name","CN")
}
CNVdatatest<-CNVdata

MAFdata<-list()
MAFdatahigh<-list()
MAFdatalow<-list()
for(i in 1:length(cancername)){
  print(i)
  MAFdata[[cancername[i]]] = read.maf(paste("/",cancername[i],'.maf',sep = ''))
  datas<-MAFdata[[cancername[i]]]@data
  datas<-as.data.frame(datas)
  datas$Tumor_Sample_Barcode<-gsub('-','.',datas$Tumor_Sample_Barcode)
  datas$Tumor_Sample_Barcode<-gsub('............$','',datas$Tumor_Sample_Barcode)
  MAFdatahigh[[cancername[i]]] = datas[datas$Tumor_Sample_Barcode %in% rownames(pan_group[[i]])[which(pan_group[[i]]$High == 1)],]
  MAFdatahigh[[cancername[i]]] = read.maf(MAFdatahigh[[cancername[i]]])
  MAFdatalow[[cancername[i]]] = datas[datas$Tumor_Sample_Barcode %in% rownames(pan_group[[i]])[which(pan_group[[i]]$High == 0)],]
  MAFdatalow[[cancername[i]]] = read.maf(MAFdatalow[[cancername[i]]])
}

CNVdatahigh<-list()
CNVdatalow<-list()
MAFdatahigh_cn<-list()
MAFdatalow_cn<-list()
datahigh<-data.frame()
datalow<-data.frame()

for (i in 2:30){
  print(i)
  datahigh<-MAFdatahigh[[cancername[i]]]@data
  datahigh<-as.data.frame(datahigh)
  datalow<-MAFdatalow[[cancername[i]]]@data
  datalow<-as.data.frame(datalow)
  CNVdatahigh[[cancername[i]]] = CNVdata[[cancername[i]]][CNVdata[[cancername[i]]]$Sample_name %in% datahigh$Tumor_Sample_Barcode,]
  CNVdatahigh[[cancername[i]]] = CNVdatahigh[[cancername[i]]][CNVdatahigh[[cancername[i]]]$Gene %in% datahigh$Hugo_Symbol,]
  CNVdatalow[[cancername[i]]] = CNVdata[[cancername[i]]][CNVdata[[cancername[i]]]$Sample_name %in% datalow$Tumor_Sample_Barcode,]
  CNVdatalow[[cancername[i]]] = CNVdatalow[[cancername[i]]][CNVdatalow[[cancername[i]]]$Gene %in% datalow$Hugo_Symbol,]
  MAFdatahigh_cn[[cancername[i]]]<-read.maf(datahigh,cnTable = CNVdatahigh[[cancername[i]]])
  MAFdatalow_cn[[cancername[i]]]<-read.maf(datalow,cnTable = CNVdatalow[[cancername[i]]])
}


MAFdata_compare_cn<-list()
for (i in 1:length(cancername)){
  print(i)
  subsethigh<-subsetMaf(MAFdatahigh_cn[[cancername[i]]], tsb = NULL, genes = NULL, fields = NULL,
                        query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE,
                        dropLevels = TRUE, restrictTo = "cnv")
  subsetlow<-subsetMaf(MAFdatalow_cn[[cancername[i]]], tsb = NULL, genes = NULL, fields = NULL,
                       query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE,
                       dropLevels = TRUE, restrictTo = "cnv")
  MAFdata_compare_cn[[cancername[i]]] <- mafCompare(m1 = subsethigh, m2 = subsetlow, 
                                                    m1Name = 'High', m2Name = 'Low', minMut = 5,useCNV = TRUE)#最小mut数为3的基???
}


MAFdata_compare_cn_test<-MAFdata_compare_cn

MAFdata_compare_cn_summary<-lapply(MAFdata_compare_cn,function(df){
  df<-as.data.frame(df$results)
  return(df)
})


MAFdata_compare_cn_summary<-lapply(MAFdata_compare_cn_summary,function(df){
  df<-df%>%filter(df$adjPval<0.05)
})


for (i in 1:30){
  print(i)
  MAFdata_compare_cn_summary[[i]]<-within(MAFdata_compare_cn_summary[[i]],{
    group<-NA
    group[or > 1] = 'high'
    group[or < 1] = 'low'
  })
}
MAFdata_compare_cn_summary<-lapply(MAFdata_compare_cn_summary,function(df){
  df$group<-as.factor(df$group)
  return(df)
})

MAFdata_compare_cn_summarybar<-data.frame()
for (i in 1:30){
  print(i)
  MAFdata_compare_cn_summarybar[2*i-1,1]<-cancername[i]
  MAFdata_compare_cn_summarybar[2*i-1,2]<-summary(MAFdata_compare_cn_summary[[i]]$group)[1]
  MAFdata_compare_cn_summarybar[2*i-1,3]<-'High'
  MAFdata_compare_cn_summarybar[,2][is.na(MAFdata_compare_cn_summarybar[,2])] <- 0
  MAFdata_compare_cn_summarybar[2*i,1]<-cancername[i]
  MAFdata_compare_cn_summarybar[2*i,2]<-summary(MAFdata_compare_cn_summary[[i]]$group)[2]
  MAFdata_compare_cn_summarybar[2*i,3]<-'Low'
  MAFdata_compare_cn_summarybar[,2][is.na(MAFdata_compare_cn_summarybar[,2])] <- 0
  MAFdata_compare_cn_summarybar[2*i-1,4]<-MAFdata_compare_cn_summarybar[2*i-1,2]+MAFdata_compare_cn_summarybar[2*i,2]
  MAFdata_compare_cn_summarybar[2*i,4]<-MAFdata_compare_cn_summarybar[2*i-1,2]+MAFdata_compare_cn_summarybar[2*i,2]
}
colnames(MAFdata_compare_cn_summarybar)<-c('cancer','value','Group','sum')
MAFdata_compare_cn_summarybar$rank<-as.numeric(DEGbar$sum)

MAFdata_compare_cn_summarybar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% 
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ 
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))

