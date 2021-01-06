GDSC_gene_drug_cor<-GDSC_gene_drug_cor_test
GDSC_gene_drug_cor$gene<-rownames(GDSC_gene_drug_cor)
Mispre<-melt(GDSC_gene_drug_cor,id = "gene",value.name = "value")
Mispre<-Mispre%>%filter(abs(value)>=0.3)
Mispre<-Mispre[Mispre$gene %in% hx$ID,]
Mispre<-within(Mispre,{
  group<-NA
  group[value >= 0] = "positive"
  group[value <= 0] = "negative"
})

GDSCtarget<-read.csv("/GDSCdrugtarget_filter.csv")
GDSCtarget<-GDSCtarget[,c(1,4)]
GDSCtarget<-GDSCtarget_test
length(summary(GDSCtarget$target_pathway,maxsum = 200))
length(summary(GDSCtarget$drug_name,maxsum = 200))
length(summary(Mispre$variable,maxsum = 200))

GDSCtarget<-GDSCtarget[match(Mispre$variable,GDSCtarget$drug_name),]
Mispre<-cbind(Mispre,GDSCtarget)

Mispre_origin<-read.csv("/Mispre.csv")
Mispre_origin<-Mispre_origin[Mispre_origin$gene %in% Mispre$gene,]
Mispre<-Mispre_origin
length(summary(Mispre$variable,maxsum = 200))
length(summary(Mispre$gene,maxsum = 400))
length(summary(Mispre$target_pathway,maxsum = 400))
Mispre$target_pathway<-as.factor(Mispre$target_pathway)

hx<-data.frame(rownames(GDSC_gene_drug_cor))
names(hx)<-"name"
hx$group<-rep("gene",685)
hx<-hx[hx$name %in% Mispre$gene,]

kakaka<-data.frame(colnames(GDSC_gene_drug_cor))
names(kakaka)<-"name"
kakaka$group<-rep("drug",168)
kakaka<-kakaka[kakaka$name %in% Mispre$variable,]

MisLinks<-rbind(hx,kakaka)
MisLinks$id<-0:268
for (i in 1:5){
  print(i)
  MisLinks<-rbind(MisLinks,MisLinks)
}
MisLinks_source<-MisLinks[match(Mispre$gene,MisLinks$name),]
MisLinks_target<-MisLinks[match(Mispre$variable,MisLinks$name),]
MisLinks<-data.frame(source=MisLinks_source$id,target=MisLinks_target$id,value=Mispre$value,color=Mispre$group)
MisLinks$color<-ifelse(MisLinks$color == "positive", "red", "blue")
MisLinks$value<-MisLinks$value*10
MisLinks<-MisLinks[order(MisLinks$source),]
MisLinks$source<-as.numeric(MisLinks$source)
MisLinks$target<-as.numeric(MisLinks$target)

#MisNodes
hx$group<-rep("gene",198)
kakaka$group<-rep("drug",71)
Mispre_kakaka<-Mispre[match(kakaka$name,Mispre$variable),]
kakaka$group<-Mispre_kakaka$target_pathway

MisNodes<-rbind(hx,kakaka)
MisNodes$size<-10
Html <- forceNetwork(Links = MisLinks, 
                     Nodes = MisNodes, 
                     Source = "source",
                     Target = "target", 
                     Value = "value",
                     NodeID = "name",
                     Group = "group",
                     Nodesize = "size" ,
                     fontSize = 20, 
                     linkColour=MisLinks$color, 
                     charge = -200,）    
                     opacity = 0.9,  
                     legend=F, 
                     arrows=F,
                     bounded=F,
                     opacityNoHover=TRUE, 
                     zoom = T
) 