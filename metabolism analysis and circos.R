metabdata<-read.csv("/metabolomics.csv")
rownames(metabdata)<-metabdata$X
metabdata<-metabdata[,-1]
colnames(metabdata)<-paste("TCGA.",colnames(metabdata),".01A",sep = "")

BRCAgroup<-pan_group$BRCA[rownames(pan_group$BRCA) %in% colnames(metabdata),]
metabdata<-metabdata[colnames(metabdata) %in% rownames(BRCAgroup)]


metabdata<-as.matrix(metabdata)
BRCAgroup<-as.matrix(BRCAgroup)
fit <- lmFit(metabdata,BRCAgroup)
contrast.matrix <- makeContrasts(High - Low,levels=BRCAgroup)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
all_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 30000,sort.by = 'logFC')
head(all_diff)
all_diff$'logP'<- -log10(all_diff$'P.Value')

all_diff<-all_diff[match(rownames(metabdata),rownames(all_diff)),]
metabdata<-t(metabdata)
metabdata<-as.data.frame(metabdata)
BRCAscore<-summarylist_rfs$BRCA
BRCAscore<-BRCAscore[match(rownames(metabdata),BRCAscore$id),]
metabdata$score<-BRCAscore$score
all_diff$cor<-NA
for (i in 1:399){
  print(i)
  all_diff[i,9]<-cor(metabdata[i],metabdata$score)
}
all_diff$cor_abs<-abs(all_diff$cor)

all_diff<-within(all_diff,{
  cor_class<-NA
  cor_class[cor < 0] = 'Negative'
  cor_class[cor > 0] = 'Positive'
})


ggscatter(all_diff,x="logFC",y="logP",size = "cor_abs")+theme_base()
all_diff$Group = "not-significant"
all_diff$Group[which((all_diff$P.Val<0.30) & (all_diff$logFC> 0.5))] = "up-regulated in High-score patients"
all_diff$Group[which((all_diff$P.Val<0.30) & (all_diff$logFC< -0.5))] = "up-regulated in Low-score patients"
table(all_diff$Group)
all_diff$Group<-as.factor(all_diff$Group)

ggscatter(all_diff,x="logFC",y="logP",size = "cor_abs",
          color = "Group",
          palette = c('gray',"#FF3333","#00CCFF")
)+theme_base()



names(all_diff)[10]<-"Correlation.index"
names(all_diff)[11]<-"Correlation"

ggscatter(all_diff,x="logFC",y="logP",size = "Correlation.index",
          color = "Group",
          palette = c('gray',"#FF3333","#0099FF"),alpha = 0.6,
          shape = "Correlation",
          label = all_diff$Label,
          font.label = 10,
          repel = T)+theme_base()+
  geom_hline(yintercept = 0.5228787,linetype = "dashed") + 
  geom_vline(xintercept = c(-0.5,0.5),linetype = "dashed")

metabdata<-read.csv("/metabolomics.csv")
rownames(metabdata)<-metabdata$X
metabdata<-metabdata[,-1]
colnames(metabdata)<-paste("TCGA.",colnames(metabdata),".01A",sep = "")
BRCAgroup<-pan_group$BRCA[rownames(pan_group$BRCA) %in% colnames(metabdata),]
metabdata<-metabdata[colnames(metabdata) %in% rownames(BRCAgroup)]
metabdata<-as.matrix(metabdata)
BRCAgroup<-as.matrix(BRCAgroup)
fit <- lmFit(metabdata,BRCAgroup)
contrast.matrix <- makeContrasts(High - Low,levels=BRCAgroup)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
all_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 30000,sort.by = 'logFC')
head(all_diff)
all_diff$'logP'<- -log10(all_diff$'P.Value')
all_diff<-all_diff[match(rownames(metabdata),rownames(all_diff)),]
metabdata<-t(metabdata)
metabdata<-as.data.frame(metabdata)
BRCAscore<-summarylist_rfs$BRCA
BRCAscore<-BRCAscore[match(rownames(metabdata),BRCAscore$id),]
metabdata$score<-BRCAscore$score
all_diff$cor<-NA
all_diff$cor_p<-NA
for (i in 1:399){
  print(i)
  all_diff[i,9]<-cor.test(metabdata[,i],metabdata$score)$estimate
  all_diff[i,10]<-cor.test(metabdata[,i],metabdata$score)$p.value
}
all_diff$cor_abs<-abs(all_diff$cor)

all_diff<-within(all_diff,{
  cor_class<-NA
  cor_class[cor < 0] = 'Negative'
  cor_class[cor > 0] = 'Positive'
})
all_diff$ID<-rownames(all_diff)
all_diff<-all_diff[,-c(2:8)]
metab_diff_posi<-all_diff%>%filter((cor>0) & (cor_p<0.05))
metab_diff_posi<-metab_diff_posi[order(-metab_diff_posi$cor),]

metab_diff_nega<-all_diff%>%filter((cor<0) & (cor_p<0.05))
metab_diff_nega<-metab_diff_nega[order(metab_diff_nega$cor),]
metab_diff_nega<-metab_diff_nega[c(1:10),]
metab_diff_nega<-metab_diff_nega[order(-metab_diff_nega$cor),]
library(circlize)
circos.clear()
category = metab_diff_nega$ID 
percent = metab_diff_nega$cor*(-100)
color = rev(rainbow(length(percent))) 

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100)) 
circos.track(ylim = c(0.5, length(percent)+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], 10), 1:10,
                               rep(xlim[2], 10), 1:10,
                               col = "#CCCCCC")
               circos.rect(rep(0, 10), 1:10 - 0.45, percent, 1:10 + 0.45,
                           col = color, border = "white")
               circos.text(rep(xlim[1], 10), 1:10, 
                           paste0(category, ""), 
                           facing = "downward", adj = c(1.05, 0.5), cex = 0.8) 
               breaks = seq(0, 85, by = 5)
               circos.axis(h = "top", major.at = breaks, labels = paste0("0.",breaks), 
                           labels.cex = 0.6)
             })
circos.clear()

library(circlize)
circos.clear()
category = metab_diff_posi$ID 
percent = metab_diff_posi$cor*(100)
color = rev(rainbow(length(percent))) 

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100)) 
circos.track(ylim = c(0.5, length(percent)+0.5), track.height = 0.8, 
             bg.border = NA, panel.fun = function(x, y) {
               xlim = CELL_META$xlim
               circos.segments(rep(xlim[1], 2), 1:2,
                               rep(xlim[2], 2), 1:2,
                               col = "#CCCCCC")
               circos.rect(rep(0, 2), 1:2 - 0.45, percent, 1:2 + 0.45,
                           col = color, border = "white")
               circos.text(rep(xlim[1], 2), 1:2, 
                           paste0(category, ""), 
                           facing = "downward", adj = c(1.05, 0.5), cex = 0.8) 
               breaks = seq(0, 85, by = 5)
               circos.axis(h = "top", major.at = breaks, labels = paste0("0.",breaks), 
                           labels.cex = 0.6)
             })
circos.clear()
