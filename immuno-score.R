mmune_cor<-data.frame(rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30),rep(NA,30))
colnames(immune_cor)<-c('cancer','Immunescore_cor','Immunescore_pvalue','CD8_cor','CD8_pvalue',
                        'IFN_cor','IFN_pvalue','TILs_cor','TILs_pvalue','CYT_cor','CYT_pvalue','Hypoxia_cor','Hypoxia_pvalue')
for (i in 1:30){
  print(i)
  immune_cor$cancer[i]<-cancername[i]
  immune_cor$Immunescore_cor[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$ImmuneScore,gsvalist_new[[i]]$up.15,)$estimate
  immune_cor$Immunescore_pvalue[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$ImmuneScore,gsvalist_new[[i]]$up.15)$p.value
  immune_cor$CD8_cor[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$CD8,gsvalist_new[[i]]$up.15)$estimate
  immune_cor$CD8_pvalue[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$CD8,gsvalist_new[[i]]$up.15)$p.value
  immune_cor$IFN_cor[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$IFN,gsvalist_new[[i]]$up.15)$estimate
  immune_cor$IFN_pvalue[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$IFN,gsvalist_new[[i]]$up.15)$p.value
  immune_cor$TILs_cor[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$TILs,gsvalist_new[[i]]$up.15)$estimate
  immune_cor$TILs_pvalue[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$TILs,gsvalist_new[[i]]$up.15)$p.value
  immune_cor$CYT_cor[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$CYT,gsvalist_new[[i]]$up.15)$estimate
  immune_cor$CYT_pvalue[i]<-cor.test(ESTIMATElist_ssgsea[[i]]$CYT,gsvalist_new[[i]]$up.15)$p.value
  immune_cor$Hypoxia_cor[i]<-cor.test(Hypoxia_ssgsea[[i]]$Hypoxia,gsvalist_new[[i]]$up.15)$estimate
  immune_cor$Hypoxia_pvalue[i]<-cor.test(Hypoxia_ssgsea[[i]]$Hypoxia,gsvalist_new[[i]]$up.15)$p.value
  
}
immune_cor_melt<-read.csv("/immune_cor.csv")
immune_cor_melt<-within(immune_cor_melt,{
  sign<-NA
  sign[Correlation > 0] = 'Positive'
  sign[Correlation < 0] = 'Negative'  
})

immune_cor_melt$'-logP_value'<--log10(immune_cor_melt$p_value)

immune_cor_melt$variable<-factor(immune_cor_melt$variable,levels = c('Hypoxia','APM','Mutation_load',"CD8","CYT","IFN","TILs","Immunescore"))#这样来改变显示在图中的顺???
ggplot(immune_cor_melt)+geom_point(aes(x=cancer,y=variable,size=`-logP_value`,color = Correlation,alpha = 0.7))+
  scale_color_gradient(low = "#0099FF",high = "red")+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5))
