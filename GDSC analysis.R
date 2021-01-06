GDSC_IC50<-read.csv('/GDSC_IC50.csv')
GDSC_IC50$Cosmic.sample.Id<-paste0("DATA.",GDSC_IC50$Cosmic.sample.Id)
GDSC_IC50<-GDSC_IC50[GDSC_IC50$TCGA.classification %in% cancername,]
#GDSC表达数据
GDSCexp<-read.table('/GDSCexpdata.txt',
                    sep = '\t',
                    header = T,
                    quote = '',
                    fill = T, 
                    stringsAsFactors = FALSE)

GDSCexp <- GDSCexp %>% 
  select(GENE_SYMBOLS, everything()) %>% 
  mutate(rowMean =rowMeans(.[grep("DATA", names(.))])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(GENE_SYMBOLS,.keep_all = T) %>% 
  select(-rowMean) 

GDSCexp<-read.csv('/GDSCexp整理.csv')

rownames(GDSCexp)<-GDSCexp$GENE_SYMBOLS
GDSCexp<-GDSCexp[,-(1:2)]

GDSCexp<-GDSCexp[,colnames(GDSCexp) %in% GDSC_IC50$Cosmic.sample.Id]
GDSC_IC50<-GDSC_IC50[GDSC_IC50$Cosmic.sample.Id %in% colnames(GDSCexp),] 


GDSCexp<-as.data.frame(t(GDSCexp))
GDSCsummary<-GDSC_IC50[,c(1,4,9)] 
GDSCsummary<-dcast(GDSCsummary,Cosmic.sample.Id ~ Drug.name, 
                   value.var = "AUC",fun.aggregate = mean) 
rownames(GDSCsummary)<-GDSCsummary$Cosmic.sample.Id
GDSCsummary<-GDSCsummary[,-1]
GDSCexp<-GDSCexp[match(rownames(GDSCsummary),rownames(GDSCexp)),]
GDSCexp<-as.matrix(GDSCexp)
GDSCexp<-t(GDSCexp)

GDSCgsva <- gsva(GDSCexp,invasive_signature, 
                 min.sz = 1,
                 max.sz = Inf,
                 mx.diff=TRUE, 
                 verbose=FALSE, 
                 parallel.sz=0,
                 method = "gsva")
GDSCgsva<-as.data.frame(t(GDSCgsva))
GDSCgsva<-GDSCgsva[match(rownames(GDSCsummary),rownames(GDSCgsva)),]


GDSCsummary$score<-GDSCgsva$MSigDB.24

GDSC_IC50<-GDSC_IC50[match(rownames(GDSCsummary),GDSC_IC50$Cosmic.sample.Id),]
GDSCsummary$cancer<-GDSC_IC50$TCGA.classification

GDSCdata<-list()
for (i in 1:length(cancername)){
  print(i)
  GDSCdata[[cancername[i]]]<-GDSCsummary[GDSCsummary$cancer %in% cancername[i],]
}

GDSCcor_all<-as.data.frame(colnames(GDSCsummary))
GDSCcor_all$cor<-NA
GDSCcor_all$p.value<-NA
GDSCcor_all<-GDSCcor_all[-c(168:169),]
for (i in 1:167){
  print(i)
  hahaha<-cor.test(GDSCsummary[,i],GDSCsummary[,168],na.action = "na.omit")
  GDSCcor_all[i,2]<-hahaha$estimate
  GDSCcor_all[i,3]<-hahaha$p.value
}
summary(GDSCcor_all)
GDSCcor_all<-GDSCcor_all%>%filter(GDSCcor_all$p.value<0.05)
GDSCcor_all<-GDSCcor_all%>%filter(abs(GDSCcor_all$cor)>0.22)


GDSC_IC50$Cosmic.sample.Id<-as.factor(GDSC_IC50$Cosmic.sample.Id)
test<-as.data.frame(summary(GDSC_IC50$Cosmic.sample.Id,maxsum = 10000))
GDSC_IC50$Drug.name<-as.factor(GDSC_IC50$Drug.name)
test<-as.data.frame(summary(GDSC_IC50$Drug.name,maxsum = 10000))


ggplot(GDSCsummary, aes(x=GDSCsummary$`KU-55933`, y=GDSCsummary$score)) +
  geom_point(color = "#0000FF",shape = 20, alpha = 0.6,size =5)+
  geom_smooth(method=lm, se=TRUE, linetype= 1,
              color="#FF3333")+
  theme_light()+
  stat_cor(data=GDSCsummary, method = "pearson",size = 12)+
  xlab("KU55933")+ylab("MSI score")+
  theme(axis.title.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5))

ggplot(GDSCsummary, aes(x=GDSCsummary$SB216763, y=GDSCsummary$score)) +
  geom_point(color = "#0000FF",shape = 20, alpha = 0.6,size =5)+
  geom_smooth(method=lm, se=TRUE, linetype= 1,
              color="#FF3333")+
  theme_light()+
  stat_cor(data=GDSCsummary, method = "pearson",size = 12,na.action = "na.omit")+
  xlab("SB216763")+ylab("MSI score")+
  theme(axis.title.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5))

ggplot(GDSCsummary, aes(x=GDSCsummary$NU7441, y=GDSCsummary$score)) +
  geom_point(color = "#0000FF",shape = 20, alpha = 0.6,size =5)+
  geom_smooth(method=lm, se=TRUE, linetype= 1,
              color="#FF3333")+
  theme_light()+
  stat_cor(data=GDSCsummary, method = "pearson",size = 12,na.action = "na.omit")+
  xlab("NU7441")+ylab("MSI score")+
  theme(axis.title.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5))