ICI<-as.data.frame(matrix(nrow = 8,ncol = 31))
colnames(ICI)<-c(cancername,'group')
ICI$group<-c("CTLA4-high","CTLA4-low",'PDCD1-high','PDCD1-low','LAG3-high','LAG3-low','CD274-high','CD274-low')

datalist_new<-lapply(datalist_new,function(df){
  df<-t(df)
})
datalist_new<-lapply(datalist_new,as.data.frame)
datalist_ici<-lapply(datalist_new,function(df){
  df[,c('CTLA4','PDCD1','LAG3','CD274')]
})
#分组
for (i in 1:30){
  datalist_ici[[i]]<-datalist_ici[[i]][match(rownames(summarylist_rfs[[i]]),rownames(datalist_ici[[i]])),]
  datalist_ici[[i]]$group<-summarylist_rfs[[i]]$score_group
}

for (i in 1:30){
  print(i)
  ICI[1,i]<-tapply(datalist_ici[[i]]$CTLA4,datalist_ici[[i]]$group,mean)[1]
  ICI[2,i]<-tapply(datalist_ici[[i]]$CTLA4,datalist_ici[[i]]$group,mean)[2]
  ICI[3,i]<-tapply(datalist_ici[[i]]$PDCD1,datalist_ici[[i]]$group,mean)[1]
  ICI[4,i]<-tapply(datalist_ici[[i]]$PDCD1,datalist_ici[[i]]$group,mean)[2]
  ICI[5,i]<-tapply(datalist_ici[[i]]$LAG3,datalist_ici[[i]]$group,mean)[1]
  ICI[6,i]<-tapply(datalist_ici[[i]]$LAG3,datalist_ici[[i]]$group,mean)[2]
  ICI[7,i]<-tapply(datalist_ici[[i]]$CD274,datalist_ici[[i]]$group,mean)[1]
  ICI[8,i]<-tapply(datalist_ici[[i]]$CD274,datalist_ici[[i]]$group,mean)[2]
}

for (i in 1:30){
  print(i)
  ICI[1,i]<-tapply(datalist_ici[[i]]$CTLA4,datalist_ici[[i]]$group,mean)[1]
  ICI[2,i]<-tapply(datalist_ici[[i]]$CTLA4,datalist_ici[[i]]$group,mean)[2]
  ICI[3,i]<-tapply(datalist_ici[[i]]$PDCD1,datalist_ici[[i]]$group,mean)[1]
  ICI[4,i]<-tapply(datalist_ici[[i]]$PDCD1,datalist_ici[[i]]$group,mean)[2]
  ICI[5,i]<-tapply(datalist_ici[[i]]$LAG3,datalist_ici[[i]]$group,mean)[1]
  ICI[6,i]<-tapply(datalist_ici[[i]]$LAG3,datalist_ici[[i]]$group,mean)[2]
  ICI[7,i]<-tapply(datalist_ici[[i]]$CD274,datalist_ici[[i]]$group,mean)[1]
  ICI[8,i]<-tapply(datalist_ici[[i]]$CD274,datalist_ici[[i]]$group,mean)[2]
}