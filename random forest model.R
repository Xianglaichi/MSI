gsvalist<-read.csv("/pangsva.csv")
rownames(gsvalist)<-gsvalist$X
gsvalist$X<-NULL
gsvalist<-t(gsvalist)
gsvalist<-as.data.frame(gsvalist)

pangroup<- read.csv("/pangroup.csv")
rownames(pangroup)<-pangroup[,1]
pangroup<-pangroup[,-1]
pangroup<-pangroup[match(rownames(gsvalist),rownames(pangroup)),]

gsvalist$group<-pangroup$group
colnames<-colnames(gsvalist)
colnames<-gsub(' ','.',colnames)
colnames(gsvalist)<-colnames

train <- sample(1:nrow(gsvalist), round(nrow(gsvalist) * 0.70))
gsvalist_train<-gsvalist[train,]
gsvalist_test<-gsvalist[-train,]

n<-length(names(gsvalist_train))     
rate=1     
for(i in 1:(n-1)){
  set.seed(100)
  rf_train<-randomForest(as.factor(gsvalist_train$group)~.,data=gsvalist_train,mtry=i,ntree=1000)
  rate[i]<-mean(rf_train$err.rate)  
  print(rf_train)    
}
rate     
plot(rate)


rf_train<-randomForest(as.factor(gsvalist_train$group)~.,data=gsvalist_train
                       ,mtry=5,importance=TRUE,
                       na.action=na.omit,ntree = 5000)
plot(rf_train)
rf_train<-randomForest(as.factor(gsvalist_train$group)~.,data=gsvalist_train,
                       mtry=5,ntree=5000,
                       importance=TRUE,proximity=TRUE,
                       na.action = na.omit)    
rf_train
summary(rf_train)
plot(rf_train)
saveRDS(rf_train, file = "/model_rf.RDS")
rf_train<-readRDS("/model_rf.RDS")
importance(rf_train)
rf_TCGA_importance<-as.data.frame(rf_TCGA_importance)
write.csv(rf_TCGA_importance,file = "/rf_TCGA_importance.csv")
MDSplot(rf_train,gsvalist_train$group,palette=c(groupA="#00CCFF",groupB="#FF3333"),pch=as.numeric(gsvalist_train$group))
pred1<-predict(rf_train,gsvalist_train)
Freq1<-table(pred1,gsvalist_train$group)
sum(diag(Freq1))/sum(Freq1)

pred2<-predict(rf_train,gsvalist_test)
Freq2<-table(pred2,gsvalist_test$group)
sum(diag(Freq2))/sum(Freq2)

class(rf_train)
