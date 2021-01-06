miRdata<-list()
for(i in 1:length(cancername)){
  print(i)
  miRdata[[cancername[i]]] = read.table(paste("/",cancername[i],'_miR.tsv',sep = ''),
                                        sep = '\t',
                                        header = T,
                                        quote = '',
                                        fill = T, 
                                        comment.char = "!",
                                        stringsAsFactors = FALSE)
}

miRdata<-lapply(miRdata,BIrownames)

miRdata<-lapply(miRdata,function(df){
  colnames(df)<-paste0(colnames(df),"A")
  return(df)
})

miRids<-read.csv("/miRids.csv")
miRids<-miRids[!apply(miRids == "", 1, all),]
miRids<-miRids[miRids$ID %in% colnames(miR_heapmap),]


miRdata<-miRdata_test

miRdata<-lapply(miRdata,function(df){
  df<-df[rownames(df) %in% miRids$CODE,]
  return(df)
})
miRids<-miRids[miRids$CODE %in% rownames(miRdata$ACC),]

miRdata<-lapply(miRdata,function(df){
  df<-df[match(miRids$CODE,rownames(df)),]
  rownames(df)<-miRids$ID
  return(df)
})

miRdata<-lapply(miRdata,function(df){
  df<-t(df)
  df<-as.data.frame(df)
})
miRdata_test<-miRdata

miRdata_test<-lapply(miRdata_test, function(df){
  df$ID<-rownames(df)
  df$ID<-substr(df$ID,9,12)
  return(df)
})

drug_response<-read.csv("/drug_response.csv")
drug_spearman<-read.csv("/drug_spearman.csv")
drug_spearman<-drug_spearman%>%filter(drug_spearman$Spearman.Correlation.All>0.3)
drug_response<-drug_response[drug_response$X %in% drug_spearman$Drug,]
rownames(drug_response)<-drug_response$X
drug_response$X<-NULL

#数据整理与对???
drug_response<-as.data.frame(t(drug_response))
drug_response_list<-list()
drug_response$ID<-rownames(drug_response)
for (i in 1:30){
  print(i)
  drug_response_list[[cancername[i]]]<-drug_response[drug_response$ID %in% miRdata_test[[cancername[i]]]$ID,]
}
for (i in 1:30){
  print(i)
  miRdata_test[[cancername[i]]]<-miRdata_test[[cancername[i]]][miRdata_test[[cancername[i]]]$ID %in% drug_response_list[[cancername[i]]]$ID,]
}
for (i in 1:30){
  print(i)
  miRdata_test[[cancername[i]]]<-miRdata_test[[cancername[i]]][match(drug_response_list[[cancername[i]]]$ID,miRdata_test[[cancername[i]]]$ID),]
}
#去除ID一???
miRdata_test<-lapply(miRdata_test, function(df){
  df<-df[c(1:21)]
  return(df)
})
drug_response_list<-lapply(drug_response_list, function(df){
  df<-df[c(1:138)]
  return(df)
})

#出相关系数矩???
miR_drug_cor<-as.data.frame(matrix(data = NA,nrow = 138,ncol = 21))
miR_drug_cor<-list()

for (i in 5:30){
  print(i)
  if(length(miRdata_test[[cancername[i]]]$`let-7a-2-3p`) != 0)
  {miR_drug_cor[[cancername[i]]]<-as.data.frame(matrix(data = NA,nrow = 138,ncol = 21))
  rownames(miR_drug_cor[[cancername[i]]])<-colnames(drug_response_list[[cancername[i]]])
  colnames(miR_drug_cor[[cancername[i]]])<-colnames(miRdata_test[[cancername[i]]])
  for (j in 1:21){
    print(j)
    for (k in 1:138){
      print(k)
      miR_drug_cor[[cancername[i]]][k,j]<-cor(miRdata_test[[cancername[i]]][j],drug_response_list[[cancername[i]]][k],use = "complete.obs")
    }
  } 
  }else print("sucessful")
}


miR_drug_cor<-miR_drug_cor_test
miR_drug_cor<-lapply(miR_drug_cor,function(df){
  apply(df,2,function(x){
    x<-ifelse(abs(x)<=0.3,NA,x) 
  })
})

miR_drug_cor<-lapply(miR_drug_cor,function(df){
  df<-as.data.frame(t(df))
})

miR_drug_cor_filter<-lapply(miR_drug_cor,function(df){
  df<-apply(df,2,function(x){
    if(summary(is.na(x))[[2]]>=3) print("sucessful") else x<-NA
  }) 
})

miR_drug_cor_filter<-lapply(miR_drug_cor_filter,as.data.frame)
miR_drug_cor_filter<-lapply(miR_drug_cor_filter,function(df){
  df$ID<-rownames(df)
  return(df)
})
miR_drug_cor_filter<-lapply(miR_drug_cor_filter,function(df){
  df=df%>%filter(is.na(df[1])==FALSE)
  return(df)
})