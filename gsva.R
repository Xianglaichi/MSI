TCGAids = read.table("/TCGAids.txt",
                     sep = '\t',
                     header = T,
                     quote = '',
                     fill = T, 
                     comment.char = "!",
                     stringsAsFactors = FALSE)
TCGAids<-TCGAids[c(1:2)]
{
  LUAD
  LUADdata = read.table("/TCGA-LUAD.htseq_fpkm.tsv",
                        sep = '\t',
                        header = T,
                        quote = '',
                        fill = T, 
                        comment.char = "!",
                        stringsAsFactors = FALSE)
  names(LUADdata)[1] <- names(TCGAids)[1]
  LUADdata$id <- as.character(LUADdata$id)
  TCGAids<-TCGAids[match(LUADdata$id,TCGAids$id),]
  LUADdata <- LUADdata %>% 
    inner_join(TCGAids,by="id") %>% 
    select(-id) %>% 
    select(gene, everything()) %>%
    mutate(rowMean =rowMeans(.[grep("TCGA", names(.))])) %>% 
    arrange(desc(rowMean)) %>% 
    distinct(gene,.keep_all = T) %>% 
    select(-rowMean) %>% 
    tibble::column_to_rownames(colnames(.)[1])
  LUADdata<-select(LUADdata,ends_with("01A"))
  LUADpheno<-read.table("/TCGA-LUAD.GDC_phenotype.tsv",
                        sep = '\t',
                        header = T,
                        quote = '',
                        fill = T, 
                        comment.char = "!",
                        stringsAsFactors = FALSE)
  LUADpheno$submitter_id.samples<-gsub('-','.',LUADpheno[,1])
  names(LUADpheno)[1]<-"id"
  LUADpheno<-LUADpheno[match(colnames(LUADdata),LUADpheno$id),]  
  LUADpheno<-LUADpheno%>%filter(is.na(LUADpheno$id) == FALSE)
  LUADsurvival<-read.table("/TCGA-LUAD.survival.tsv",
                           sep = '\t',
                           header = T,
                           quote = '',
                           fill = T, 
                           comment.char = "!",
                           stringsAsFactors = FALSE)
  LUADsurvival[1]<-gsub('-','.',LUADsurvival[,1])
  names(LUADsurvival)[1]<-"id"
  LUADsurvival<-LUADsurvival[match(colnames(LUADdata),LUADsurvival$id),]  
  LUADsurvival<-LUADsurvival%>%filter(is.na(LUADsurvival$id) == FALSE)
for (i in 1:length(cancername)) {   
    print(i)
    path<-paste("/",cancername[i],"data.csv",sep = "")
    assign(cancernamedata[i],read.csv(path[1]))
  }
  datalist <- setNames(lapply(ls(pattern="[A-Z]{1,5}data"), function(x) if (class(get(x)) == "data.frame")  get(x)),ls(pattern="[A-Z]{1,5}data"))
  datalist <- list(CHOLdata,DLBCdata) 
  datalist<-lapply(datalist,function(df){
    rownames(df) = df$X
    df$X = NULL
    return(df)
  })
  group_list=c(rep('before',18),rep('after',18))
  group_list=factor(group_list)
  datalistnormalize<-lapply(datalist,function(df){
    df<-normalizeBetweenArrays(df)
  })
  for (df in 1:length(cancernamedata)){
    print(df)
    datalist[[df]]<-normalizeBetweenArrays(datalist[[df]])
  }
  datalist_new<-lapply(datalist,function(df){
    df<-as.data.frame(df)
  })
  datalist_new<-lapply(datalist_new,function(df){ 
    df<-as.matrix(df)
  })
  gsvalist<-lapply(datalist_new,function(df){
    gsva(df,MSI_signature, 
         min.sz = 1,
         max.sz = Inf,
         mx.diff=TRUE, 
         verbose=FALSE, 
         parallel.sz=0,
         method = "gsva")
  })
  gsvalist<-lapply(gsvalist,function(df){ 
    df<-as.data.frame(df)
  })
  rm(datalist)
