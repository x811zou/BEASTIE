get_geneID <- function(chr,pos){
  tmp<-test_input_full%>%dplyr::filter(chr==chr)%>%dplyr::filter(start<pos)%>%dplyr::filter(end>pos)
  if(nrow(tmp)!=0){
    ref_geneID<-tmp[1,4]
    gene_symbol<-tmp[1,3]
    return(ref_geneID)  
}
  else{
    return("NA")
  }
}

Find_error <- function(example_sub){
  error=0
  for( i in 1:(dim(example_sub)[1]-1)){
  #print(example_sub[i,])
  #print(example_sub[i+1,])
  if (example_sub[i,1][[1]]==example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]==example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]==example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]==example_sub[i+1,4][[1]]){
      error=error
    }else{
      error=error+1     
      #print(paste0("error : ",error))
      #print(example_sub[i,])
      #print(example_sub[i+1,])
    }
  }else if (example_sub[i,1][[1]]!=example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]!=example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]!=example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]!=example_sub[i+1,4][[1]]){
      error=error
    }else{
      error=error+1    
      #print(paste0("error : ",error))
      #print(example_sub[i,])
      #print(example_sub[i+1,])
    }
  }
  }
  return(error)
}

Report_total_error<-function(data){
  data$PS<-as.character(data$PS)
  data<-data%>%arrange(pos,PS)
  phased_set<-data%>%group_by(PS)%>%summarise(Freq=n())
  phased_set<-phased_set%>%mutate(pair=Freq-1)
  print(paste0("Total pair is:",sum(phased_set$pair))) # 422 pairs
  error_sum=0
  j=1
  for (PS1 in unique(data$PS)){
    #PS1=15965961
    example<-data%>%filter(PS==PS1)
    example_df<-example%>%select(t_paternal,t_maternal,e_paternal,e_maternal)
    error=Find_error(example_df)
    error_sum=error_sum+error
    #print(j)
    j=j+1      
  }  
  print(paste0("Total error number is: ",error_sum))
  print(paste0("Phasing error rate is: ",round(error_sum/sum(phased_set$pair)*100,2),"%"))
}

Report_pair_error<-function(data2){
  data2$PS<-as.character(data2$PS)
  data2<-data2%>%arrange(pos,PS)
  data2$PS<-as.integer(data2$PS)
  phased_set<-data2%>%group_by(PS)%>%summarise(Freq=n())
  phased_set<-phased_set%>%mutate(pair=Freq-1)%>%filter(Freq>=2)
  data3<-data2%>%filter(PS%in%phased_set$PS)
  print(paste0("Total pair is:",sum(phased_set$pair))) # 422 pairs
  error_sum=0
  error_list<-list()
  j=1
  for (PS1 in unique(data3$PS)){
    #print(paste0("phase set:",PS1))
    #PS1=753405
    example<-data3%>%filter(PS==PS1)
    #print(PS1)
    example_df<-example%>%select(t_paternal,t_maternal,e_paternal,e_maternal)
    error=Find_error(example_df)
    errorlist=Find_errorlist(example_df)
    error_list <- c(error_list,errorlist)
    error_sum=error_sum+error
    #print(j)
    j=j+1
  }  
  print(paste0("Total error number is: ",error_sum)) 
  print(paste0("Phasing error rate is: ",round(error_sum/sum(phased_set$pair)*100,2),"%"))
  data3$error<-error_list
  data3<-data3%>%group_by(PS)%>%mutate(distance=pos-lag(pos))
  data3$error<-as.integer(data3$error)
  return(data3)
}

Find_distancelist <- function(data){
  distance<-list()
  distance[[1]]<-0  
  
  error[[i+1]]<-0
}

Find_errorlist <- function(example_sub){
  error<-list()
  error[[1]]<-0
  for( i in 1:(dim(example_sub)[1]-1)){
    if (example_sub[i,1][[1]]==example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]==example_sub[i+1,2][[1]]){
      
      if (example_sub[i,3][[1]]==example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]==example_sub[i+1,4][[1]]){
        #print("combination 1: inphase & inphase")
        #print("    1|0 or 0|1")
        #print("    1|0 or 0|1")
        error[[i+1]]<-0
      }else{
        #print("error")
        error[[i+1]]<-1     
      }
    }else if (example_sub[i,1][[1]]!=example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]!=example_sub[i+1,2][[1]]){
      if (example_sub[i,3][[1]]!=example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]!=example_sub[i+1,4][[1]]){
        #print("combination 2:antiphase & antiphase")
        #print("    1|0 or 0|1")
        #print("    0|1 or 1|0")
        error[[i+1]]<-0
      }else{
        #print("error")
        error[[i+1]]<-1      
      }
    }
  }
  if (length(error)==dim(example_sub)[1]){
    #print("> check")
  }else{
    #print("> length incorrect")
  }
  return(error)
}

extract_PS <- function(string){
  Split <- str_split(string, ":")
  PS=Split[[1]][length(Split[[1]])]
  return(PS)
}

Read_whatshap <- function(directory,sample,N){
  directory=directory
  filename=paste(sample, "_chr",N,".vcf",sep = "")
  whatshap=read.delim(file.path(directory,filename),header=FALSE,sep="\t")
  whatshap$phased <- lapply(strsplit(as.character(whatshap$V10), ":"), "[", 1)
  whatshap$t_paternal <- lapply(strsplit(as.character(whatshap$phased), "\\|"), "[", 1)
  whatshap$t_maternal <- lapply(strsplit(as.character(whatshap$phased), "\\|"), "[", 2)
  whatshap <-whatshap%>%rowwise()%>%mutate(PS=extract_PS(V10))
  #whatshap$PS <- lapply(strsplit(as.character(whatshap$V10), ":"), "[", 2)
  whatshap<-whatshap%>%select(-V10)
  colnames(whatshap) <- c("chr","pos","rsid","ref","alt","qual","filter","infor","format","t_phased","t_paternal","t_maternal","PS")
  whatshap$chr<-as.character(whatshap$chr)
  whatshap$rsid<-as.character(whatshap$rsid)
  whatshap$ref<-as.character(whatshap$ref)
  whatshap$alt<-as.character(whatshap$alt)
  whatshap$pos<-as.integer(whatshap$pos)
  whatshap$paternal<-as.integer(whatshap$t_paternal)
  whatshap$maternal<-as.integer(whatshap$t_maternal)
  whatshap<-whatshap%>%select(chr,pos,rsid,ref,alt,t_paternal,t_maternal,PS)  
  whatshap<-unique(whatshap)
  return(whatshap)
}

Combine_whatshap <- function(sample,N1,N2,directory){
  #directory=paste0("/Users/scarlett/Desktop/example/",sample,"/r_input")
  directory=directory
  df0=Read_whatshap(directory=directory,sample=sample,N1)
  for(i in N1+1:N2-1){
    print(i)
    df=Read_whatshap(directory=directory,sample=sample,i)
    df0 <- rbind(df0,df)
  }
  return(df0)
}

Read_shapeit2 <- function(directory,sample,N){
  directory=paste(directory, "/chr",N,"/phase_withoutseq",sep = "")
  filename=paste(sample, "_chr",N,".phased.with.ref.haps",sep = "")
  phase_withoutseq=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  colnames(phase_withoutseq) <- c("chr","rsid","pos","ref","alt","e_paternal","e_maternal")
  phase_withoutseq<-phase_withoutseq%>%select(chr,pos,rsid,ref,alt,e_paternal,e_maternal)
  phase_withoutseq$chr<-as.character(phase_withoutseq$chr)
  phase_withoutseq$rsid<-as.character(phase_withoutseq$rsid)
  phase_withoutseq$ref<-as.character(phase_withoutseq$ref)
  phase_withoutseq$alt<-as.character(phase_withoutseq$alt)
  phase_withoutseq$pos<-as.integer(phase_withoutseq$pos)
  phase_withoutseq$e_maternal<-as.integer(phase_withoutseq$e_maternal)
  phase_withoutseq$e_paternal<-as.integer(phase_withoutseq$e_paternal)
  phase_withoutseq<-unique(phase_withoutseq)
  return(phase_withoutseq)
}

Combine_shapeit2 <- function(sample,N1,N2,directory){
  directory=directory
  df0=Read_shapeit2(directory=directory,sample=sample,N1)
  for(i in N1+1:N2-1){
    print(i)
    df=Read_shapeit2(directory=directory,sample=sample,i)
    df0 <- rbind(df0,df)
    print(dim(df0))
  }
  return(df0)
}

Report_chr <- function(N){
  sample="NA12878"
  directory="/data/allenlab/scarlett/data/shapeit2/GIAB/NA12878"
  N1=N
  chr1_S=Read_shapeit2(directory=directory,sample=sample,N)
  directory="/Users/scarlett/Desktop/example/HG00096/r_input"
  chr1_W=Read_whatshap(directory=directory,sample=sample,N)
  data<-inner_join(chr1_S,chr1_W,by=c("chr","pos","rsid","ref","alt"))
  data$PS<-as.character(data$PS)
  data2<-unique(data)
  data_new<-Report_pair_error(data2) 
  return(data_new)
}

Report_chr_error <- function(data,source){
  error<-rep(NA,22)
  error_per<-rep(NA,22)
  if(source=="whatshap"){
    print("1. whatshap")
    for(i in seq(22)){
      chromosome=paste0("chr",i)
      data_input<-data%>%filter(chr==chromosome)
      cat(i,chromosome)
      results = Report_pair_error_whatshap(data_input)
      error[i]=results[[1]]
      error_per[i]=results[[2]]
    }
  }else{
    print("2. shapeit")
    for(i in seq(22)){
      chromosome=paste0("chr",i)
      data_input<-data%>%filter(chr==chromosome)
      results = Report_pair_error_shapeit2(data_input)
      cat(i,chromosome)
      error[i]=results[[1]]
      error_per[i]=results[[2]]
    }
  }
  return(list(error,error_per))
}

Read_GIAB <- function(path){
  #path="/Users/scarlett/Desktop/HARDAC/scarlett"
  directory="/data/shapeit2/GIAB"
  directory<-paste0(path,directory,sep="")
  filename="NA12878.SNPs.phased.hets.r_input.vcf"
  df=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  df2<-df
  df2$phased <- lapply(strsplit(as.character(df2$V6), ":"), "[", 1)
  df2<-df2%>%filter(phased!="1/1")
  df2$t_paternal <- lapply(strsplit(as.character(df2$phased), "\\|"), "[", 1)
  df2$t_maternal <- lapply(strsplit(as.character(df2$phased), "\\|"), "[", 2)
  df2<-df2%>%select(-V6)
  colnames(df2) <- c("chr","pos","rsid","ref","alt","phased","paternal","maternal")
  df2$chr<-as.character(df2$chr)
  df2$rsid<-as.character(df2$rsid)
  df2$ref<-as.character(df2$ref)
  df2$alt<-as.character(df2$alt)
  df2$pos<-as.integer(df2$pos)
  df2$t_paternal<-as.integer(df2$t_paternal)
  df2$t_maternal<-as.integer(df2$t_maternal)
  df2<-unique(df2)
  return(df2)
}

Find_errorlist <- function(example_sub){
  error<-list()
  error[[1]]<-0
  for( i in 1:(dim(example_sub)[1]-1)){
    if (example_sub[i,1][[1]]==example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]==example_sub[i+1,2][[1]]){
      
      if (example_sub[i,3][[1]]==example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]==example_sub[i+1,4][[1]]){
        #print("combination 1: inphase & inphase")
        #print("    1|0 or 0|1")
        #print("    1|0 or 0|1")
        error[[i+1]]<-0
      }else{
        #print("error")
        error[[i+1]]<-1     
      }
    }else if (example_sub[i,1][[1]]!=example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]!=example_sub[i+1,2][[1]]){
      if (example_sub[i,3][[1]]!=example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]!=example_sub[i+1,4][[1]]){
        #print("combination 2:antiphase & antiphase")
        #print("    1|0 or 0|1")
        #print("    0|1 or 1|0")
        error[[i+1]]<-0
      }else{
        #print("error")
        error[[i+1]]<-1      
      }
    }
  }
  if (length(error)==dim(example_sub)[1]){
    #print("> check")
  }else{
    #print("> length incorrect")
  }
  return(error)
}

Find_errorlist_shapeit2 <- function(example_sub){
  error<-list()
  i=1
  if (example_sub[i,1][[1]]==example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]==example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]==example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]==example_sub[i+1,4][[1]]){
      error[[i+1]]<-0
    }else{
      error[[i+1]]<-1     
    }
  }else if (example_sub[i,1][[1]]!=example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]!=example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]!=example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]!=example_sub[i+1,4][[1]]){
      error[[i+1]]<-0
    }else{
      error[[i+1]]<-1      
    }
  }
  return(error)
}

Report_pair_error_shapeit2 <- function(data){
  error_sum=0
  error_list<-list()
  error_list[[1]]<-0
  data$pos<-as.integer(data$pos)
  for (i in 2:dim(data)[1]){
    #print(i)
    example <- data[c(i-1,i),]
    example_df<-example%>%select(t_paternal,t_maternal,e_paternal,e_maternal)
    error=Find_error(example_df)
    errorlist=Find_errorlist_shapeit2(example_df)
    error_list[[i]]<-error
    error_sum=error_sum+error
  }  
  print(paste0("Total error number is: ",error_sum," out of ",dim(data)[1])) 
  print(paste0("Phasing error rate is: ",round(error_sum/dim(data)[1]*100,2),"%"))
  data$error<-error_list
  data<-data%>%group_by(chr)%>%mutate(distance=pos-dplyr::lag(pos))
  data$error<-as.integer(data$error)
  #if (return_item=="data"){
  #  return(data)
  #}else{
  return(list(error_sum,round(error_sum/dim(data)[1]*100,2)))
  #}
}

Find_error <- function(example_sub){
  error=0
  for( i in 1:(dim(example_sub)[1]-1)){
  #print(example_sub[i,])
  #print(example_sub[i+1,])
  if (example_sub[i,1][[1]]==example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]==example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]==example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]==example_sub[i+1,4][[1]]){
      error=error
    }else{
      error=error+1     
      #print(paste0("error : ",error))
      #print(example_sub[i,])
      #print(example_sub[i+1,])
    }
  }else if (example_sub[i,1][[1]]!=example_sub[i+1,1][[1]]&&example_sub[i,2][[1]]!=example_sub[i+1,2][[1]]){
    if (example_sub[i,3][[1]]!=example_sub[i+1,3][[1]]&&example_sub[i,4][[1]]!=example_sub[i+1,4][[1]]){
      error=error
    }else{
      error=error+1    
      #print(paste0("error : ",error))
      #print(example_sub[i,])
      #print(example_sub[i+1,])
    }
  }
  }
  return(error)
}

Report_pair_error_whatshap<-function(data2,return_item){
  data2$PS<-as.character(data2$PS)
  data2<-data2%>%arrange(pos,PS)
  data2$PS<-as.integer(data2$PS)
  data2$paternal<-as.integer(data2$paternal)
  data2$maternal<-as.integer(data2$maternal)
  data2$t_paternal<-as.integer(data2$t_paternal)
  data2$t_maternal<-as.integer(data2$t_maternal)
  phased_set<-data2%>%group_by(PS)%>%summarise(Freq=n())
  phased_set<-phased_set%>%mutate(pair=Freq-1)%>%filter(Freq>=2)
  data3<-data2%>%filter(PS%in%phased_set$PS)
  print(paste0("Total pair is:",sum(phased_set$pair))) # 422 pairs
  error_sum=0
  error_list<-list()
  j=1
  for (PS1 in unique(data3$PS)){
    #print(paste0("phase set:",PS1))
    #PS1=941284
    example<-data3%>%filter(PS==PS1)
    #print(PS1)
    example_df<-example%>%select(paternal,maternal,t_paternal,t_maternal)
    error=Find_error(example_df)
    errorlist=Find_errorlist(example_df)
    error_list <- c(error_list,errorlist)
    error_sum=error_sum+error
    #print(j)
    j=j+1
  }  
  print(paste0("Total error number is: ",error_sum)) 
  print(paste0("Phasing error rate is: ",round(error_sum/sum(phased_set$pair)*100,2),"%"))
  data3$error<-error_list
  data3<-data3%>%group_by(PS)%>%mutate(distance=pos-lag(pos))
  data3$error<-as.integer(data3$error)
  if (return_item=="data"){
    return(data3)
  }else{
  return(list(error_sum,round(error_sum/sum(phased_set$pair)*100,2)))
  }
}

Read_AF <- function(directory,sample,N){
  directory=paste(directory, "/chr",N,"/vcf",sep = "")
  filename=paste("1KGP_AF_",sample, "_biSNPs_chr",N,"_phased.vcf",sep = "")
  AF=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  colnames(AF) <- c("chr","rsid","pos","ref","alt","AF_info")
  AF$AF_extract <- lapply(strsplit(as.character(AF$AF_info), ";"), "[", 2)
  AF$AF <- lapply(strsplit(as.character(AF$AF_extract), "AF="), "[", 2)
  AF$AF<-as.numeric(AF$AF)
  AF$EUR_AF_extract <- lapply(strsplit(as.character(AF$AF_info), ";"), "[", 9)
  AF$EUR_AF <- lapply(strsplit(as.character(AF$EUR_AF_extract), "EUR_AF="), "[", 2)
  AF$EUR_AF<-as.numeric(AF$EUR_AF)  
  AF$chr<-as.character(AF$chr)
  AF$rsid<-as.character(AF$rsid)
  AF$pos<-as.integer(AF$pos)
  AF$ref<-as.character(AF$ref)
  AF$alt<-as.character(AF$alt)
  AF<-AF%>%select(-AF_extract,-EUR_AF_extract,-AF_info)
  return(AF)
}

Combine_AF <- function(sample,N1,N2){
  directory="/Users/scarlett/Desktop/HARDAC/scarlett/data/shapeit2/GIAB/NA12878"
  df0=Read_AF(directory=directory,sample=sample,N1)
  for(i in N1+1:N2-1){
    print(i)
    df=Read_AF(directory=directory,sample=sample,i)
    df0 <- rbind(df0,df)
  }
  return(df0)
}

save_file <- function(file,filename,directory){
  directory=directory
  filename=paste0(directory,filename,".tsv")
  file=as.matrix(file)
  write.table(file,filename,sep = "\t", row.names = FALSE, col.names = TRUE)
}



Read_triovcf <- function(path,directory,filename){
  directory<-paste0(path,directory,sep="")
  # clean
  # only kept 1|0 & 0|1
  df=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  #935566 0|1
  #947443 1|0
  df2<-df%>%select(-V6,-V7,-V8,-V9)
  colnames(df2)<-c("chr","pos","rsid","ref","alt","NA19238","NA19239","NA19240")
  df2$paternal <- lapply(strsplit(as.character(df2$NA19240), "\\|"), "[", 1)
  df2$maternal <- lapply(strsplit(as.character(df2$NA19240), "\\|"), "[", 2)
  
  df2$chr<-as.character(df2$chr)
  df2$rsid<-as.character(df2$rsid)
  df2$ref<-as.character(df2$ref)
  df2$alt<-as.character(df2$alt)
  df2$pos<-as.integer(df2$pos)
  df2$paternal<-as.integer(df2$paternal)
  df2$maternal<-as.integer(df2$maternal)
  df2<-unique(df2) # 1883009
  return(df2)
}

Combine_shapeit2 <- function(sample,N1,N2,directory){
  directory=directory
  df0=Read_shapeit2(directory=directory,sample=sample,N1)
  for(i in N1+1:N2-1){
    #print(i)
    df=Read_shapeit2(directory=directory,sample=sample,i)
    df0 <- rbind(df0,df)
    #print(dim(df0))
  }
  print(paste0(directory," reference panel has ",dim(df0)[1]," phased haplotype"))
  return(df0)
}


Read_triovcf <- function(path,directory,filename){
  directory<-paste0(path,directory,sep="")
  # clean
  # only kept 1|0 & 0|1
  df=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  #935566 0|1
  #947443 1|0
  df2<-df%>%select(-V6,-V7,-V8,-V9)
  colnames(df2)<-c("chr","pos","rsid","ref","alt","NA19238","NA19239","NA19240")
  df2$paternal <- lapply(strsplit(as.character(df2$NA19240), "\\|"), "[", 1)
  df2$maternal <- lapply(strsplit(as.character(df2$NA19240), "\\|"), "[", 2)
  
  df2$chr<-as.character(df2$chr)
  df2$rsid<-as.character(df2$rsid)
  df2$ref<-as.character(df2$ref)
  df2$alt<-as.character(df2$alt)
  df2$pos<-as.integer(df2$pos)
  df2$paternal<-as.integer(df2$paternal)
  df2$maternal<-as.integer(df2$maternal)
  df2<-unique(df2) # 1883009
  return(df2)
}

Combine_shapeit2 <- function(sample,N1,N2,directory){
  directory=directory
  df0=Read_shapeit2(directory=directory,sample=sample,N1)
  for(i in N1+1:N2-1){
    #print(i)
    df=Read_shapeit2(directory=directory,sample=sample,i)
    df0 <- rbind(df0,df)
    #print(dim(df0))
  }
  print(paste0(directory," reference panel has ",dim(df0)[1]," phased haplotype"))
  return(df0)
}

Find_exonstart <- function(pos,chrNum){
  exons$pos<-as.integer(exons$pos)
  exons$exon_start<-as.integer(exons$exon_start)
  exons$exon_end<-as.integer(exons$exon_end)
  pos<-as.integer(pos)
  data<-exons%>%dplyr::filter(chr==chrNum)%>%dplyr::select(exon_start<pos)%>%unique()
  start<-data%>%select(exon_start)%>%unique()
  rsid_pos <- as.integer(pos)
  #rsid_pos<-19930639 #,19933261,19935449,19950062
  #print(rsid_pos)
  start<-start[which(start<=rsid_pos)]
  if (length(start)>0){
    idx = which(abs(start-rsid_pos)==min(abs(rsid_pos-start)))[1]
    identified_start<-start[idx]
    data<-data%>%filter(exon_start==identified_start)
    data<-data[1,]
    if(data$exon_end>rsid_pos){
      data2<-data%>%filter(exon_start==identified_start)
      Start<-data2$exon_start 
      End<-data2$exon_end
      Gene<-data2$geneID
      return(list(Start,End,Gene))
    }else{
      return(list(NA,NA,NA))
      }
    }else{
    return(list(NA,NA,NA))
  }
}