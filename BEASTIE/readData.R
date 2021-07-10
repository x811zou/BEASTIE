Read_GIAB <- function(filename){
  path="/Users/scarlett/Desktop/HARDAC/scarlett"
  directory="/data/shapeit2/GIAB"
  directory<-paste0(path,directory,sep="")
  #filename="NA12878.SNPs.phased.hets.r_input.vcf"
  df=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  df2<-df
  df2$phased <- lapply(strsplit(as.character(df2$V6), ":"), "[", 1)
  df2<-df2%>%filter(phased!="1/1")
  df2$t_paternal <- lapply(strsplit(as.character(df2$phased), "\\|"), "[", 1)
  df2$t_maternal <- lapply(strsplit(as.character(df2$phased), "\\|"), "[", 2)
  df2<-df2%>%select(-V6)
  colnames(df2) <- c("chr","pos","rsid","ref","alt","phased","t_paternal","t_maternal")
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

#
Read_giab_shapeitphased <- function(ref,sample){
  path="/Users/scarlett/Desktop/HARDAC/scarlett/data/shapeit2/GIAB/"
  directory=paste0(path,ref,"/",sample)
  statistical<-Combine_shapeit2(sample,1,22,directory)               #15,084 (1-22)
  statistical$chr<-as.factor(statistical$chr)
  reorderFactors(statistical,"chr", desired_level_order = c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12", "chr13", "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
  statistical<-statistical%>%arrange(chr,pos)%>%unique()
  return(statistical)
}

Combine_AF <- function(sample,N1,N2,ancestry){
  directory="/Users/scarlett/Desktop/HARDAC/scarlett/data/shapeit2/GIAB/NA12878"
  df0=Read_AF(directory=directory,sample=sample,N1,ancestry)
  for(i in N1+1:N2-1){
    print(i)
    df=Read_AF(directory=directory,sample=sample,i,ancestry)
    df0 <- rbind(df0,df)
  }
  return(df0)
}

Read_AF <- function(directory,sample,N,ancestry){
  directory=paste(directory, "/chr",N,"/vcf",sep = "")
  filename=paste("1KGP_AF_",sample, "_biSNPs_chr",N,"_phased.vcf",sep = "")
  AF=read.delim(file.path(directory,filename),header=FALSE,sep=" ")
  colnames(AF) <- c("chr","rsid","pos","ref","alt","AF_info")
  AF$AF_extract <- lapply(strsplit(as.character(AF$AF_info), ";"), "[", 2)
  AF$AF <- lapply(strsplit(as.character(AF$AF_extract), "AF="), "[", 2)
  AF$AF<-as.numeric(AF$AF)
  AF$ancestry_AF_extract <- lapply(strsplit(as.character(AF$AF_info), ";"), "[", 9)
  Ancestry=paste0(ancestry,"_AF=",sep="")
  AF$ancestry_AF <- lapply(strsplit(as.character(AF$ancestry_AF_extract), Ancestry), "[", 2)
  AF$ancestry_AF<-as.numeric(AF$ancestry_AF)  
  AF$chr<-as.character(AF$chr)
  AF$rsid<-as.character(AF$rsid)
  AF$pos<-as.integer(AF$pos)
  AF$ref<-as.character(AF$ref)
  AF$alt<-as.character(AF$alt)
  AF<-AF%>%select(-AF_extract,-ancestry_AF_extract,-AF_info)
  return(AF)
}

save_file <- function(file,filename,directory){
  directory=directory
  filename=paste0(directory,filename,".tsv")
  file=as.matrix(file)
  write.table(file,filename,sep = "\t", row.names = FALSE, col.names = TRUE)
}

clean_exon <- function(exons){
  exons<-exons%>%select(chr,geneID,exon_start,exon_end)
  exons$chr<-factor(exons$chr,levels = c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12", "chr13", "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")) 
  exons$exon_start<-as.integer(exons$exon_start)
  exons$exon_end<-as.integer(exons$exon_end)
  exons<-exons%>%arrange(chr,exon_start,exon_end)
  exons_uniqe<-exons%>%group_by(chr,geneID,exon_start,exon_end)%>%slice(1)%>%ungroup()
  return(exons_uniqe)
}

process_data<-function(data){
  #data<-data%>%select(chr,pos)
  data$chr<-factor(data$chr,levels = c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12", "chr13", "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))   
  data$pos<-as.integer(data$pos)
  #data_unique<-data%>%arrange(chr,pos)%>%unique()
  data<-data%>%group_by(chr,pos)%>%slice(1)
  return(data)
}

clean_for_LD<-function(data){
  data4<-data
  data4$chr<-as.character(data4$chr)
  data4$pos<-as.integer(data4$pos)
  #data4$LD_rsid<-as.character(data4$LD_rsid)
  data4$chr<-factor(data4$chr,levels = c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12", "chr13", "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
  data4$pos<-as.integer(data4$pos)
  data4<-data4%>%arrange(chr,pos)
  return(data4)
}

