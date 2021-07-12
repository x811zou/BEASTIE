
save_file <- function(file,filename,directory){
  directory=directory
  filename=paste0(directory,filename,".tsv")
  file=as.matrix(file)
  write.table(file,filename,sep = "\t", row.names = FALSE, col.names = TRUE)
}

process_data<-function(data){
  data$chr<-factor(data$chr,levels = c("chr1","chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12", "chr13", "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))   
  data$pos<-as.integer(data$pos)
  data<-data%>%group_by(chr,pos)%>%slice(1)
  return(data)
}


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

Combine_LD <- function(N1,N2,directory){
  filename<-paste0("out","1",".results",sep="")
  df0<-read.delim(file.path(directory,filename),header=TRUE,sep="\t")  
  df0<-df0%>%select(-i)
  print(1)
  df0<-df0%>%unique()
  for(i in N1+2:N2-1){
    print(i)
    filename<-paste0("out",i,".results",sep="")
    df=read.delim(file.path(directory,filename),header=TRUE,sep="\t")
    df<-df%>%select(-i)
    df<-df%>%unique()
    df0 <- rbind(df0,df)
  }
  return(df0)
}

Get_LD <- function(snp,lag_snp,mytoken,pop,r2d){
    if(is.na(snp)|is.na(lag_snp)){
      return(NA)
    }
    else{
      mat <- LDlinkR::LDmatrix(snps = c(lag_snp,snp), pop = pop, r2d = r2d, token = mytoken, file = FALSE)
      string<-mat[1,3]
      if(is.null(string)){
          for(i in 1:10){
              mat <- LDlinkR::LDmatrix(snps = c(lag_snp,snp), pop = pop, r2d = r2d, token = mytoken, file = FALSE)
              string <- mat[1,3]
              if(!is.null(string)){
                  return(as.numeric(string))
              }
          }
          return(NA)
      }else{
          return(as.numeric(string))
      }
    }
}

loop_chr <-function(outpath,data,ancestry,batch_size,mytoken,chr_start,chr_end){
  final_df_filename=paste0(outpath,prefix,"_chr_",chr_start,"_",chr_end,"_LD.tsv",sep="")
  if (!file.exists(final_df_filename)) {
    counter=0
    first_counter=0
    for(j in c(chr_start:chr_end)){
      print(paste0(">>>>>>>>>> start chr",j,"!"))
      data_chr <- data%>%filter(chrN==j)
      if(dim(data_chr)[1]>0){
        counter=counter+1
        df_LD_chr <- generate_LD(j,data_chr,ancestry,batch_size,mytoken,outpath)
        if(counter==1){
          first_counter=j
          last_counter=first_counter
          df_LD<-df_LD_chr
        }else{
          df_LD<-rbind(df_LD,df_LD_chr)
          ff=last_counter+1
          last_counter=ff
        }
        print(paste0(">>>>>>>>>> finish chr",j,"!"))
      }else{
        print(paste0(">>>>>>>>>> no content on chr",j,"!")) 
        break
      }   
      print(paste0("chr",j," done!",sep=""))
    }
    tmp_df_filename=paste0(outpath,prefix,"_chr_",first_counter,"_",last_counter,"_LD.tsv",sep="")
    write.table(df_LD, tmp_df_filename, row.names = FALSE,quote=FALSE, sep='\t')
    df_LD=read.table(file = tmp_df_filename, sep = '\t', header = TRUE)
    if((dim(df_LD)[1]>1)&&(first_counter==chr_start)&&(last_counter==chr_end)){
      print(paste0("file save to ",final_df_filename,sep=""))
      for(file in list.files(path = outpath)){
        if(grepl("temp_chr", file)){
          file.remove(paste0(outpath,file,sep=""))
      }else{
        print(paste0("saved ",tmp_df_filename,"not matching to ",final_df_filename," PLEASE CHECK!",sep=""))
      }
    }
    }
  }else{
    print(paste0("(file exists!) ",final_df_filename,sep=""))
    df_LD=read.table(file = final_df_filename, sep = '\t', header = TRUE)
  }
  return(df_LD)
}


calculate_LD_V2 <- function(df,pop,r2d,mytoken){
    df_out = LDlinkR::LDmatrix(snps = df$snp, pop = pop, r2d = r2d, token = mytoken, file = FALSE)
    if(dim(df_out)[1]==3){
      print(df_out)
    }
    df_snps = reshape2::melt(df_out)
    df_snps$snp_pair = paste(df_snps$RS_number,df_snps$variable, sep="-")
    mat_out = merge(df, df_snps, by="snp_pair")
    return(mat_out)
  }

generate_LD <- function(chr,data,ancestry,batch_size,mytoken,outpath){
  final_df_filename=paste0(outpath,"/temp_chr",chr,"_LD.tsv",sep="")
  if (file.exists(final_df_filename)) {
    final_df=read.table(file = final_df_filename, sep = '\t', header = TRUE)
  }else{
    print(paste0(prefix," chr",chr," input data size: ",dim(data)[1],sep=""))
    rownames(data) <- NULL
    start_time <- Sys.time()
    df_list_d = list()
    count = 1
    print(" >>>> start processing!")
    while(length(df_list_d)!=ceiling(nrow(data)/batch_size)){
      x <- for (i in seq(1, nrow(data), batch_size)){
        start = i
        end = min(i+batch_size-1, nrow(data))
        temp_d_i_name=paste0(outpath,"/temp_chr",chr,"_d_",count,".tsv",sep="")
        if(!file.exists(temp_d_i_name)){
          if(dim(data[start:end,])[1]<2){
            temp_d_i = data[start:end,]
            temp_d_i <- temp_d_i%>%rowwise()%>%dplyr::mutate(d=Get_LD(snp,lag_snp,mytoken,ancestry,"d"))
          }else{
            temp_d_i = calculate_LD_V2(data[start:end,],ancestry,"d",mytoken)
          }
          write.table(temp_d_i, temp_d_i_name, row.names = FALSE,quote=FALSE, sep='\t')
          print(paste0(prefix," chr",chr," : d - batch ",count," LD >>>> just processed: ",round(end/nrow(data)*100,2)," %",sep=""))
        }else{
          temp_d_i=read.table(file = temp_d_i_name, sep = '\t', header = TRUE)
          print(paste0(prefix," chr",chr," : d - batch ",count," LD >>>> already processed: ",round(end/nrow(data)*100,2)," %",sep=""))
        }
        df_list_d[[count]] = temp_d_i
        count = count+1
        }
    }
    final_df_d = bind_rows(df_list_d, .id = "column_label")
    final_df_d<-final_df_d%>%select(geneID,chr,chrN,pos,lag_pos,rsid,snp,lag_snp,lag_rs,snp_pair,value)
    colnames(final_df_d)[ncol(final_df_d)]<-"d"

    df_list_r2 = list()
    count = 1
    while(length(df_list_r2)!=ceiling(nrow(data)/batch_size)){
      x <- for (i in seq(1, nrow(data), batch_size)){
        start = i
        end = min(i+batch_size-1, nrow(data))
        temp_r2_i_name=paste0(outpath,"/temp_chr",chr,"_r2_",count,".tsv",sep="")
        if(!file.exists(temp_r2_i_name)){
          if(dim(data[start:end,])[1]<2){
            temp_r2_i = data[start:end,]
            temp_r2_i <- temp_r2_i%>%rowwise()%>%dplyr::mutate(d=Get_LD(snp,lag_snp,mytoken,ancestry,"r2"))
          }else{
            temp_r2_i = calculate_LD_V2(data[start:end,],ancestry,"r2",mytoken)
          }
          write.table(temp_r2_i, temp_r2_i_name, row.names = FALSE,quote=FALSE, sep='\t')
          print(paste0(prefix," chr",chr," : r2 - batch ",count," LD >>>> just processed: ",round(end/nrow(data)*100,2)," %",sep=""))
        }else{
          temp_r2_i=read.table(file = temp_r2_i_name, sep = '\t', header = TRUE)
          print(paste0(prefix," chr",chr," : r2 - batch ",count," LD >>>> already processed: ",round(end/nrow(data)*100,2)," %",sep=""))
        }
        df_list_r2[[count]] = temp_r2_i
        count = count+1
        }
    }
    final_df_r2 = bind_rows(df_list_r2, .id = "column_label")
    final_df_r2<-final_df_r2%>%select(geneID,chr,chrN,pos,lag_pos,rsid,snp,lag_snp,lag_rs,snp_pair,value)
    colnames(final_df_r2)[ncol(final_df_r2)]<-"r2"
  ###########    
    final_df<-inner_join(final_df_d,final_df_r2)
    end_time <- Sys.time()
    runtime = round(end_time - start_time,2)
    print(paste0("batch size: ",batch_size," uses ",runtime," min for ",dim(data)[1]," SNPs!",sep=""))
    print(paste0("output complete data size: ",dim(final_df)[1],sep=""))
    write.table(final_df, final_df_filename, row.names = FALSE,quote=FALSE, sep='\t')
    print(paste0(final_df_filename," saved!",sep=""))
    tmp_r2=paste0("temp_chr",chr,"_r2_",sep="")
    tmp_d=paste0("temp_chr",chr,"_d_",sep="")
    for(file in list.files(path = outpath)){
      if(grepl(tmp_r2, file)){
        file.remove(paste0(outpath,file,sep=""))
      }
      if(grepl(tmp_d, file)){
        file.remove(paste0(outpath,file,sep=""))
      }
    }
}
  return(final_df)
}


