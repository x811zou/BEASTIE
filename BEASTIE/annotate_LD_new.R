suppressMessages(library(LDlinkR))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library("dplyr"))
library("LDlinkR")

args = commandArgs(trailingOnly=TRUE)
prefix=args[1]
ancestry=args[2]
infile=args[3]
outpath=args[4]
mytoken=args[5]
chr_start=args[6]
chr_end=args[7]
meta=args[8]
beastie_wd=args[9]

source(file.path(beastie_wd, "Get_LD.R"))

if (!file.exists(meta)){
#================================================ specify input
    df = read.table(infile, header=T, sep = '\t')
    batch_size = 1000
    in_data=df

    #================================================ start running
    if (ancestry=="EUR"){
        ancestry="CEU"
    }
    if (ancestry=="AFR"){
        ancestry="YRI"
    }
    print(paste0("start annotating LD/AF information for ",prefix," in ",ancestry,sep=""))
    processed_exon<-process_data(in_data)
    data<-processed_exon
    data$lag_pos=lag(data$pos)
    data$lag_rs=lag(data$rsid)
    data<-data%>%dplyr::select(geneID,chr,chrN,pos,lag_pos,rsid,lag_rs)%>%
        dplyr::arrange(chr,pos)%>%
        dplyr::mutate(snp=paste(chr,pos,sep=":"))%>%
        dplyr::mutate(snp_pair=paste(rsid,lag_rs,sep="-"))
    data$lag_snp=lag(data$snp)
    data<-data%>%group_by(chr)%>%dplyr::mutate(lag_pos=ifelse(pos==dplyr::first(pos),NA,lag_pos))%>%ungroup()
    data<-data%>%group_by(geneID)%>%dplyr::mutate(lag_pos=ifelse(pos==dplyr::first(pos),NA,lag_pos))%>%ungroup()

    #================================================ loop through each chromosome
    data_LD=loop_chr(outpath,data,ancestry,batch_size,mytoken,chr_start,chr_end)
    #================================================ re-look at NAs
    data_w_LD<-left_join(processed_exon,data_LD)

    data_final_reginput<-data_w_LD%>%
        mutate(distance=pos-lag_pos)%>%
        mutate(log10_distance=log10(distance))%>%
        mutate(MAF=ifelse(AF>0.5,1-AF,AF))                 # the minimum AF for the current variant

    data_final_reginput$lag_MAF=lag(data_final_reginput$MAF)

    data_final_reginput<-data_final_reginput%>%            # MAF for the previous variant
        mutate(min_MAF=ifelse(MAF>lag_MAF,lag_MAF,MAF))%>% # the mininum between the pair
        mutate(diff_MAF=abs(MAF-lag_MAF))                  # the difference between the pair

    data_final_reginput<-data_final_reginput%>%
        group_by(geneID)%>%
        arrange(chr,pos)%>%
        mutate(diff_pos=pos-lag_pos)%>%
        ungroup()

    data_final_reginput<-data_final_reginput%>%select(geneID,chr,pos,lag_pos,rsid,log10_distance,d,r2,MAF,lag_MAF,min_MAF,diff_MAF)

    write.table(data_final_reginput,meta,sep = "\t", row.names = FALSE, col.names = TRUE)
}else{
    print("file exists!")
}
