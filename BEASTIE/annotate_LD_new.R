# prefix="HG00096_chr20"
# ancestry="EUR"
###################################################
#Rscript --vanilla annotated_LD_new.R "NA12878" "EUR" hetSNP_AF_file "/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/NA12878/TEMP/"
###################################################
#Rscript annotated_LD_new.R "NA19247" "AFR" "/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/NA19247/TEMP/NA19247_hetSnps_AF.tsv" "/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/NA19247/TEMP/"
###################################################
#Rscript annotated_LD_new.R "HG00096" "EUR" "/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/HG00096/TEMP/HG00096_hetSnps_AF.tsv" "/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/HG00096/TEMP/"
###################################################
suppressMessages(library(LDlinkR))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library("dplyr"))
library("LDlinkR")
source("readData.R")
source("Get_LD.R")

args = commandArgs(trailingOnly=TRUE)
prefix=args[1]#"HG00096_chr20"
ancestry=args[2]#
infile=args[3] #
#infile=paste0("/Users/scarlett/Documents/Allen_lab/BEASTIE_example/HG00096_chr20/output/TEMP/",prefix,"_hetSNP_intersect_unique_AF.tsv",sep="")
outpath=args[4]
#outpath="/Users/scarlett/Documents/Allen_lab/BEASTIE_example/HG00096_chr20/output/TEMP/"
mytoken=args[5]
chr_start=args[6]
chr_end=args[7]

#in_data<-read.delim(infile,header=TRUE,row.names = NULL,"\t")

outData = paste0(prefix,"_logisticReg_input.tsv",sep="")
target_data=paste0(outpath,outData,sep="")

if (!file.exists(target_data)){
    #================================================ specify input
    #df = read.table("/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/example/NA12878/TEMP/NA12878_hetSnps_AF.tsv", header=T, sep = '\t', row.names = 'X')
df = read.table(infile, header=T, sep = '\t', row.names = 'X')
batch_size = 1000
in_data=df#[1:100,]
# if(prefix=="NA12878"){
#     mytoken="1a2b434466cb"
# }
# if(prefix=="NA19247"){
#     mytoken="c313799c13c3"
# }
# if(prefix=="HG00096_chr20"){
#     mytoken="6a89f99c6b1b"
# }
# if(prefix=="HG00096"){
#     mytoken="027b629de521"
# }

#================================================ start running
if (ancestry=="EUR"){
    ancestry="CEU"
}
if (ancestry=="AFR"){
    ancestry="YRI"
}
# if(prefix=="NA12878"){
#     in_data<-in_data%>%select(-rsid_y)
#     colnames(in_data)[5]<-"rsid"
# }
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

#[1] "batch size: 100 uses 1.78512 s for 93 SNPs!"
data_w_LD<-left_join(processed_exon,data_LD)
    #outData = paste0(prefix,"_hetSnps_AF_LD_",ancestry,"_",dim(data_w_LD)[1],sep="")
    #save_file(data_w_LD,outData,outpath)

    # prepare_model_input <- function(){
    #     m
    # }

data_final_reginput<-data_w_LD%>%
    mutate(distance=pos-lag_pos)%>%
    mutate(log10_distance=log10(distance))%>%
    mutate(MAF=ifelse(AF>0.5,1-AF,AF))    # the minimum AF for the current variant

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
    outData = paste0(prefix,"_logisticReg_input",sep="")
    save_file(data_final_reginput,outData,outpath)
    print(paste0("output saved in ",outData,sep=""))

    #================================================ save data
    # data_w_LD_NA <- data_w_LD%>%filter(is.na(d)|is.na(r2))%>%
    #     filter(!is.na(snp))%>%filter(!is.na(lag_snp))
    # print(" >>>> re-looking the NAs!")
    # data_w_LD_NA<-data_w_LD_NA%>%rowwise()%>%
    #     dplyr::mutate(d=Get_LD(snp,lag_snp,mytoken,ancestry,"d"))%>%
    #     dplyr::mutate(r2=Get_LD(snp,lag_snp,mytoken,ancestry,"r2"))
    # print(" >>>> finish looking at the NAs!")

    # data_final<-left_join(data_w_LD,data_w_LD_NA)
    # save_file(data_final,outData,outpath)

    # data_final_reginput<-data_final%>%
    #     mutate(distance=pos-lag_pos)%>%
    #     mutate(log10_distance=log10(distance))%>%
    #     mutate(MAF=ifelse(AF>0.5,1-AF,AF))    # the minimum AF for the current variant

    # data_final_reginput$lag_MAF=lag(data_final_reginput$MAF)

    # data_final_reginput<-data_final_reginput%>%            # MAF for the previous variant
    #     mutate(min_MAF=ifelse(MAF>lag_MAF,lag_MAF,MAF))%>% # the mininum between the pair
    #     mutate(diff_MAF=abs(MAF-lag_MAF))                  # the difference between the pair

    # data_final_reginput<-data_final_reginput%>%
    #     group_by(geneID)%>%
    #     arrange(chr,pos)%>%
    #     mutate(diff_pos=pos-lag_pos)%>%
    #     ungroup()

    # data_final_reginput<-data_final_reginput%>%select(geneID,chr,transcript_pos,pos,lag_pos,rsid,log10_distance,d,r2,MAF,lag_MAF,min_MAF,diff_MAF)
    # outData = paste0(prefix,"_logisticReg_input",sep="")
    # save_file(data_final_reginput,outData,outpath)
    # print(paste0("output saved in ",outData,sep=""))

}else{
    print("file exists!")
}
