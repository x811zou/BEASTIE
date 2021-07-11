#!/usr/bin/env Rscript

#python style read in parameters:
#https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/#

# cmd="Rscript --vanilla predict_lambda_phasingError.R %s %s %s"%(out,prefix,alpha)
# os.system(cmd)
# #input: NA12878_parsed_mpileup_exonicHetSnps_model_input.forLambdaPred.csv
# #output: NA12878_parsed_mpileup_exonicHetSnps_model_input.LambdaPredicted.csv
# #      geneID      |tota_reads| num_hets | predicted lambda
# #"ENSG00000177757.1"     2	        1	   1.3978282039415

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

#### call library
suppressMessages(library("pasilla"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("LDlinkR"))
library(glmnetUtils)
#mytoken="c313799c13c3"
#suppressMessages(library(foreach))
#suppressMessages(library(doParallel))
#suppressMessages(library(parallel))
#### read in self-defined functions
#### read in parameters
alpha=as.numeric(args[1])                          #alpha=0.05
tmp=paste0(args[2],"/",sep="")                     #tmp="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/"
sample=args[3]                                     #sample="HG00096_chr20"
model=args[4]                                      #model="iBEASTIE2"
hetSNP_intersect_unique=args[5]                    #hetSNP_intersect_unique="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/HG00096_chr20_hetSNP_intersect_unique.tsv"
hetSNP_intersect_unique_forlambda_file=args[6]     #hetSNP_intersect_unique_forlambda_file="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/HG00096_chr20_hetSNP_intersect_unique_forLambda.tsv"
hetSNP_intersect_unique_lambdaPredicted_file=args[7]#hetSNP_intersect_unique_lambdaPredicted_file="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/HG00096_chr20_hetSNP_intersect_unique_lambdaPredicted.tsv"
meta=args[8]                                       #meta="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/HG00096_chr20_logisticReg_input.tsv"
meta_error=args[9]                                 #meta_error="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/HG00096_chr20_logisticReg_input_w_error.tsv"

source("Get_phasing_error_rate.R")
source("Get_LD.R")

predict_lambda_realdata <- function(alpha,in_data,out_data){
  #colnames(in_data)<-c("gene_ID","total_reads","num_hets")
  data<-in_data%>%dplyr::mutate(predicted_lambda=(log(alpha/(1-alpha)) -(as.numeric(lambda.fit.simulation$coefficients[1])+as.numeric(lambda.fit.simulation$coefficients[3])*as.integer(totalCount)))/as.numeric(lambda.fit.simulation$coefficients[2]))
  #data<-in_data%>%mutate(predicted_lambda=(log(alpha/(1-alpha)) -(15.587909+-0.006483*total_reads))/-13.248682)
  write.table(data,file = out_data,row.names=FALSE,col.names = FALSE,sep="\t")
  print(paste0("model input with predicted lambda saved to ",out_data,sep=""))
  return(data)
}

############################################## 1. predict lambda #####################################################
if(grepl("iBEASTIE", model, fixed=TRUE)){
  lambda.fit.simulation<- readRDS("LinearReg_iBEASTIE_fitted_model_lambda.rds")
}else{
  lambda.fit.simulation<- readRDS("LinearReg_BEASTIE_fitted_model_lambda.rds")
}

in_data<-read.delim(file.path(hetSNP_intersect_unique_forlambda_file),header=TRUE,sep="\t") 
out_data<-hetSNP_intersect_unique_lambdaPredicted_file
if (!file.exists(out_data)) {
  predicted_df=predict_lambda_realdata(alpha,in_data,out_data)
}else{
  print("lambda prediction file exists!")
}

################################################ 2.1 predict phasing error ###################################
if (!file.exists(meta_error)) {
  sample_info<-read.delim(file.path(meta),header=TRUE,sep="\t")
  cv.glmnet.fit.GIAB<- readRDS("LogisticReg_GIAB_fitted_phasing_error.rds")
  x.test<-as.matrix(sample_info%>%dplyr::select(min_MAF,diff_MAF,log10_distance,r2,d)%>%
                      dplyr::mutate(min_MAF_diff_MAF=min_MAF*diff_MAF)%>%
                      mutate(min_MAF_log10_distance=log10_distance*min_MAF)%>%
                      mutate(min_MAF_r2=min_MAF*r2)%>%
                      mutate(min_MAF_d=min_MAF*d)%>%
                      mutate(diff_MAF_log10_distance=log10_distance*diff_MAF)%>%
                      mutate(diff_MAF_r2=diff_MAF*r2)%>%
                      mutate(diff_MAF_d=diff_MAF*d)%>%
                      mutate(log10_distance_r2=log10_distance*r2)%>%
                      mutate(log10_distance_d=log10_distance*d)%>%
                      mutate(r2_d=r2*d)%>%
                      mutate(min_MAF_diff_MAF_log10_distance=min_MAF*diff_MAF*log10_distance)%>%
                      mutate(min_MAF_diff_MAF_r2=min_MAF*diff_MAF*r2)%>%
                      mutate(min_MAF_diff_MAF_d=min_MAF*diff_MAF*d)%>%
                      mutate(min_MAF_log10_distance_r2=min_MAF*log10_distance*r2)%>%
                      mutate(min_MAF_log10_distance_d=min_MAF*log10_distance*d)%>%
                      mutate(min_MAF_r2_d=min_MAF*r2*d)%>%
                      mutate(diff_MAF_log10_distance_r2=diff_MAF*log10_distance*r2)%>%
                      mutate(diff_MAF_log10_distance_d=diff_MAF*log10_distance*d)%>%
                      mutate(log10_distance_r2_d=log10_distance*r2*d)%>%
                      mutate(min_MAF_diff_MAF_log10_distance_r2=min_MAF*diff_MAF*log10_distance*r2)%>%
                      mutate(min_MAF_diff_MAF_log10_distance_d=min_MAF*diff_MAF*log10_distance*d)%>%
                      mutate(min_MAF_diff_MAF_r2_d=min_MAF*diff_MAF*r2*d)%>%
                      mutate(min_MAF_log10_distance_r2_d=min_MAF*log10_distance*r2*d)%>%
                      mutate(min_MAF_diff_MAF_log10_distance_r2_d=min_MAF*diff_MAF*log10_distance*r2*d))
  sample_info$pred_error_GIAB <- as.numeric(predict(cv.glmnet.fit.GIAB,alpha=0.7163,lambda=0.000271232500776062,newx=x.test,type='response'))

  sample_info<-sample_info%>%group_by(geneID)%>%arrange(chr,pos)%>%dplyr::mutate(pred_error_GIAB=ifelse(pos==dplyr::first(pos),NA,pred_error_GIAB))%>%ungroup()

  write.table(sample_info,meta_error,sep = "\t", row.names = FALSE, col.names = TRUE)
  print(paste0("phasing error sample information saved to ",meta_error,sep=""))

}else{
  print("phasing error prediction file exists!")
}
