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
suppressMessages(library("DESeq2"))
suppressMessages(library("pasilla"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("LDlinkR"))
library(glmnetUtils)
mytoken="c313799c13c3"
#suppressMessages(library(foreach))
#suppressMessages(library(doParallel))
#suppressMessages(library(parallel))
#### read in self-defined functions
source("/Users/scarlett/Box/Allen_lab/Backup/R/Get_phasing_error_rate.R")
source("/Users/scarlett/Box/Allen_lab/Backup/R/Get_LD.R")
#### read in parameters
alpha=as.numeric(args[1])
out=args[2]
sample=args[3]
model=args[4]

############################################## 1. predict lambda #####################################################
lambda.fit.simulation<- readRDS("/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/LinearReg_BEASTIE_fitted_model_lambda.rds")
predict_lambda_realdata <- function(alpha,local_dir,sample,in_data){
  in_data<-read.delim(file.path(paste0(out,sample,"_parsed_mpileup_overlapped_exonicHetSnps_forLambda.tsv",sep="")),header=TRUE,sep="\t")
  out_data<-paste0(out,sample,"_parsed_mpileup_overlapped_exonicHetSnps_","Lambda_",as.character(as.numeric(alpha)),"_prediced.tsv",sep="")  
  #colnames(in_data)<-c("gene_ID","total_reads","num_hets")
  data<-in_data%>%dplyr::mutate(predicted_lambda=(log(alpha/(1-alpha)) -(as.numeric(lambda.fit.simulation$coefficients[1])+as.numeric(lambda.fit.simulation$coefficients[3])*as.integer(totalCount)))/as.numeric(lambda.fit.simulation$coefficients[2]))
  #data<-in_data%>%mutate(predicted_lambda=(log(alpha/(1-alpha)) -(15.587909+-0.006483*total_reads))/-13.248682)
  write.table(data,file = out_data,row.names=FALSE,col.names = FALSE,sep="\t")
  print(paste0("model input with predicted lambda saved to ",out_data,sep=""))
  return(data)
}

predicted_df=predict_lambda_realdata(alpha,local_dir,sample,lambda_pred)

################################################ 2.1 predict phasing error ###################################
sample_info<-read.delim(file.path(paste0(out,sample,"_parsed_mpileup_overlapped_exonicHetSnps.tsv",sep="")),header=TRUE,sep="\t")
cv.glmnet.fit.GIAB<- readRDS("/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/LogisticReg_GIAB_fitted_phasing_error.rds")
x.test<-as.matrix(sample_info%>%dplyr::select(min_EUR_MAF,diff_EUR_MAF,log10_distance,r2,d)%>%
                     dplyr::mutate(min_EUR_MAF_diff_EUR_MAF=min_EUR_MAF*diff_EUR_MAF)%>%
                     mutate(min_EUR_MAF_log10_distance=log10_distance*min_EUR_MAF)%>%
                     mutate(min_EUR_MAF_r2=min_EUR_MAF*r2)%>%
                     mutate(min_EUR_MAF_d=min_EUR_MAF*d)%>%
                     mutate(diff_EUR_MAF_log10_distance=log10_distance*diff_EUR_MAF)%>%
                     mutate(diff_EUR_MAF_r2=diff_EUR_MAF*r2)%>%
                     mutate(diff_EUR_MAF_d=diff_EUR_MAF*d)%>%
                     mutate(log10_distance_r2=log10_distance*r2)%>%
                     mutate(log10_distance_d=log10_distance*d)%>%
                     mutate(r2_d=r2*d)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_log10_distance=min_EUR_MAF*diff_EUR_MAF*log10_distance)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_r2=min_EUR_MAF*diff_EUR_MAF*r2)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_d=min_EUR_MAF*diff_EUR_MAF*d)%>%
                     mutate(min_EUR_MAF_log10_distance_r2=min_EUR_MAF*log10_distance*r2)%>%
                     mutate(min_EUR_MAF_log10_distance_d=min_EUR_MAF*log10_distance*d)%>%
                     mutate(min_EUR_MAF_r2_d=min_EUR_MAF*r2*d)%>%
                     mutate(diff_EUR_MAF_log10_distance_r2=diff_EUR_MAF*log10_distance*r2)%>%
                     mutate(diff_EUR_MAF_log10_distance_d=diff_EUR_MAF*log10_distance*d)%>%
                     mutate(log10_distance_r2_d=log10_distance*r2*d)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_log10_distance_r2=min_EUR_MAF*diff_EUR_MAF*log10_distance*r2)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_log10_distance_d=min_EUR_MAF*diff_EUR_MAF*log10_distance*d)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_r2_d=min_EUR_MAF*diff_EUR_MAF*r2*d)%>%
                     mutate(min_EUR_MAF_log10_distance_r2_d=min_EUR_MAF*log10_distance*r2*d)%>%
                     mutate(min_EUR_MAF_diff_EUR_MAF_log10_distance_r2_d=min_EUR_MAF*diff_EUR_MAF*log10_distance*r2*d))
sample_info$pred_error_GIAB <- as.numeric(predict(cv.glmnet.fit.GIAB,alpha=0.7163,lambda=0.0002712325,newx=x.test,type='response'))
#set.seed(123)               
random<-runif(n=dim(sample_info)[1])
sample_info$random<-random
sample_info<-sample_info%>%mutate(pred_label_GIAB=ifelse(random<=pred_error_GIAB,1,0))
# -1 for those NAs
out_data=paste0(out,sample,"_logisticRegression_input_phasingErrorpredicted.tsv",sep="")
write.table(sample_info,out_data,sep = "\t", row.names = FALSE, col.names = TRUE)
print(paste0("phasing error sample information saved to ",out_data,sep=""))
