#!/usr/bin/env Rscript

# python style read in parameters:
# https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/#

# cmd="Rscript --vanilla predict_lambda_phasingError.R %s %s %s"%(out,prefix,alpha)
# os.system(cmd)
# #input: NA12878_parsed_mpileup_exonicHetSnps_model_input.forLambdaPred.csv
# #output: NA12878_parsed_mpileup_exonicHetSnps_model_input.LambdaPredicted.csv
# #      geneID      |tota_reads| num_hets | predicted lambda
# #"ENSG00000177757.1"     2	        1	   1.3978282039415

args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
  # default output file
  args[2] <- "out.txt"
}

#### call library
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library(glmnetUtils))
#suppressMessages(library(glmnet))

tmp <- paste0(args[1], "/", sep = "")
sample <- args[2]
hetSNP_intersect_unique <- args[3]
meta <- args[4]
meta_error <- args[5]
beastie_wd <- args[6]
phasing_method <- args[7]
source(file.path(beastie_wd, "Get_phasing_error_rate.R"))

################################################ 2 predict phasing error ###################################
if (!file.exists(meta_error) && (phasing_method != "nophasing")) {
  print("meta file with phasing error not exists, has to run logistic reg!")
  sample_info <- read.delim(file.path(meta), header = TRUE, sep = "\t")
  data <- sample_info
  data$chr <- factor(data$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
  data$pos <- as.integer(data$pos)
  data$lag_pos <- lag(data$pos)
  data <- data %>%
    group_by(chr) %>%
    dplyr::mutate(lag_pos = ifelse(pos == dplyr::first(pos), NA, lag_pos)) %>%
    ungroup()
  data <- data %>%
    group_by(geneID) %>%
    dplyr::mutate(lag_pos = ifelse(pos == dplyr::first(pos), NA, lag_pos)) %>%
    ungroup()
  # data<-data%>%group_by(chr,pos)%>%slice(1)
  # data_w_LD<-left_join(processed_exon,data_LD)
  data_final_reginput <- data %>%
    mutate(distance = pos - lag_pos) %>%
    mutate(log10_distance = log10(distance)) %>%
    mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF)) # the minimum AF for the current variant

  data_final_reginput$lag_MAF <- lag(data_final_reginput$MAF)

  data_final_reginput <- data_final_reginput %>% # MAF for the previous variant
    mutate(min_MAF = ifelse(MAF > lag_MAF, lag_MAF, MAF)) %>% # the mininum between the pair
    mutate(diff_MAF = abs(MAF - lag_MAF)) # the difference between the pair

  data_final_reginput <- data_final_reginput %>%
    group_by(geneID) %>%
    arrange(chr, pos) %>%
    mutate(diff_pos = pos - lag_pos) %>%
    ungroup()

  sample_info <- data_final_reginput %>% select(geneID, chr, pos, pair_pos, lag_pos, rsid, log10_distance, d, r2, MAF, lag_MAF, min_MAF, diff_MAF)
  if (phasing_method == "VCF") {
    print("VCF phasing is used! we are not predicting phasing error !")
    sample_info$pred_error_GIAB <- NA
  } else {
    print("shapeit2 phasing is used! we are predicting phasing error !")
    cv.glmnet.fit.GIAB <- readRDS(file.path(beastie_wd, "LogisticReg_GIAB_fitted_phasing_error.rds"))
    x.test <- as.matrix(sample_info %>% dplyr::select(min_MAF, diff_MAF, log10_distance, r2, d) %>%
      dplyr::mutate(min_MAF_diff_MAF = min_MAF * diff_MAF) %>%
      mutate(min_MAF_log10_distance = log10_distance * min_MAF) %>%
      mutate(min_MAF_r2 = min_MAF * r2) %>%
      mutate(min_MAF_d = min_MAF * d) %>%
      mutate(diff_MAF_log10_distance = log10_distance * diff_MAF) %>%
      mutate(diff_MAF_r2 = diff_MAF * r2) %>%
      mutate(diff_MAF_d = diff_MAF * d) %>%
      mutate(log10_distance_r2 = log10_distance * r2) %>%
      mutate(log10_distance_d = log10_distance * d) %>%
      mutate(r2_d = r2 * d) %>%
      mutate(min_MAF_diff_MAF_log10_distance = min_MAF * diff_MAF * log10_distance) %>%
      mutate(min_MAF_diff_MAF_r2 = min_MAF * diff_MAF * r2) %>%
      mutate(min_MAF_diff_MAF_d = min_MAF * diff_MAF * d) %>%
      mutate(min_MAF_log10_distance_r2 = min_MAF * log10_distance * r2) %>%
      mutate(min_MAF_log10_distance_d = min_MAF * log10_distance * d) %>%
      mutate(min_MAF_r2_d = min_MAF * r2 * d) %>%
      mutate(diff_MAF_log10_distance_r2 = diff_MAF * log10_distance * r2) %>%
      mutate(diff_MAF_log10_distance_d = diff_MAF * log10_distance * d) %>%
      mutate(log10_distance_r2_d = log10_distance * r2 * d) %>%
      mutate(min_MAF_diff_MAF_log10_distance_r2 = min_MAF * diff_MAF * log10_distance * r2) %>%
      mutate(min_MAF_diff_MAF_log10_distance_d = min_MAF * diff_MAF * log10_distance * d) %>%
      mutate(min_MAF_diff_MAF_r2_d = min_MAF * diff_MAF * r2 * d) %>%
      mutate(min_MAF_log10_distance_r2_d = min_MAF * log10_distance * r2 * d) %>%
      mutate(min_MAF_diff_MAF_log10_distance_r2_d = min_MAF * diff_MAF * log10_distance * r2 * d))
    sample_info$pred_error_GIAB <- as.numeric(predict(cv.glmnet.fit.GIAB, alpha = 0.7163, lambda = 0.000271232500776062, newx = x.test, type = "response"))
    sample_info <- sample_info %>%
      group_by(geneID) %>%
      arrange(chr, pos) %>%
      dplyr::mutate(pred_error_GIAB = ifelse(pos == dplyr::first(pos), NA, pred_error_GIAB)) %>%
      ungroup()
  }
  write.table(sample_info, meta_error, sep = "\t", row.names = FALSE, col.names = TRUE)
  print(paste0("phasing error sample information saved to ", meta_error, sep = ""))
} else {
  print("phasing error prediction file exists!")
}
