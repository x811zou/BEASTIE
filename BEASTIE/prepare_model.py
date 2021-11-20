#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================

import re
import os
import logging
import scipy.stats
import pandas as pd
import numpy as np
from math import floor, log10,pi

def generate_modelCount(filename):
    logging.info('filename: {0}'.format(filename))
    gene_df=pd.read_csv(filename,sep="\t",header=0,index_col=False)
    base_out=os.path.splitext(filename)[0]
    out_modelinput='{0}.modelinput.tsv'.format(base_out)
    logging.info('out_modelinput: {0}'.format(out_modelinput))
    #gene_df=pd.read_csv(filename,sep="\t",header=0,index_col=False)
    if not os.path.isfile(out_modelinput):
        unique_gene = gene_df['geneID'].unique()
        rst = ""
        for each in unique_gene:
            idx = each
            gene_lst = [idx, sum(gene_df["geneID"] == each)] + list(np.ravel(gene_df[gene_df["geneID"] == each][["altCount","refCount"]].values))
            rst += "\t".join([str(x) for x in gene_lst]) +"\n"
        file1 = open(out_modelinput,"w")
        file1.write(rst)
        file1.close()
        logging.info('..... model input saved to {0}'.format(out_modelinput))
    else:
        logging.info('..... model input exists at {0}'.format(out_modelinput))

def count_reads(fields):
    if(len(fields)>=4):
        base_thetas=[]
        Mreps=int(fields[1])
        total_A=0
        total_R=0
        for rep in range(Mreps):
            A = float(fields[2+rep*2])
            R = float(fields[3+(rep)*2])
            base = (A+1)/(R+1)
            base_thetas.append(base)
            total_A+=A
            total_R+=R
        total_reads=total_A+total_R
    return total_reads

def update_model_input_lambda_phasing(pred_prob_column,base_modelin,base_modelin_error,meta_error):
    pred_prob_column="pred_error_GIAB"
    outfile = base_modelin_error
    model_input=base_modelin
    phasing_error = pd.read_csv(meta_error,header=0,sep="\t")
    updated_line = ""
    counter = 0
    with open(model_input, "r") as stream_in:
        for i, line in enumerate(stream_in):             # start reading in my pileup results line by line
            line=re.sub('\n','',line)
            counter += 1
            if counter % 1000 == 0:
                logging.info('{0} line processed'.format(counter))
            geneID = line.split("\t")[0]
            nhet= int(line.split("\t")[1])
            phasing_error_selected = phasing_error[(phasing_error["geneID"]==geneID)]
            pred_prob=phasing_error_selected[str(pred_prob_column)].tolist()
            if len(pred_prob)!=0:
                del pred_prob[0]
            nNAs=np.count_nonzero(np.isnan(pred_prob))
            pred_prob=[-1 if x!=x else x for x in pred_prob] # -1 if missing
            PI_pred=[]
            if nhet==1:
                updated_line += '%s\t%d\n'%(line,0) # 0 for 1 het gene
            else:
                for m in range(int(nhet)-1): #0,1,2
                    if m==0:
                        PI_pred="%d\t%s" % (nNAs,pred_prob[m])
                    else:
                        PI_pred="%s\t%s" % (PI_pred,pred_prob[m])
                updated_line += '%s\t%s\n' % (line,PI_pred)
    file1 = open(outfile,"w")
    file1.write(updated_line)
    file1.close()

def significant_genes(df_ibeastie,df_binomial,df_adm,outfilename,outfilename_ase,cutoff,hetSNP_intersect_unique_lambdaPredicted_file):
    data_modeloutput = pd.read_csv(hetSNP_intersect_unique_lambdaPredicted_file,header=0,sep="\t") #data_modeloutput = pd.read_csv("/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr20/output/TEMP/HG00096_chr20_hetSNP_intersect_unique_alpha0.05_lambdaPredicted.tsv",sep="\t")
    # data_modeloutput.columns = [
    #     "geneID",
    #     "median_altratio",
    #     "num_hets",
    #     "totalRef",
    #     "totalAlt",
    #     "total_reads",
    #     "predicted_lambda",
    # ]
    #print(data_modeloutput)
    logging.debug('size of hetSNP_intersect_unique_lambdaPredicted_file file {}'.format(len(data_modeloutput)))
    logging.debug('size of df_beastie {}'.format(len(df_ibeastie)))
    

    df_output = pd.merge(data_modeloutput,df_ibeastie,on=['geneID'], how="inner")
    logging.debug('size of model output is {0} ; size of hetSNP_intersect_unique_lambdaPredicted_file is {1}; intersection size is {2}'.format(len(df_ibeastie),len(data_modeloutput),len(df_output)))
    ncount = df_output[df_output["posterior_mass_support_ALT"] > cutoff].count()[8]
    logging.info('{} genes with ASE out of total genes {} ({}%) at @ {} > ASE cutoff {}'.format(ncount,len(df_output),round((ncount/len(df_output))*100,3),"posterior_mass_support_ALT",cutoff))
    df_output["foldLog2MedSq_over_var"] = (abs(np.log2(df_output["posterior_median"])))**2/df_output["posterior_variance"]
    #df_output["foldLog2MedSq_over_var"]=df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    df_output['median.altRatio']=df_output['median.altRatio'].round(3)
    #df_output["extreme_val"] = df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    extremes = scipy.stats.chi2.ppf(0.95, 1)
    ncount2 = df_output[df_output["foldLog2MedSq_over_var"] > extremes].count()[8]
    logging.debug('{} genes with ASE out of total genes {} ({}%) at @ {} > extreme values: {} = QCHISQ(p=0.95,df=1)'.format(ncount2,len(df_output),round((ncount2/len(data_modeloutput))*100,3),"foldLog2MedSq_over_var",extremes))
    df_output = df_output.assign(
        extreme_val = lambda dataframe: dataframe['foldLog2MedSq_over_var'].map(lambda foldLog2MedSq_over_var: "Y" if foldLog2MedSq_over_var > extremes else "N")
    )
    df_output["ASE"] = df_output["posterior_mass_support_ALT"]
    #print(">> merged V2")
    #print(df_output.head(n=5))
    df_output = df_output.assign(
        ASE = lambda dataframe: dataframe['posterior_mass_support_ALT'].map(lambda posterior_mass_support_ALT: "Y" if posterior_mass_support_ALT > cutoff else "N")
    )
    df_output_bi = pd.merge(df_output,df_binomial,on=['geneID'], how="inner") 

    if df_output_bi.shape[0]==0:
        print(">> merged binomial with beastie is empty")
        df_output_bi=df_output
    else:
        df_output_bi = df_output_bi.drop(['FirstSite_pval','FirstSite_esti','NaiveSum_esti','Pseudo_pval','Pseudo_esti','MajorSite_esti'], axis=1)
    df_output_bi_adm = pd.merge(df_output_bi,df_adm,on=['geneID'], how="inner") 
    #print(">> merged binomial with beastie and ADM")
    if df_output_bi_adm.shape[0]==0:
        print(">> merged binomial with beastie and ADM is empty")
        df_output_bi_adm=df_output
    else:
        df_output_bi_adm = df_output_bi_adm.drop(['ADM_esti'], axis=1)
    #print(df_output_bi_adm.head(n=5))

    logging.info('..... done with dropping')
    df_output=df_output_bi_adm
    #outfilename="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_all.tsv"
    #outfilename_ase="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_cutoff_0.5_filtered.tsv"
    df_output_bi.to_csv(outfilename,sep="\t",header=True,index=False)
    df_output_ase=df_output[df_output["ASE"]=='Y']
    df_output_ase.to_csv(outfilename_ase,sep="\t",header=True,index=False)