#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import pandas as pd
import re
import pickle
import numpy as np
import os
import logging

def generate_modelCount(filename): 
    base_out=os.path.splitext(filename)[0]
    out_modelinput='{0}_modelinput.tsv'.format(base_out)
    gene_df=pd.read_csv(filename,sep="\t",header=0,index_col=False)
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
            pred_prob=[-1 if x!=x else x for x in pred_prob]
            PI_pred=[]
            if nhet==1:
                updated_line += '%s\t%d\n'%(line,0)
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


def significant_genes(prefix,out,modeloutput_dir,outname1,outname2,cutoff,hetSNP_intersect_unique_lambdaPredicted_file):
    outfilename=out+"/"+prefix+"_ASE_all.tsv"
    outfilename_ase=out+"/"+prefix+"_ASE_genes.tsv"
    if (os.path.isfile(outfilename)) and (os.path.isfile(outfilename_ase)):
        logging.info('.... Already Processed {0}'.format(outfilename))
        logging.info('.... Already Processed {0}'.format(outfilename_ase))
    else:
        data_beastie_medtheta = pickle.load(
            open(
                modeloutput_dir
                + "model_theta_med/"
                + outname1,
                "rb",
            )
        )
        data_beastie_maxtail = pickle.load(
            open(
                modeloutput_dir
                + "model_prob/"
                + outname1,
                "rb",
            )
        )

        data_beastie_predictedlambda = pickle.load(
            open(
                modeloutput_dir
                + "model_prob_sum_lambda_predicted/"
                + outname2,
                "rb",
            )
        )
        #
        data_modeloutput = pd.read_csv(hetSNP_intersect_unique_lambdaPredicted_file,
            header=None,
            sep="\t",
        )
        data_modeloutput.columns = [
            "gene_ID",
            "median_altratio",
            "num_hets",
            "totalRef",
            "totalAlt",
            "total_reads",
            "predicted_lambda",
        ]
        data_modeloutput.columns=['gene_ID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda']
        data_modeloutput['BEASTIE_med_theta']=data_beastie_medtheta
        data_modeloutput['BEASTIE_maxtail']=data_beastie_maxtail
        data_modeloutput['BEASTIE_sumtail_lambda_pred']=data_beastie_predictedlambda

        logging.info('ASE gene cut off is {0}'.format(cutoff))

        df_sub = data_modeloutput[['BEASTIE_sumtail_lambda_pred']]
        columns=list(df_sub.columns.values)
        for index, method in enumerate(columns):
            data = df_sub.iloc[:, index]
            if index == 0:
                ncount = data_modeloutput[
                    data_modeloutput["BEASTIE_sumtail_lambda_pred"] > 0.5
                ].count()[8]
            else:
                ncount = data_modeloutput[
                    data_modeloutput["BEASTIE_sumtail_lambda_pred"] > 0.5
                ].count()[8]
            logging.info('Num significant genes @ {}: {}'.format(method,ncount))
        data_modeloutput["ASE"] = data_modeloutput["BEASTIE_sumtail_lambda_pred"]
        data_modeloutput = data_modeloutput.assign(
            ASE=data_modeloutput.apply(ASE_judge, axis=1)
        )

        data_modeloutput.to_csv(outfilename,sep="\t",header=True)
        data_modeloutput_ase=data_modeloutput[data_modeloutput['ASE']=='Y']
        data_modeloutput_ase.to_csv(outfilename_ase,sep="\t",header=True)

def ASE_judge(x):
    if (x['BEASTIE_sumtail_lambda_pred']>0.5):
        return 'Y'
    else:
        return 'N'