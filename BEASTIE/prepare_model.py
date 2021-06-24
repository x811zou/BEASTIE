#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================
import pandas as pd
import re
import pickle


def count_reads(fields):
    if len(fields) >= 4:
        base_thetas = []
        Mreps = int(fields[1])
        total_A = 0
        total_R = 0
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])
            R = float(fields[3 + (rep) * 2])
            base = (A + 1) / (R + 1)
            base_thetas.append(base)
            total_A += A
            total_R += R
        total_reads = total_A + total_R
    return total_reads


# def wrap_for_lambda_prediction(out,sample):
#     nhets_list=[]
#     total_reaDepth=[]
#     IDs = []
#     inFile=out+sample+"_parsed_mpileup_exonicHetSnps_model_input.txt"
#     with open(inFile, "r") as IN:
#         for line in IN:             # start reading in my pileup results line by line
#             #print(line)
#             fields=line.rstrip().split()
#             ID=fields[0]
#             IDs.append(ID)
#             nhets=fields[1]
#             total_reads=count_reads(fields)
#             total_reaDepth.append(total_reads)
#             nhets_list.append(nhets)
#     dict = {'gene_sample': IDs,'total_reads':total_reaDepth,'num_hets':nhets_list}
#     #print(dict)
#     df = pd.DataFrame(dict)
#     outFile=inFile.replace('.txt','')+'.forLambdaPred.csv'
#     df.to_csv(outFile, header=False, index=False)
#     print("convert model input into a format for lambda prediction size %s saved to %s"%(len(outFile),outFile))


def update_model_input_lambda_phasing(out, prefix, pred_prob_column):
    pred_prob_column = "pred_error_GIAB"
    outfile = (
        out
        + str(prefix)
        + "_parsed_mpileup_exonicHetSnps_model_input_w_phasingError.txt"
    )
    # input:  NA12878_parsed_mpileup_exonicHetSnps_model_input.txt; NA12878_logisticRegression_input_phasingErrorpredicted.tsv
    # output: NA12878_parsed_mpileup_exonicHetSnps_model_input.phasingProb.txt
    model_input = out + str(prefix) + "_parsed_mpileup_exonicHetSnps_model_input.txt"
    phasing_error = pd.read_csv(
        out + str(prefix) + "_logisticRegression_input_phasingErrorpredicted.tsv",
        header=0,
        sep="\t",
    )
    # ENSG00000177757.1	3	1	0	1	0	1	1 prob1 prob2
    # make sure these two data sets are consistent = nrow(model_input) = num of genes in phasing_error
    updated_line = ""
    # obtain a overlapped dataframe between these two data input: df
    # N=len(model_input)[0] # num of genes
    counter = 0
    with open(model_input, "r") as stream_in:
        for i, line in enumerate(
            stream_in
        ):  # start reading in my pileup results line by line
            line = re.sub("\n", "", line)
            counter += 1
            if counter % 1000 == 0:
                print("%d processed." % counter)
                # break
            # ENSG00000099331	4	184	225	317	355	351	301	294	304
            geneID = line.split("\t")[0]
            nhet = int(line.split("\t")[1])
            # print("Number of hets %s"%(nhet))
            phasing_error_selected = phasing_error[(phasing_error["geneID"] == geneID)]
            pred_prob = phasing_error_selected[str(pred_prob_column)].tolist()
            PI_pred = []
            # updated_line+=line
            # print(updated_line)
            if nhet == 1:
                updated_line += "%s\t%d\n" % (line, -1)
            else:
                for m in range(int(nhet) - 1):  # 0,1,2
                    # print("m : %s"%(m))
                    if m == 0:
                        PI_pred = "%s" % (pred_prob[m])
                    else:
                        PI_pred = "%s\t%s" % (PI_pred, pred_prob[m])
                # print(updated_line)
                # print(PI_pred)
                updated_line += "%s\t%s\n" % (line, PI_pred)
                # print(updated_line)
    # print(updated_line)
    file1 = open(outfile, "w")
    file1.write(updated_line)
    file1.close()
    print("saved to %s" % (outfile))
    # return updated_line


def significant_genes(prefix, out, cutoff):
    model_dir = out
    modeloutput_dir = model_dir + "output_pkl/"
    data_beastie_medtheta = pickle.load(
        open(
            modeloutput_dir
            + "model_theta_med/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5.pickle",
            "rb",
        )
    )
    data_beastie_maxtail = pickle.load(
        open(
            modeloutput_dir
            + "model_prob/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5.pickle",
            "rb",
        )
    )
    data_beastie_lambda12 = pickle.load(
        open(
            modeloutput_dir
            + "model_prob_sum_lambda1.2/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5.pickle",
            "rb",
        )
    )
    data_beastie_lambda14 = pickle.load(
        open(
            modeloutput_dir
            + "model_prob_sum_lambda1.4/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5.pickle",
            "rb",
        )
    )
    data_beastie_lambda16 = pickle.load(
        open(
            modeloutput_dir
            + "model_prob_sum_lambda1.6/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5.pickle",
            "rb",
        )
    )
    data_beastie_predictedlambda = pickle.load(
        open(
            modeloutput_dir
            + "model_prob_sum_lambda_predicted/"
            + str(prefix)
            + "_parsed_mpileup_exonicHetSnps_model_input_s-0.5_a-0.05.pickle",
            "rb",
        )
    )
    #
    data_modeloutput = pd.read_csv(
        model_dir
        + prefix
        + "_parsed_mpileup_overlapped_exonicHetSnps_Lambda_0.05_prediced.tsv",
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
    data_modeloutput["BEASTIE_med_theta"] = data_beastie_medtheta
    data_modeloutput["BEASTIE_maxtail"] = data_beastie_maxtail
    data_modeloutput["BEASTIE_sumtail_lambda_12"] = data_beastie_lambda12
    data_modeloutput["BEASTIE_sumtail_lambda_14"] = data_beastie_lambda14
    data_modeloutput["BEASTIE_sumtail_lambda_16"] = data_beastie_lambda16
    data_modeloutput["BEASTIE_sumtail_lambda_pred"] = data_beastie_predictedlambda

    sigGenes = {}
    print("Cutoff : %s" % (cutoff))
    df_sub = data_modeloutput[["BEASTIE_sumtail_lambda_pred"]]
    columns = list(df_sub.columns.values)
    data_modeloutput = data_modeloutput.drop(
        [
            "BEASTIE_maxtail",
            "BEASTIE_sumtail_lambda_12",
            "BEASTIE_sumtail_lambda_14",
            "BEASTIE_sumtail_lambda_16",
        ],
        axis=1,
    )
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
        print("Num significant genes @ {}: {}".format(method, ncount))
        # sigGenes[method] = data[(data<=cutoff)].count()
    data_modeloutput["ASE"] = data_modeloutput["BEASTIE_sumtail_lambda_pred"]
    data_modeloutput = data_modeloutput.assign(
        ASE=data_modeloutput.apply(ASE_judge, axis=1)
    )
    outfilename = model_dir + prefix + "_ASE_all.tsv"
    data_modeloutput.to_csv(outfilename, sep="\t", header=True)
    data_modeloutput_ase = data_modeloutput[data_modeloutput["ASE"] == "Y"]
    outfilename_ase = model_dir + prefix + "_ASE_genes.tsv"
    data_modeloutput_ase.to_csv(outfilename_ase, sep="\t", header=True)
    # return sigGenes#,data_modeloutput#


def ASE_judge(x):
    if x["BEASTIE_sumtail_lambda_pred"] > 0.5:
        return "Y"
    else:
        return "N"
