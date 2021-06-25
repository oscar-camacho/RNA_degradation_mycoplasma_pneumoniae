import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import operator
import scipy
from scipy.stats import f_oneway

def correlation(numerical1, numerical2):
    df = pd.read_csv("transcript_annotation.csv",  encoding = "unicode_escape")
    numlist_1 = []
    numlist_2 = []

    df = df[df[numerical1].notna()]
    df = df[df[numerical2].notna()]

    for index, row in df.iterrows():
        if row["half_lifeCDS"] > 0:
            #if row["ORF"][0:2] != "nc":
             #   if row["ORF"][0:4] != "MPNt":
                    if row["r_squaredCDS"] > 0.75:
                        #if row["r_squared_local"] > 0.70:   ## Change when doing gene_level
                            #if row["Classification"] == "1 CDS. 1 TSS":
                                numlist_1.append(row[numerical1])
                                numlist_2.append(row[numerical2])


    #numlist_1 = np.log10(numlist_1)


    results = (scipy.stats.pearsonr(numlist_1, numlist_2))
    if results[1] < 0.05:
        print(numerical1)
        print("p-value:    " + str(results[1]))    
        print("Pearson's: " + str(results[0]))
        print("N:         " + str(len(numlist_1)))
        print("\n\n")
        plt.scatter(numlist_1,numlist_2)
        plt.xlabel("Half life (min) - Including all time points")
        plt.ylabel("Half life (min) - Not including time point 0 and 2")
        plt.axis('square')

        plt.show()
        plt.savefig("nada.png", dpi=100)



#features = ["log2_foldchange(three/five)", "log2_foldchange(maximum/minimum)", "glucose_starvation", "rnase3_downreg", "gene_number_operon", "gene_level_avg_ribo_density_scaled", "half_lifeCDS", "first_riboseq", "last_riboseq", "first_riboseq1", "last_riboseq1", "CDS_length", "transcript_length", "fiveUTR_length", "threeUTR_length", "GC_content", "protein_copy_number", "protein_halflife", "CDS_energy", "fiveUTR_energy", "threeUTR_energy", "fiveCDS_energy", "threeCDS_energy", "first_energy3", "last_energy3", "first_energy5", "last_energy5", "RNA_copy_number_6h", "coverage", "CAI_index", "fiveUTR_energy?", "RBS_internal", "RBS_upstream"]
#features = ["rnase3_downreg", "first_riboseq1", "first_riboseq2", "last_riboseq1", "last_riboseq2", "CDS_energy", "fiveCDS_energy", "threeCDS_energy", "first_energy1", "last_energy1", "first_energy3", "last_energy3", "first_energy5", "last_energy5"]

#features = ["normhalf_life", "normGC_content", "normenergy", "normcoverage", "normribo-seq_RP4", "normribo-seq_both", "normrnase3_downreg", "normCAI_index"]
#features = ["half_life_local", "GC_content_local", "energy_local", "coverage_local", "ribo-seq_RP4_local", "ribo-seq_both_local", "rnase3_downreg_local", "CAI_index_local"]

#for feature in features:
    #correlation(feature, "half_lifeCDS")

correlation("half_lifeCDS_menos0", "half_lifeCDS")


