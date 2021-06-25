import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import operator
import scipy
from scipy.stats import f_oneway
import statistics as stat
import seaborn as sns

def coeff_analysis(numerical1, numerical2):
    gene_dict = dict()
    df = pd.read_csv("operon_windows2.csv")

    df = df[df[numerical1].notna()]
    df = df[df[numerical2].notna()]

    for index, row in df.iterrows():
        if row["r_squaredCDS"] > 0.75:    ## Change to 0.7?
            try:
                if row["ORF"][0:2] != "nc":
                    if row["ORF"][0:4] != "MPNt":
                        if row["ORF"] not in gene_dict:
                            gene_dict[row["ORF"]] = [[],[]]
                        if row["r_squared_local"] > 0.7:   ## Change?
                            gene_dict[row["ORF"]][0].append(row[numerical1])   #
                            gene_dict[row["ORF"]][1].append(row[numerical2])   #
            except:
                pass

    coeff_dict = dict()
    coeff_list_sig = []
    coeff_list_notsig = []
    pforty_list = []
    pthirty_list = []
    psixty_list = []
    peighty_list = []
    nforty_list = []
    nthirty_list = []
    nsixty_list = []
    neighty_list = []
    for gene, values in gene_dict.items():
        if len(values[0]) > 10:    ## genes with at least 10 windows
            numlist_1 = values[0]
            numlist_2 = values[1]
            results = (scipy.stats.pearsonr(numlist_1, numlist_2))
            if results[1]:
                coeff_dict[gene] = results[0]
                coeff_list_notsig.append(results[0])
                #print("GC_content") #
                #print("p-value:    " + str(results[1]))    
                #print("Pearson's: " + str(results[0]))
                #print("N:         " + str(len(numlist_1)))
                #print("\n\n")
                #plt.scatter(numlist_1,numlist_2)
                #plt.show()
                #if gene == "MPN681":
                #    print(results[0], results[1])
            if results[1] < 0.05:
                coeff_list_sig.append(results[0]) 
                if results[0] > 0.3:
                    pthirty_list.append([gene, results[0]])
                if results[0] > 0.4:
                    pforty_list.append([gene, results[0]])
                if results[0] > 0.6:
                    psixty_list.append([gene, results[0]])
                if results[0] > 0.8:
                    peighty_list.append([gene, results[0]])
                    #plt.scatter(numlist_1,numlist_2)
                    #plt.show()
                if results[0] < -0.3:
                    nthirty_list.append([gene, results[0]])
                if results[0] < -0.4:
                    nforty_list.append([gene, results[0]])
                if results[0] < -0.6:
                    nsixty_list.append([gene, results[0]])
                if results[0] < -0.8:
                    neighty_list.append([gene, results[0]])
                    #plt.scatter(numlist_1,numlist_2)
                    #plt.show()

                #plt.title("MPN664  C: " + str(results[0]))
                #plt.show()

    print(numerical1)  #
    print("Sig. Mean C:\t" + str(stat.mean(coeff_list_sig)))
    print("Not Sig. Mean C:\t" + str(stat.mean(coeff_list_notsig)))
    print("N total:\t"                + str(len(gene_dict)))
    print("N total Sig:\t"            + str(len(coeff_list_sig)))
    print("N over 0.3:\t"            + str(len(pthirty_list)) + ". (" + str(len(pthirty_list)/len(gene_dict)*100) + "%)")
    print("N over 0.4:\t"            + str(len(pforty_list)) + ". (" + str(len(pforty_list)/len(gene_dict)*100) + "%)")
    print("N over 0.6:\t"            + str(len(psixty_list)) + ". (" + str(len(psixty_list)/len(gene_dict)*100) + "%)")
    print("N over 0.8:\t"            + str(len(peighty_list)) + ". (" + str(len(peighty_list)/len(gene_dict)*100) + "%)")
    print("N under -0.3:\t"          + str(len(nthirty_list)) + ". (" + str(len(nthirty_list)/len(gene_dict)*100) + "%)")
    print("N under -0.4:\t"          + str(len(nforty_list)) + ". (" + str(len(nforty_list)/len(gene_dict)*100) + "%)")
    print("N under -0.6:\t"          + str(len(nsixty_list)) + ". (" + str(len(nsixty_list)/len(gene_dict)*100) + "%)")
    print("N under -0.8:\t"          + str(len(neighty_list)) + ". (" + str(len(neighty_list)/len(gene_dict)*100) + "%)")
    peighty_list.sort(key = lambda x: x[1], reverse = True)
    neighty_list.sort(key = lambda x: x[1])
    ptop_genes = [i[0] for i in peighty_list]
    ntop_genes = [i[0] for i in neighty_list]
    print("Top positive Coeff: "    + str(ptop_genes))
    print("Top negative Coeff: "    + str(ntop_genes))
    print("===========================================\n\n")


    ## PLOT
    sns.set_theme(style="whitegrid")
    data =  {"Pearson's coefficient": [coeff_list_sig]}
    df = pd.DataFrame(coeff_list_sig, columns = ["Pearson's coefficient"])  
    sns.histplot(data=df, x="Pearson's coefficient", bins = 20)
    #plt.show()

#coeff_analysis("coverage_local", "half_life_local")   # SINGLE-ANALYSIS
# ALL-FEATURES ANALYSIS
#features = ["half_life_local", "GC_content_local", "energy_local", "coverage_local", "ribo-seq_RP4_local", "ribo-seq_both_local", "rnase3_downreg_local", "CAI_index_local"]
#for feature in features:
#    coeff_analysis(feature, "half_life_local")

features = ["half_life_local", "GC_content_local", "energy_local", "coverage_local", "ribo-seq_RP4_local", "ribo-seq_both_local", "rnase3_downreg_local"]
for feature in features:
    coeff_analysis(feature, "half_life_local")


# Importnat parameters: 
#    Quality r-squared (0.7)
#    Min number of windows of gene (10)
