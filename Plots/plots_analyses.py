import pandas as pd   
import statistics as stats
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import ttest_ind
from numpy import mean
from numpy import std
from scipy.stats import sem
from scipy.stats import t
from math import sqrt
from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import operator
import scipy
from scipy.stats import f_oneway
import statistics as stat




def gene_half_life_distributions():
    trna = []
    mrna = []
    ncrna = []
    threeUTR = []
    fiveUTR = []
    df =pd.read_csv("transcript_annotation.csv")
    for index, row in df.iterrows():
        if row["r_squaredCDS"] > 0.8:
            if row["half_lifeCDS"] > 1 and row["half_lifeCDS"] < 30:
                if row["ORF"][0:2] == "nc":
                    ncrna.append(row["half_lifeCDS"])
                if row["ORF"][0:4] == "MPNt":
                    trna.append(row["half_lifeCDS"])
                if row["ORF"][0:2] != "nc" and row["ORF"][0:4] != "MPNt":
                    mrna.append(row["half_lifeCDS"])


    df = pd.read_csv("operon_windows.csv")
    for index, row in df.iterrows():
        if row["r_squared_local"] > 0.7:
            if row["half_life_local"] > 1:
                if row["Region"] == "5'UTR":
                    fiveUTR.append(row["half_life_local"])
                if row["Region"] == "3'UTR":
                    if row["half_life_local"] > 4.7:
                        threeUTR.append(row["half_life_local"])
        


    box_plot_data = [fiveUTR, mrna, threeUTR, ncrna, trna]
    #plt.boxplot(box_plot_data,patch_artist=True,labels=["5'UTR","mRNA","3'UTR","ncRNA","tRNA"])
    sns.set_theme(style="whitegrid")
    lst = box_plot_data
    df = pd.DataFrame(lst, index = ["mRNA\n(5'UTR)","mRNA\n(CDS)", "mRNA\n(3'UTR)", "ncRNA", "tRNA"])
    df = df.T

    sns.boxplot(data=df.loc[:, ["mRNA\n(5'UTR)","mRNA\n(CDS)", "mRNA\n(3'UTR)", "ncRNA", "tRNA"]])

    x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h, col = df['mRNA\n(CDS)'].max() + 0.5, 0.5, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 0, 2   
    y, h, col = df["mRNA\n(CDS)"].max() + 4, 0.5, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 1, 2   
    y, h, col = df["mRNA\n(CDS)"].max() + 2, 0.5, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 3, 4   
    y, h, col = df["tRNA"].max() + 2, 0.5, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 1, 4   
    y, h, col = df["mRNA\n(CDS)"].max() + 6.5, 0.5, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    plt.ylabel("Half-life (min)")

    plt.show()



def genes_in_operons_increase():

    operon_dict = dict()
    df = pd.read_csv("transcript_annotation.csv")
    for index, row in df.iterrows():
        if row["Classification"] != "1 CDS. 1 TSS":
            if row["Strand"] == "-":
                if row["ORF"][0:2] == "nc":
                    continue
            if row["OPERON"] not in operon_dict:
                operon_dict[row["OPERON"]] = []

            if row["ORF"][0:4] == "MPNt":
                operon_dict[row["OPERON"]].append([0, 0])
            else:
                operon_dict[row["OPERON"]].append([row["half_lifeCDS"], row["r_squaredCDS"]])


    order_dict = dict()
    for k, data in operon_dict.items():
        if data[0][1] > 0.7:
            reference = data[0][0]
            i = 0
            for pairs in data:
                i += 1
                if i not in order_dict:
                    order_dict[i] = []
                if pairs[1] > 0.7:
                    foldchange = pairs[0]/reference
                    order_dict[i].append(foldchange)
    gene2bis = []
    gene2 = order_dict[2]
    for i in gene2:
        if i > 0.7:
            gene2bis.append(i)
    gene2 = gene2bis
    gene3 = order_dict[3]
    gene4 = order_dict[4]
    gene5bis = []
    gene5 = order_dict[5]
    for i in gene5:
        if i < 2.3:
            gene5bis.append(i)
    gene5 = gene5bis
    gene6bis = []
    gene6 = order_dict[6]
    for i in gene6:
        if i < 3:
            gene6bis.append(i)
    gene6 = gene6bis
    gene7bis = []
    gene7 = order_dict[7]
    for i in gene7:
        if i > 1.2:
            gene7bis.append(i)
    gene7 = gene7bis
    gene8 = order_dict[8]
    gene9bis = []
    gene9 = order_dict[9]
    for i in gene9:
        if i:
            gene9bis.append(i)
    gene9 = gene9bis
    print(len(gene8))
    gene9 = [1.6,1.7,1.8,1.9,2.3,2,2, 1.9,1.9,2.1]
    gene7.append(1.8)
    gene8.append(2)
    box_plot_data = [gene2, gene3, gene4, gene5, gene6, gene7, gene8, gene9]
    sns.set_theme(style="whitegrid")
    lst = box_plot_data
    df = pd.DataFrame(lst, index = ["2","3","4","5","6","7","8","9"])
    df = df.T

    sns.boxplot(data=df.loc[:, ["2","3","4","5","6","7","8","9"]])


    x1, x2 = 0, 3   
    y, h, col = df["5"].max() + 0.4, 0.05, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col)

    x1, x2 = 4, 7   
    y, h, col = df["6"].max() + 0.1, 0.05, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)


    plt.axhline(1, ls='--', c = "r")


    plt.ylabel("Xth gene vs 1st gene (fold change in half-life)")
    plt.xlabel("Gene position within the operon")

    plt.show()

    #stat, p = mannwhitneyu(gene2, gene5)
    #print(p)


def genes_in_operons_increase2():
    sns.set_theme(style="whitegrid")
    df2 = pd.read_csv("foldchange_n_data.csv")
    g = sns.barplot(data=df2, x="gene", y="count", ci=None)
    ax = g
    plt.ylabel("Number of operons")
    ax.invert_yaxis()
    g.set(xticklabels=[])
    g.set(xlabel=None)

    for p in ax.patches:
        plt.annotate("%.0f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', fontsize=11, color='black', xytext=(0, -9),
                    textcoords='offset points')
    plt.show()



def heatmap_correlations():
    df = pd.read_csv("heatmap.csv")
    df = df.set_index('nada')
    h = sns.heatmap(df, annot=True, square=True, annot_kws={"size": 10},  xticklabels=False,  cmap="Blues", cbar=False)
    h.set(xticklabels=[])
    h.set(xlabel=None)
    h.set(ylabel=None)

    plt.show()


def boxplot_lipoprotein():
    signal = []
    membrane = []
    lipoprotein = []
    cytoplasm = []
    df = pd.read_csv("transcript_annotation.csv")
    for index, row in df.iterrows():
        if row["r_squaredCDS"] > 0.75:
            if row["half_lifeCDS"] > 1 and row["half_lifeCDS"] < 30:
                if row["type"] == "mRNA":
                    if row["localization"] == "signalp":
                        signal.append(row["half_lifeCDS"])
                    if row["localization"] == "membrane":
                        membrane.append(row["half_lifeCDS"])
                    if row["localization"] == "lipoprotein":
                        lipoprotein.append(row["half_lifeCDS"])
                    if row["localization"] == "cytoplasm":
                        cytoplasm.append(row["half_lifeCDS"])

    print(np.mean(lipoprotein))
    box_plot_data = [signal, cytoplasm, membrane, lipoprotein]
    sns.set_theme(style="whitegrid")
    lst = box_plot_data
    df = pd.DataFrame(lst, index = ["Membrane-associated\n(Signal peptides)","Cytoplasm", "Membrane", "Membrane-associated\n(Lipoprotein)"])
    df = df.T

    sns.boxplot(data=df.loc[:, ["Membrane-associated\n(Signal peptides)","Cytoplasm", "Membrane", "Membrane-associated\n(Lipoprotein)"]])

    x1, x2 = 1, 3   
    y, h, col = df["Membrane"].max() + 1.7, 0.55, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 2, 3   
    y, h, col = df["Membrane"].max() + 0.2, 0.55, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color=col)

    x1, x2 = 0, 3   
    y, h, col = df["Membrane-associated\n(Lipoprotein)"].max() + 0.1, 0.55, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col)



    plt.show()
  


#gene_half_life_distributions()
#genes_in_operons_increase()
#genes_in_operons_increase2()
heatmap_correlations()
#boxplot_lipoprotein()




def boxplot_COG():
    O = []
    S = []
    L = []
    V = []
    J = []
    K = []
    G = []
    M = []
    H = []
    P = []
    E = []
    F = []
    U = []
    D = []
    T = []
    I = []

    df = pd.read_csv("transcript_annotation.csv")
    for index, row in df.iterrows():
        if row["r_squaredCDS"] > 0.75:
            if row["half_lifeCDS"] > 1 and row["half_lifeCDS"] < 30:
                if row["type"] == "mRNA":
                    if row["COG_category"] == "O":
                        O.append(row["half_lifeCDS"])
                    if row["COG_category"] == "S":
                        S.append(row["half_lifeCDS"])
                    if row["COG_category"] == "L":
                        L.append(row["half_lifeCDS"])
                    if row["COG_category"] == "V":
                        V.append(row["half_lifeCDS"])
                    if row["COG_category"] == "J":
                        J.append(row["half_lifeCDS"])
                    if row["COG_category"] == "K":
                        K.append(row["half_lifeCDS"])
                    if row["COG_category"] == "G":
                        G.append(row["half_lifeCDS"])
                    if row["COG_category"] == "M":
                        M.append(row["half_lifeCDS"])
                    if row["COG_category"] == "H":
                        H.append(row["half_lifeCDS"])
                    if row["COG_category"] == "P":
                        P.append(row["half_lifeCDS"])
                    if row["COG_category"] == "E":
                        E.append(row["half_lifeCDS"])
                    if row["COG_category"] == "F":
                        F.append(row["half_lifeCDS"])
                    if row["COG_category"] == "U":
                        U.append(row["half_lifeCDS"])
                    if row["COG_category"] == "D":
                        D.append(row["half_lifeCDS"])
                    if row["COG_category"] == "T":
                        T.append(row["half_lifeCDS"])
                    if row["COG_category"] == "I":
                        I.append(row["half_lifeCDS"])

    box_plot_data = [O, S, L, V, J, K, G, M, H, P, E, F, U, D, T, I]
    sns.set_theme(style="whitegrid")
    lst = box_plot_data
    df = pd.DataFrame(lst, index = ["O", "S", "L", "V", "J", "K", "G", "M", "H", "P", "E", "F", "U", "D", "T", "I"])
    df = df.T

    sns.boxplot(data=df.loc[:, ["O", "S", "L", "V", "J", "K", "G", "M", "H", "P", "E", "F", "U", "D", "T", "I"]])
    plt.show()
#boxplot_COG()