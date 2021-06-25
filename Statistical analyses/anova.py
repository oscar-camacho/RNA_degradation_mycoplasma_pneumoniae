import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import operator
import scipy
from scipy.stats import f_oneway
from statsmodels.formula.api import ols
import statsmodels.api as sm
import statistics
import statsmodels.stats.multicomp as mc
from bioinfokit.analys import stat


df = pd.read_csv("transcript_annotation.csv",  encoding = "unicode_escape")


for index, row in df.iterrows():
    if row["ORF"][0:2] == "nc":
        df = df.drop([index])

for index, row in df.iterrows():
    if row["ORF"][0:4] == "MPNt":
        df = df.drop([index])

for index, row in df.iterrows():
    if row["r_squared CDS"] < 0.7:
        df = df.drop([index])

df = df[df['COG_category'].notna()]   ######



model = ols('half_lifeCDS ~ C(COG_category)', data=df).fit()  #######
aov_table = sm.stats.anova_lm(model, typ=2)
#print(aov_table)

def anova_table(aov):
    aov['mean_sq'] = aov[:]['sum_sq']/aov[:]['df']

    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])

    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*aov['mean_sq'][-1]))/(sum(aov['sum_sq'])+aov['mean_sq'][-1])

    cols = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[cols]
    print(aov)

anova_table(aov_table)

print("\n")
mc = mc.MultiComparison(df['half_lifeCDS'],df['COG_category'])    #################
mc_results = mc.tukeyhsd()
print(mc_results)

group_dict = dict()
for index, row in df.iterrows():
    list = row["COG_category"].split(", ")
    for i in list:
        if i not in group_dict:
            group_dict[i] = []    #################3
        group_dict[i].append(row["half_lifeCDS"])


print("\n")
for k,v in group_dict.items():
    print((str(k)) + ". N: " + str(len(v)) + ". Mean: " + str(statistics.mean(v)) + ". Std: " + str(np.std(v)))
