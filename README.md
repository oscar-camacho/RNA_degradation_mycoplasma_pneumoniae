# Genome-wide analysis of RNA degradation in Mycoplasma pneumonaie

Oscar Camacho<sup>1,2,\*</sup>, Samuel Miravet-Verde<sup>1</sup>, Luis Serrano<sup>1,2,3,\*</sup> and Marc Weber<sup>1</sup>

This is a GitHub repository with Supplementary Data and Python scripts for the article entitled "Genome-wide analysis of RNA degradation in Mycoplasma pneumoniae". 



## Supplementary Tables

- supplementary_table_S5-S7.csv:
  -  Table S5: Values and statistics determined for each gene in MPN. Such values include RNA half-life and coefficient of determination (R squared) from the fit of the expression values to a negative exponential decay for WT and RNase III KO strains, fold change in half-life in RNase III KO, normalized gene expression for WT and RNAse III KO strains, and log2 fold change in gene expression in RNase III KO, treatment with novobiocin at 10 ug/ml and overexpression RNases R, J1, J2 and Y.





- Half-life, coefficient of determination for the RNA decays and fold change in expression due novobiocin treatment (10 µg/ml) for each annotated gene in M. pneumoniae genome.
- supplementary_table_S6.csv: Fold changes in gene expression due to RNase III knockout and RNase R, J1, J2 and Y overexpression for each gene in M. pneumoniae genome.
- supplementary_table_S7.csv: Annotated coordinates and all the calculated features for each gene (ORF, ncRNA, tRNA, rRNA) in M. pneumoniae genome.
- supplementary_table_S8.csv: Annotated coordinates and all the calculated features for each sliding window in which the M. pneumoniae genome was divided.

## Python scripts

### Half-life calculations

- half-life_genes.py: Python script used to normalize the RNA-Seq data from the different time-points and calculate the RNA decay constant and half-life for each gene in the genome.
- half-life_windows.py: Python script used to normalize the RNA-Seq data from the different time-points, divide the genome into sliding windows and calculate the RNA decay constant and half-life for each window in the genome.

### Random forest classifier 
- rf_classifier.ipynb: Jupyter notebook used to perform the classification of half-life of genes based on different gene/genomic features, using a random forest classifier algorithm.

### Clustering of profiles
- clustering.py: Python script used to classify the RNA decay profiles into different groups/clusters using a k-means algorithm.

### Statistical analyses

- anova.py: Python script used to calculate the statistical values in association analyses using 1 numerical + 1 categorical variable.
- correlation_genes.py: Python script used to calculate the statistical values in association analyses at gene-level, using 2 numerical variables.
- correlation_windows.py: Python script used to calculate the statistical values in association analyses at windows-level, using 2 numerical variables.

### Plots

- plots_analyses.py: Python script used to generate all the plots included in the study.
- plot_profile.ipynb: Jupyter notebook used to plot the half-life profiles, based on the half-life calculated for each window.
