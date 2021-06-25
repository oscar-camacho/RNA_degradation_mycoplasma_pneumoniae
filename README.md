# Genome-wide analysis of RNA degradation in Mycoplasma pneumonaie

This is a GitHub repository with Supplementary Materials and Python scripts from the master's thesis entitled "Genome-wide analysis of RNA degradation in Mycoplasma pneumoniae". 

- Author: Oscar Camacho Martínez
- Supervisor: Professor Luis Serrano Pubul
- Centre for Genomic Regulation
- Universitat Pompeu Fabra


## Supplementary Materials

Supplementary Table 4 is included. It consists of a table with half-life, R² of the adjustment of RNA decay to the exponential decay curve, fold change in expression due to glucose starvation and novobiocin treatment (10 µg/ml), and transcription rate for each annotated gene in M. pneumoniae.

## Python scripts and datasets

### Genome annotation datasets + features

- transcript_annotation.csv: Dataset with the annotated coordinates for each gene (ORF, ncRNA, tRNA, rRNA), with all the calculated features for each one used in association analyses with half-life.
- windows_annotation.csv: Dataset with the annotated coordinates for each sliding window in which the genome was divided, with all the calculated features for each one used in association analyses with half-life.

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
