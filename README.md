# Traits-Genetic-Correlations
These scripts were used to explore genetic correlations and residual correlations for 
There are three files:
1) The Linux script to run the multivariate mixed linear model implemented in GIBBS3F90 (gibbs_runopt1.sh) 
2) There is also a sample of the parameter file used to run GIBBS3F90 (renf90.par), which was generated after running renumf90 using the parameter file renum.txt
3) the R script for the phenotypic, genetic and residual correlations (Corr_script.R). 
The R script used the postgibbs_samples file obtained after running blupf90. It contains the code for creating correlation tables with significant levels, PCA biplots, and undirected graphical models. It also include code to put plots for phenotypic, genetic and residuals correlations in singles plots
