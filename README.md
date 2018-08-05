################################################################################
## File:             ReadMe.txt                                               ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.07.2018                                               ##
##                                                                            ##
## Contains description of the reproducing scripts for the manuscript         ##
## 'Nonparametric imputation by data depth'                                   ##
## by Pavlo Mozharovskyi, Julie Josse, and Francois Husson.                   ##
##                                                                            ##
################################################################################

Current folder contains the sources of R-package 'imputeDepth' and files 
required to reproduce all the experiments and data-connected figures of the 
manuscript 'Nonparametric imputation by data depth' (2017) 
by Pavlo Mozharovskyi, Julie Josse, and Francois Husson. Detail descriptions 
for utilisation of the files follow.

Folder 'imputeDepth' as well as the archive 'imputeDepth_0.1.0.tar.gz' contain 
sources of the R-package 'imputeDepth', whose functions implement the 
imputation methods from the manuscript, and which are used in the experiments.

Folder 'data' contains files of the data sets 'banknotes.dat', 
'bloodtransfusion_gp.dat', 'cows.dat, 'glass.dat'; see the manuscript and 
references therein for details.

File 'impute.functions.R' contains routines for competing imputation methods 
implemented by the authors and a function for introducing of missing data into 
a complete data set.

Files 'exp.impute.*.R' contain R-sources of the experiments, namely:

 - files 'exp.impute.StudentT_MCAR.R', 'exp.impute.StudentToutl_MCAR.R', 
   'exp.impute.StudentT_MAR.R', 'exp.impute.lowrank.R', 'exp.impute.large.R', 
   'exp.impute.skewed.R', 'exp.impute.moon.R' contain scripts for simulation 
   experiments with missing completely at random (MCAR), MCAR with outliers, 
   missing at random (MAR), lowrank models, contamination in higher 
   dimension, skewed normal distribution and distribution wiht non-convex 
   support, respectively, see Section 4.2.
   
 - files 'exp.impute.banknotes.R', 'exp.impute.bloodtransfusion.R', 
   'exp.impute.cows.R', 'exp.impute.glass.R' contain scripts for real-data 
   experiments; see Section 5.3. The data sets are supposed to be in the folder 
   'data', if not the path in the R-file should be changed accordingly.
   
 - file 'exp.impute.quantile.R' contains the R-script for reproducing the 
   quantile experiment for multiple imputation; see Section 5.3.
   
 - file 'exp.impute.multiple.R' contains the R-script for reproducing the 
   complete multiple imputation experiment for a regression model, 
   see Section 5.2.

Files 'plot.intro.R', 'plot.impute.skewed.R', 'plot.impute.moon.R', 
'plot3d.optim.R', 'plot.one.density.R' , 'plot.impute.Tukey.R' reproduce 
Figures 1 and 2 (Section 1), Figures 4 and 5 (Section 4.2.6), 
Figures 1--3 (supplement), respectively.

File 'plot3d.realdata.R' produces pairs() and pairs3d() plots used to get 
intuition by visualising the real data.
