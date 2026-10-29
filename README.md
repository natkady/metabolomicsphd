# metabolomicsphd
Supplementary material for my metabolomics PhD thesis

File id_vs_id_assess_criteria.csv has the data to calculate the means of R2X (PCA), R2X (PLS-DA), R2Y, and Q2Y for the group comparisons infant sex, foetal versus maternal (fvm) and inner versus outer (ivo). This data is essential for the calculation of the values in chapter 3, section 3.2.2.2, table 3.3. This will be used as a supplementary file for the thesis.

# Thesis R Scripts

## YPOPS

The files for YPOPS data analysis are as follows:

# PREPROCESSING

Metabolomics_Data_Pre_Processing_Procedure_And_PCA_Plots.R

Package Installer.R

C3_transpose_after_SVR.R

PCA_of_Data_Matrix_Before_Pre_Processing.R

# MULTIVARIATE ANALYSIS

ID_V_ID_Means.R

Inner_or_Outer_Expl_Data_Analysis.R

Foetal_or_Maternal_Expl_Data_Analysis.R

Inf_Sex_Expl_Data_Analysis.R

# MUMMICHOG ANALYSIS

Mummichog_Pre_Transform_Data.R

Mummichog_Prep_Sep_Comparisons.R

Mummichog_Prep.R

Correction_For_Multiple_Testing.R



## UPBEAT

In order of how the analysis was completed, here are the files used:

1_Adiposity_Binary_Variable.Rmd

2_NMR_Visits_And_obesity.Rmd

3_Metabol_Data_Adjust.Rmd

4_Metabol_Missingness.Rmd

5_GDM_Metabol_Add.Rmd

6_Data_MVI_GDM_ADIP.Rmd

7_Primary_Metabolity_Analysis.Rmd

  7A_Functions_For_Primary_Metabolite_Analysis.Rmd
  
  7B_Rate_of_Change_in_Metabolites_Btwn_Visits.Rmd
  
8_Significant_Metabolite_Comparison_from_GAMS.Rmd

  8B_Data_Simulation_PVAL.R
  
  8C_Graph_Adjustments.Rmd
  
  8C_Interaction_GAM_function.R
  
9_GAM_Cross_Validation.Rmd

Files 7_, 7A_, and 7B_ should be run concurrently i.e. 7A_ and 7B_ should be run with 7_. Similarly, 8_, 8B_, and both 8C_ should be run concurrently.

## GEO

The GEO shiny app contains 2 R files: 

geo_analysis.R 

diff_cell.R. 

The main code for the app is within geo_analysis.R, with diff_cell.R containing the functions needed to produce analyses and graphs within geo_analysis.R. 

For testing of the GEO accession data analysis app, please use the files: helper-setup.R, test-analysis-functions.R, test-app.R and test-validation-geo2r.R

To run the RShiny app, please download the files diff_cell.R and geo_analysis.R and save them in the same folder on your machine. To run the app, open the geo_analysis.R file in RStudio and click 'Run App'. Steps to running the code file can also been seen at the start of the geo_analysis.R file. The file diff_cell.R is required to be in the same folder as the geo_analysis.R file so that the app may run its functions.

To validate the app and its functions, please download helper-setup.R, test-analysis-functions.R, test-app.R and test-validation-geo2r.R and save them in the same folder as the diff_cell.R and geo_analysis.R files.

To be able to validate the app and functions, you must run the helper-setup.R file first. Within the file, please state the directory where the GEO files have been saved with the “” in app_dir <- i.e.

app_dir <-"[directory pathway where GEO files have been saved eg. C:/Users/natkady]"

replace [directory pathway where GEO files have been saved eg. C:/Users/natkady] with the pathway to the folder. Using the example this would look like:

app_dir <- “C:/Users/natkady”

and save the file before running. Once it has been saved, you may run the test- files.



