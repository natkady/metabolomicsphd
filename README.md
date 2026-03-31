# metabolomicsphd
Supplementary material for my metabolomics PhD thesis

File id_vs_id_assess_criteria.csv has the data to calculate the means of R2X (PCA), R2X (PLS-DA), R2Y, and Q2Y for the group comparisons infant sex, foetal versus maternal (fvm) and inner versus outer (ivo). This data is essential for the calculation of the values in chapter 3, section 3.2.2.2, table 3.3. This will be used as a supplementary file for the thesis.

#Thesis R Scripts

#YPOPS



#UPBEAT

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
#GEO

The GEO shiny app contains 2 R files: geo_analysis.R and diff_cell.R. The main code for the app is within geo_analysis.R, with diff_cell.R containing the functions needed to produce analyses and graphs within geo_analysis.R. For geo_analysis.R to run, JAVA2 must be installed on the local device, and diff_cell.R must be read first.


