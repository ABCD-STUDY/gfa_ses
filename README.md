# gfa_ses
Group factor analysis for socioeconomic measures, first released in ABCD 3.0 (derived from ABCD 2.0.1 data).

Author: Marybel Robledo Gonzalez
Date: April 30, 2020

Overview
These scores will be included in the speciality_sum_score instrument with the 3.0 release.  

Variable Name	& Description
latent_factor_ss_general_ses:	General latent factor of economic, social, and physiological well-being.
latent_factor_ss_social:	Latent factor for youth perceived social support.
latent_factor_ss_perinatal:	Latent factor for perinatal health. 

These three latent factors were derived with baseline data with the ABCD 2.0.1 release for an N = 8,158 subjects who had complete available data for all of the measures of interest.  Details are described in Gonzalez et. al. (2019 preprint).

We selected 22 variables thought to be proximal measures of the broader socioeconomic context experienced by children.  We grouped the 22 variables into the categories of economic, social, community, adverse childhood experiences, physiological and perinatal well-being.  The 22 variables were subjected to a Group Factor Analysis using the R package ‘GFA’ (Klami, et al., 2015).  Code was adapted from publicly available scripts for a GFA analyses with ABCD data (Paulus, et al., 2019).

Three separate and orthogonal latent factors were extracted, overall accounting for ~26% of the variance among the 22 variables.  The variables are described in detail in the supplementary material in Gonzalez et. al. (2019 preprint).  To test the stability and robustness of the latent factors, we completed 10 different iterations of the GFA.  Robust latent factors were chosen based on latent factor loadings that met a 0.9 correlation threshold across all 10 iterations.  Robust factor loadings across all 10 GFA iterations were averaged.  Separate robust GFAs were examined in split-half samples to test replication of the latent factor loadings.  Robust GFA latent factors accounting for more than 5% of the GFA variance were chosen.   

latent_factor_ss_general_ses was positively associated with socioeconomic status (SES) as measured by income-to-needs.  There was no significant association between SES and the other latent factors.  

All three latent factors were positively associated with total cortical surface area and NIH total cognition scores.   

Key input files:

The R markdown files reads in an “.Rds” file (not provided) for the ABCD 2.0.1 NDA data containing the 22 variables of interest, demographic variables, total cortical surface area, and NIH total cognition scores. 

R Markdown Files and Scripts:

“ABCD_SocioEconomicContext-LatentFactors-GFA-variables.Rmd” -  this file reads in the ABCD 2.0.1 Rds file and selects the 22 variables of interest and recoding of variables.  The file also creates the matrices and labels necessary to run the GFA using the “ABCD_SocioEconomicContext_GFA_robust.Rmd” file.  

“robustGFAFcn.R” – This script must first be run to define the necessary functions for the robust GFA.  This script was originally implemented in Paulus, et al. (2019). 

“ABCD_SocioEconomicContext_GFA_robust.Rmd – this file is used to implement the group factor analysis after all of the necessary variables have been defined/recoded and grouped into matrices using the “ABCD_SocioEconomicContext-LatentFactors-GFA-variables.Rmd” file.  The code implements a robust GFA, producing the median loadings across 10 iterations of the GFA for latent factors that met a 0.9 correlation threshold between each iteration.  Code for producing tables and visualizations for the latent factor loadings is also provided. 

“ABCD_SocioEconomicContext-LatentFactors_GAMM_Figs_Tables.Rmd” - this file contains the code to merge the values for the GFA latent factors with demographic and outcome variables by subject into a data frame.  Code is provided for running the generalized additive mixed-effect models (GAMMs) used to associate the latent factors with income-to-needs, total cortical surface area and total cognition scores, as well as the code for producing the figures and tables in Gonzalez, et al. (2019 preprint). 




References

1.	Gonzalez, M. R., Palmer, C. E., Uban, K. A., Jernigan, T. L., Thompson, W. K., & Sowell, E. R. (2019). Economic, social, and physiological resilience predict brain structure and cognitive performance in 9-10-year-old children. bioRxiv, 852988.
https://www.biorxiv.org/content/10.1101/852988v1

2.	Klami, A., Virtanen, S., Leppäaho, E., & Kaski, S. (2015). Group factor analysis. IEEE transactions on neural networks and learning systems, 26(9), 2136-2147.

3.	Paulus, M. P., Squeglia, L. M., Bagot, K., Jacobus, J., Kuplicki, R., Breslin, F. J., ... & Tapert, S. F. (2019). Screen media activity and brain structure in youth: Evidence for diverse structural correlation networks from the ABCD study. Neuroimage, 185, 140-153.


