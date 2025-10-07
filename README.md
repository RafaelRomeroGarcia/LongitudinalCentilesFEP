# Neurobiology and Centiles in Psychosis

This repository contains code and data created in support of the project *"Atypical brain maturation in psychosis is linked to long-term cognitive and symptom progression"*. All code was written in R. Folders, files, and first steps are described below.

## **Data**

The `Data` folder contains the data required for running the analyses. Here are the files that need to be downloaded and stored in a specific location. The remaining files will be automatically generated:

-	The code used to compute centiles (`BrainCharts` folder) is available at https://github.com/brainchart/Lifespan

-	The `all_microsc_DesikanKilliany68.csv` file is available at https://github.com/netneurolab/netneurotools

-	The code used for the PCA-CCA analyses (`cca_pls_toolkit_dev-grcca` folder) is available at https://github.com/anaston/cca_pls_toolkit

## **First steps**

1.	Download `Code` folder, which contains the scripts and functions used for the analyses.

2.	Download `Data` folder, which contains data used to run the analyses.



## **Code**

The `Code` folder contains all the code required for running the analyses and generate the figures. Don't forget to update the paths. 




-	[CognitionStudy.R](CODE/CognitionStudy.R) –This script performs a mixed linear regression looking for the relationship between the overall cognitive functioning assessment of the cohort and covariates such as age, sex , medication and regional centiles.
-	[MesulmanTurkeyStudy.R](CODE/MesulmanTurkeyStudy.R) – This script performs statistical analysis and visualization on neuroimaging data.It includes data loading, data transformation, ANOVA, and Tukey's HSD test.
         The script uses generic file paths and column names for ease of adaptation to other datasets.
-	[CaseBaseStructuralStudy.R](CODE/CaseBaseStructuralStudy.R)-This script performs a linear regression looking for the relationship between the baseline centiles of the cohort and covariates such as age, sex or medication.
-	[LongitudinalStructuralStudy.R](CODE/LongitudinalStructuralStudy.R)-This script performs a mixed linear regression looking for the relationship between the centiles of the cohort and covariates such as age, sex or medication.
-	[SymptomsLongitudinalStudy.R](CODE/SymptomsLongitudinalStudy.R)-This script performs a mixed linear regression looking for the relationship between the Symptoms of the cohort and covariates such as age, sex , medication and regional centiles.
-	[mediationStudy.R](CODE/mediationStudy.R)-This script performs a mediation analyses using bootstrapping applied to Mixed Models 




## **License**

This project is licensed under the terms of the [GNU General Public License v3.0 license](LICENSE).

