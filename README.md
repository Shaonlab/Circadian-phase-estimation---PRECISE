PRECISE Documentation Suite
 
README.md
PRECISE: PREdicting CIrcadian phase from Stochastic gene Expression
PRECISE is an R-based supervised machine learning algorithm that is designed to accurately infer the circadian phase from RNA measurements of a test sample. The associated paper can be found here:
The RNA measurements could be from any experimental method, for example bulk/single-cell RNA-seq, qPCR, etc.  PRECISE builds on Gaussian Processes (GPs) to learn noisy and non-stationary oscillatory trajectories. 
This repository provides the R scripts and functions required for running PRECISE on different datasets. The R script has a modular design where each step (data input, normalization, bootstrapping, likelihood estimation and testing) are available as separate modules in the R script. 
 
Installation
If not already installed, R can be downloaded from the CRAN homepage. RStudio is available for download in this link. A detailed documentation to installation of both R and RStudio can be found at this webpage. 
PRECISE requires manual installation of ODeGP from a local folder supplied to the user.
1. Install ODeGP locally
Place the ODeGP source folder ( ODeGP-main/) anywhere on your computer. This folder contains the oscSetup.R script that will used to install local version of ODeGP package. The following command installs ODeGP (also included in the PRECISE_validation.R script (line 11) :
source("~/path_to_folder_where_ODeGP-main_present/ODeGP-main/oscSetup.R")
Then load the library:
library(ODeGP)
2. Install remaining dependencies
install.packages(c(
  "DEoptim", "tidyverse", "ggpubr", "ggpmisc", "circular"
))
3. Confirm installation
library(ODeGP)
sessionInfo()
You should see ODeGP listed among the loaded packages.
4. The PRECISE scripts
This repository contains two scripts, PRECISE_validation.R and PRECISE_prediction.R. The first can be used to reproduce Figures 3 and 4 in the main text while the second can be used for predicting phase of unknown samples.  
 

Instructions for running PRECISE validation

Input data structure required for running PRECISE_validation.R
In order to validate the performance of PRECISE and reproduce Figure 3 and 4 in the main paper, the following inputs are required: 
•	Two folders: one containing the training datasets and one containing the test datasets. These training and test datasets comprise the gene-expression data. Each folder should contain as many CSV files as  timepoints, named such that the timepoints are in increasing order (see Example below for details). This ordering is essential because the script uses index-based referencing. Each CSV file must contain: 
Column names: user-defined gene names
Entries of the columns: RNA counts for each gene in each cell/sample. 
Therefore a row represents the vector of RNA counts for different user-defined genes measured in a particular single cell/sample. 
•	Two CSV files: one containing training phases and one containing test phases corresponding to the time-point labels. The training/test phase files contain two columns:
Time: ZT or time post-Dexamethasone synchronization
Phase: the true circadian phase at that timepoint, estimated using a chosen gene (the same gene must be used for both training and test datasets). The true phase needs to be estimated using wavelet transforms; for details see Supplementary Info of the associated paper. 
Additional (Optional) Columns in the Example Datasets
In the datasets included in this repository, two extra columns are present:
•	Cell_ID_unique
Represents the field of view (FOV) and the cell number (e.g., "001_10" = FOV 001, cell 10)
•	Area_Cell
The area of the cell in pixels (corresponding to the segmentation mask)
These columns are not required to run PRECISE_validation.R. They are included because the manuscript analyzes the effect of cell area; they are provided for completeness but can be ignored for running PRECISE.
 
Example
The various datasets required for implementation of PRECISE_validation.R are provided in the docs folder.
It includes the following files:
1.	training_Data.xlsx : An Excel file containing the NIH3T3 training dataset used in the original study, with each sheet corresponding to a single timepoint. Each sheet includes five columns: Cell_ID_unique; the gene names (as column headers), where each entry represents the RNA count for that specific gene; and Area_Cell, which denotes the cell area in pixels corresponding to the full cytoplasmic mask. Since in this manuscript 4 genes (bmal1, nr1d1, nr1d2, tef) were measured hence the column names are  - “Cell_ID_unique”, “bmal1”, “nr1d1”, “nr1d2”, “tef” and  “Area_Cell”.
2.	test_Data.xlsx : an excel file for the NIH3T3 test data generated in the paper, where each sheet represents one time point. Same columns present as the training data.
3.	training_phase_nr1d1.csv : a csv file containing the circadian phases of the different time point cell populations determined using the Nr1d1 gene using wavelet analysis.
4.	test_phase_nr1d1.csv : a csv file containing the true circadian phases of the test samples. Determined using the same gene as the one used for training data i.e. Nr1d1.
Before running the script : 
1.	Create two folders:
•	training_data/
•	test_data/
2.	Convert each time point sheet from the Excel files into an individual .csv file and place them in the appropriate folder.
3.	Name the files according to their timepoints (for example, 08_final.csv, 16_final.csv, and so on until the last timepoint). Using this naming convention ensures that the files are read in the correct chronological order, which is important because the code relies on index-based access at several steps.
After setting up the folders, update the relevant file paths in the code and run the PRECISE_validation.R script. Upon completion, the script will generate a set of graphs — one for each test timepoint analyzed — showing how the PRECISE-inferred phase deviates from the true phase (estimated using wavelet analysis) as a function of the number of cells averaged per sample, starting from single-cell measurements.
 


