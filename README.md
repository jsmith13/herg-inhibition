# hERG Inhibition
An inferential model for small molecule inhibition of the human Ether-à-go-go channel.
http://bicologie.com/jakeallensmith/projects/herg-inhibition/

## Objectives
- Build an inferential model to identify molecular properties contributing to inhibition of the human Ether-à-go-go (hERG) ion channel.
- Generate an interactive visualization of the modeling process.

## Dataset
- hERG assay data collected from the PubChem bioassay repository (https://pubchem.ncbi.nlm.nih.gov/bioassay/588834).
- Molecular descriptors calculated using RDKit and KNIME.

## Contents 
- app.R  
An R Shiny application producing an interactive visualization of the model fitting process.

- hERG_descriptors.csv  
The primary dataset.

- hERG_validation_descriptors.csv  
An independent test dataset.

- model_selection.Rmd  
Variable selection and model cross-validation.

- potency_model.Rmd  
A secondary analysis fitting a predictive model for hERG IC50. 

- validation.Rmd  
Validation of the model against the test set.

- writeup.Rmd  
A final report.
