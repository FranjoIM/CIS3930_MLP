# Overview
In this study, we leverage multi-modal, longitudinal data from ABCD
study and multi-kernel machine learning approaches to accurately classify ambiguous
OCD cases and increase power to downstream genetic analyses. We leverage logistic
regression model with elastic net penalty to train each kernel independently, followed by
accuracy-weighted case probability score summed across kernels 

# DataSets
ABCD Study data were accessed through National Institute of Mental Health (NIMH) 60
Data Archives (NDA). We are using different types of Datasets for building our prediction model.
* Genetic Data : Genetic data for participants were derived from blood and saliva samples and 66
genotyping on National Institute for Drugs and Addiction (NIDA)
* NeuroImaging Data: Tabulated summaries of post-processed neuroimaging data obtained via structural 84
magnetic resonance imaging (MRI) and diffusion MRI scans were downloaded from 85
NDA for imaging kernels.
* Questionarrie and phenotype Data: Tabulated questionnaire and phenotype data were downloaded from NDA for 93
questionnaire kernel and label definition.

A total of 6 processed data sets were defined for all samples after processing the above datasets:
```yaml
1) Congnitive and Behavioural Phenotype
2) Genotypes
3) Diffusion MRI
4) Structural MRI
5) Rest State FMRI
6) Task FMRI
```

# Machine Learning Model:
## Method
We built 6 different kernels based on these datasets categories and are combined using weighted probabilistic accuracy 
kernel to make prediction on final test data-set.
We trained each kernel using Logistic Regression model with elastic net as the penalty for every batch.
## Support Method
We compared model performance using Cohen’s kappa score to Cohen’s kappa scores 
reported in psychiatry and psychiatric epidemiology literature.
## Success Measurement
The aim of this project is to develop a machine learning model to classify individuals 
with OCD status. In order to measure our success, we have utilized confusion-matrix 
based assessments like accuracy, specificity and sensitivity. At the baseline, we 
expect accuracy of 0.5, sensitivity of 0.5, and specificity of 0.94, assuming prevalence of 
0.5.


# Languages and Tools used
* R-language for data preprocessing and cleaning the data.
* Machine learning Model: Logistic Regression, Elastic Net
* Python and Scikit learn used for computation and development of model.

# Results
Performance of each kernel in terms of accuracy has been visually represented. Phenotype kernel vastly outperformed all other kernels in terms of accuracy. Additionally, plots of receiver operative curves and precision-recall curves are shown for each individual kernel. Based on AUC metrics and F1 scores, we can see that phenotype kernel outperforms all other kernels. The kernels based on imaging and genetic data performed similarly, very close to expected baseline performance.


# Contributors

- Franjo Ivankovic
- Aryan Singh
- Shashank Kumar
- Richa Gupta

 
