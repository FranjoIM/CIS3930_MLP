```bash
# Load UFRC module and begin a session
module load ufrc
srundev \
  --account=carolmathews \
  --qos=carolmathews \
  --time=06:00:00 \
  --mem=15G \
  --cpus-per-task=2 \
  --ntasks=1

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Assign paths to variable names
WD="..."
```
```bash
# Create holding directory
mkdir -p CSI
cd ${WD}/CSI

# List files to stitch
ls ${WD}/0504/ABCD.chr*.imputed.filtered.sorted.bcf.gz > \
    ${WD}/CSI/Chunks.txt

# Load in BCFtools 
module load bcftools/1.13

# Concatonate files 
bcftools concat \
    --file-list ${WD}/CSI/Chunks.txt \
    --output-type z \
    --threads 2 \
    --output ${WD}/CSI/ABCD.vcf.gz

# Filter variants to keep
awk 'NR == 1 || $9 < 0.001' ${WD}/1002/ocd_aug2017 | cut -f2 | awk 'NR>1' >  ${WD}/CSI/OCD_SigHit

bcftools view \
    --output ${WD}/CSI/ABCDFilt.vcf.gz \
    --output-type z \
    --threads 2 \
    --include ID==@OCD_SigHit \
    ${WD}/CSI/ABCD.vcf.gz

bcftools view \
    --output ${WD}/CSI/ABCDFinal.vcf \
    --output-type v \
    --threads 2 \
    --include "INFO>0.9 & AF>=0.1 & AF<=0.9" \
    ${WD}/CSI/ABCDFilt.vcf.gz

bcftools query \
    --format '%ID\t[\t%DS]\n' \
    --print-header \
    ${WD}/CSI/ABCDFinal.vcf > \
    ${WD}/CSI/DosageMatrix

bcftools query \
    -list-samples \
    ${WD}/CSI/ABCDFinal.vcf > \
    ${WD}/CSI/SampleIDs

# Finish preparation in R
module load R/4.1
R
```
```r
########################
### DATA PREPARATION ###
########################

# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import data
OCD <- read.table("./OCD.pheno", sep="\t", header=F)
GEN <- read.table("./CSI/DosageMatrix", sep="\t", header=F)
SAM <- read.table("./CSI/SampleIDs", sep="\t", header=F)

# Remove empty column
GEN <- GEN %>%
    select(-V2)

# Name Samples 
colnames(GEN) <- c("SampleID", SAM$V1)

# Transpose
GENt <- as_tibble(cbind(SampleID = names(GEN), t(GEN)))
colnames(GENt) <- GENt[1,]
GENt <- GENt[-1, ] 

# Keep only samples with information on both
KEEP <- intersect(GENt$SampleID, OCD$V1)

OCD <- OCD %>%
    filter(V1 %in% KEEP)

GENt <- GENt %>%
    filter(SampleID %in% KEEP)

# Change class to numeric where applicable and save
GENt <- GENt %>%
    mutate_at(vars(rs1806446:rs111926263), as.numeric)
saveRDS(GENt, file="./CSI/Genetic.rds")

# Update OCD colnames and save phenotype data
colnames(OCD) <- c("SampleID", "OCD")
saveRDS(OCD, file="./CSI/OCD.rds")

# Create Batches
BAT <- OCD %>%
    mutate(Experimental = ifelse(OCD==1, 1, 0))

BAT_Case <- OCD %>%
    filter(OCD==2) %>%
    pull(SampleID)

BAT_Cont <- OCD %>%
    filter(OCD==0) %>%
    pull(SampleID)

set.seed(1)
SEL <- c(sample(BAT_Case, 50), sample(BAT_Cont, 50))
BAT <- BAT %>%
    mutate(Test = ifelse(SampleID %in% SEL, 1, 0))
BAT_Case <- BAT_Case[!BAT_Case %in% SEL]
BAT_Cont <- BAT_Cont[!BAT_Cont %in% SEL]

set.seed(1)
SEL <- c(sample(BAT_Case, 50), sample(BAT_Cont, 50))
BAT <- BAT %>%
    mutate(Valid = ifelse(SampleID %in% SEL, 1, 0))
BAT_Case <- BAT_Case[!BAT_Case %in% SEL]
BAT_Cont <- BAT_Cont[!BAT_Cont %in% SEL]

i=1
set.seed(1)
while(length(BAT_Cont)>=94){
    name <- paste0("Train_",i,collapse="")
    SEL <- c(sample(BAT_Case, 100), sample(BAT_Cont, 94))
    BAT <- BAT %>%
        mutate(!!sym(name) := ifelse(SampleID %in% SEL, 1, 0))
    BAT_Cont <- BAT_Cont[!BAT_Cont %in% SEL]
    i=i+1
}

BAT <- BAT %>%
    mutate_at(vars(Experimental:Train_77), as.numeric) %>%
    rowwise() %>% 
    mutate(ALL = sum(c_across(all_of(names(BAT)[3:82]))),
        TRAIN = sum(c_across(all_of(names(BAT)[6:82]))))

saveRDS(BAT, file="./CSI/Batch.rds")

###################
### SIMULATIONS ###
###################

CODE <- select(GENt, SampleID)

# Genetic Data Simulation
SNP <- select(GEN, SampleID)
colnames(SNP) <- "RSID"
CODE$Dummy <- paste0(stri_rand_strings(8701, 4, '[A-Z]'), 
    stri_rand_strings(8701, 6,'[0-9]'))
SNP$Dummy <- paste0("rs", stri_rand_strings(5654, 10,'[0-9]'))
GEN_Dummy <- CODE %>%
    left_join(GENt, by=c("SampleID"="SampleID")) %>%
    select(-SampleID) %>%
    column_to_rownames(var="Dummy")
GEN_Dummy <- data.matrix(GEN_Dummy)
GEN_Dummy <- matrix(sample(GEN_Dummy),nrow=nrow(GEN_Dummy))
rownames(GEN_Dummy) <- CODE$Dummy
colnames(GEN_Dummy) <- SNP$Dummy
GEN_Dummy <- as.data.frame(GEN_Dummy) %>%
    rownames_to_column(var="SampleID")
saveRDS(GEN_Dummy, file="./CSI/Genetic_Dummy.rds")

# Phenotype Unlabeling
OCD_Dummy <- CODE %>%
    left_join(OCD, by=c("SampleID"="SampleID")) %>%
    select(-SampleID) %>%
    rename(SampleID=Dummy)
saveRDS(OCD_Dummy, file="./CSI/OCD_Dummy.rds")

# Batch Unlabeling
BAT_Dummy <- CODE %>%
    left_join(BAT, by=c("SampleID"="SampleID")) %>%
    select(-SampleID) %>%
    rename(SampleID=Dummy)
saveRDS(BAT_Dummy, file="./CSI/Batch_Dummy.rds")
```
