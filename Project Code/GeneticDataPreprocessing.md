# 00. Staging Environment

In this segment, we will create all directories and path variables (to be recycled in SLURM scripts). We will also upload data into respective directories.

```bash
# Load UFRC module and begin a session
module load ufrc
srundev  --account=carolmathews --qos=carolmathews --time=06:00:00

#  Create working, scripting, and log directories
mkdir -p .../WorkDir \
    .../SLURM_Scripts/LogFiles

# Assign paths to variable names
PD="..."
WD="..."
SLURM="..."
LOGS="..."
RF="..."

# Migrate into parent directory
cd ${PD}
```

## 00.01 Upload Files

The following files should be uploaded to the working directory

1. abcd_ksad01.txt (Parent KSADS-5 (v1))
2. abcd_ksad01_definitions.csv (Dictionary for parent KSADS-5 (v1))
3. ABCD_release3.0_.batch_info.txt
4. ABCD_release_3.0_QCed.fam
5. ABCD_release_3.0_QCed.bim
6. ABCD_release_3.0_QCed.bed

## 01.01 Phenotype Definitions

```bash
# Migrate into parent directory
cd ${PD}

# Load and launch R
module load R 
R
```

Run the following script to define OCD phenotypes. Recode values in the dataset for the easier undersanding:

1. 0 -> 0 = Negative (diagnosis not present)
2. 1 -> 1 = Positive (diagnosis presnt)
3. 555 -> N = Not administered in the assessment
4. 777 -> D = Declined to answer
5. 888 -> 0 = Question not asked due to answer to a previous question (equivalent to negative)
6. 999 -> K = Did not know
7. "NA" -> NA

```r
# Set Working directory to ${PD}
setwd("...")

# Load necessary packages
library(readr)
library(tidyverse)

# Load the dataset and associated dictionary
KSAD1P <- read_delim("./WorkDir/abcd_ksad01.txt", 
    delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE)[-1,]
KSAD2P <- read_delim("./WorkDir/ksads2daic_use_only_p01.txt", 
    delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE)[-1,]
DKSAD1P <- read_delim("./WorkDir/abcd_ksad01_definitions.csv", 
    delim = ",", escape_double = FALSE, col_types = "c", trim_ws = TRUE)
DKSAD2P <- read_delim("./WorkDir/ksads2daic_use_only_p01_definitions.csv", 
    delim = ",", escape_double = FALSE, col_types = "c", trim_ws = TRUE)
CBCLR <- read_delim("./WorkDir/abcd_cbcl01.txt", delim = "\t", escape_double = FALSE, col_types = "c",
    trim_ws = TRUE)[-1,]
CBCL <- read_delim("./WorkDir/abcd_cbcls01.txt", delim = "\t", escape_double = FALSE, col_types = "c",
    trim_ws = TRUE)[-1,]

# Recode KSAD1P values
KSAD1P[KSAD1P == "555"] <- "N"
KSAD1P[KSAD1P == "777"] <- "D"
KSAD1P[KSAD1P == "888"] <- "0"
KSAD1P[KSAD1P == "999"] <- "K"

KSAD2P[KSAD2P == "555"] <- "N"
KSAD2P[KSAD2P == "777"] <- "D"
KSAD2P[KSAD2P == "888"] <- "0"
KSAD2P[KSAD2P == "999"] <- "K"

CBCLR[CBCLR == "NA"] <- NA
CBCL[CBCL == "NA"] <- NA

# Remove noninformative columns, retain only useful variables in dictionaries
ColN <- function(x) apply(x, 2, function(row) unique(unique(row) == "N" & length(unique(row)) == 1))

KSAD1P[is.na(KSAD1P)] <- "N"
KSAD1P <- KSAD1P[,!ColN(KSAD1P)]
DKSAD1P <- DKSAD1P %>% filter(ElementName %in% colnames(KSAD1P))

KSAD2P[is.na(KSAD2P)] <- "N"
KSAD2P <- KSAD2P[,!ColN(KSAD2P)]
DKSAD2P <- DKSAD2P %>% filter(ElementName %in% colnames(KSAD2P))

# Subset OCD Variables
OCD1Ps <- c("ksads_11_917_p", "ksads_11_918_p", "ksads_11_919_p", "ksads_11_920_p")
OCD1P <- KSAD1P %>% dplyr::select(any_of(c(DKSAD1P$ElementName[c(1:4)],OCD1Ps)))

# Remove events with non-administered modules
OCD1P <- OCD1P %>% 
    unite(TEMP, ksads_11_917_p:ksads_11_920_p, remove = F, sep = "") %>%
    filter(!(TEMP == "NNNN")) %>%
    select(-TEMP)

# Calculate Lifetime OCD diagnoses
OCD <- OCD1P %>%
    select(c(subjectkey, ksads_11_917_p, ksads_11_918_p)) %>%
    mutate(LifetimeOCD = as.numeric(ksads_11_917_p) | as.numeric(ksads_11_918_p), .keep = "unused")  %>%
    group_by(subjectkey) %>%
    summarise(LifetimeOCD = sum(as.numeric(LifetimeOCD), na.rm=TRUE))

# Tabulate results
table(OCD$LifetimeOCD)

# Export results
write_delim(OCD, "./WorkDir/OCD.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

# Subset Tic Disorder Variables
MotorTics <- c("ksads2_17_99_p", "ksads2_17_100_p")
VocalTics <- c("ksads2_17_101_p", "ksads2_17_102_p")
TD <- KSAD2P %>% dplyr::select(any_of(c(DKSAD2P$ElementName[c(1,3,4,5)], MotorTics, VocalTics)))

# Remove events with non-administered modules
TD <- TD %>% 
    unite(TEMP, ksads2_17_99_p:ksads2_17_102_p, remove = F, sep = "") %>%
    filter(!(TEMP == "NNNN")) %>%
    select(-TEMP)

# Calculate Lifetime TD diagnoses
TD <- TD %>%
    mutate(LifetimeMotorTics = as.numeric(as.numeric(ksads2_17_99_p) | as.numeric(ksads2_17_100_p)), 
        LifetimeVocalTics = as.numeric(as.numeric(ksads2_17_101_p) | as.numeric(ksads2_17_102_p)), 
        .keep = "unused")  %>%
    rowwise() %>%
    mutate(LifetimeTD = sum(LifetimeMotorTics, LifetimeVocalTics, na.rm=TRUE))

# Tabulate results
table(TD$LifetimeTD)

# Save Tic Data
TD %>%
    select(c(subjectkey, LifetimeTD)) %>%
    write_delim("./WorkDir/TD.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

# Join in OCD Data
Pheno <- TD %>%
    left_join(OCD, by=c("subjectkey"="subjectkey"))

# Preprocess CBCL score dataset removing non-useful columns.
KEEP_CBCL <- CBCL %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    group_by(subjectkey) %>%
    summarise(n = n()) %>%
    filter(n == 3) %>%
    pull(subjectkey)

KEEP_CBCLR <- CBCLR %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    group_by(subjectkey) %>%
    summarise(n = n()) %>%
    filter(n == 3) %>%
    pull(subjectkey)

OCP <- CBCL %>%
    filter(subjectkey %in% KEEP_CBCL) %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(-c(collection_id, abcd_cbcls01_id, dataset_id, src_subject_id, interview_date, interview_age, sex, eventname, collection_title, ends_with(c("_m", "_t", "_nm")))) %>%
    mutate_at(., vars(contains("cbcl")), as.numeric) %>%
    select(subjectkey, cbcl_scr_07_ocd_r) %>%
    group_by(subjectkey) %>%
    summarise(OCP = sum(as.numeric(cbcl_scr_07_ocd_r), na.rm=TRUE)) %>%
    ungroup()

OCS <- CBCLR %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(any_of(c("subjectkey", "sex", "cbcl_q09_p", "cbcl_q66_p"))) %>%
    group_by(subjectkey) %>%
    summarise(OBS = sum(as.numeric(cbcl_q09_p), na.rm=TRUE),
        COM = sum(as.numeric(cbcl_q66_p), na.rm=TRUE)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(OCS = sum(c(OBS,COM), na.rm=TRUE))

# Join in with OCD data and run regressions
MASTER <- full_join(OCS, OCP, by=c("subjectkey"="subjectkey")) %>%
    full_join(OCD, by=c("subjectkey"="subjectkey")) %>%
    mutate(nOCD=ifelse(LifetimeOCD == 2, 1, 0),
    bnOCD=ifelse(LifetimeOCD >= 1, 1, 0))

# Save CBCL Data
MASTER %>%
    select(c(subjectkey, OBS)) %>%
    write_delim("./WorkDir/OBS.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

MASTER %>%
    select(c(subjectkey, COM)) %>%
    write_delim("./WorkDir/COM.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

MASTER %>%
    select(c(subjectkey, OCS)) %>%
    write_delim("./WorkDir/OCS.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

MASTER %>%
    select(c(subjectkey, OCP)) %>%
    write_delim("./WorkDir/OCP.pheno", delim = "\t", col_names = FALSE, quote = "none", eol = "\n")

# Exit out of R and don't save workspace
quit()
n
```
Based on the 4th release of ABCD cohort, the OCD distribution is as follows:

0. OCD Negative: 10,268
1. Broad OCD: 1,276
2. Narrow OCD: 309

```bash
# Unload R module
module unload R

# Unload UFRC module after session timeout
module unload UFRC
```
## 02.01. Check Call Rates for Samples and Markers

```bash
# Create a holding folder
mkdir -p ${WD}/0201

# Calculate missingness
plink --bfile ${WD}/ABCD_release_3.0_QCed \
    --missing \
    --out ${WD}/0201/ABCD

# Unload plink, load R and launch R module
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import genotyping frequency files
Samples <- read.csv("./0201/ABCD.imiss", sep = "")
Markers <- read.csv("./0201/ABCD.lmiss", sep = "")

# Update datasets with genotyping rates
Samples <- Samples %>%
    mutate(GF = 1 - F_MISS)
Markers <- Markers %>%
    mutate(GF = 1 - F_MISS)

# Plot a histogram of genotyping frequency in all markers and samples
PlotSamples <- ggplot(Samples, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=0.95, color="#e63946", linetype="dashed") +
    labs(title="Sample genotyping frequencies",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 11,099") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureA.png", path="./0201/",
    width=400, height=150, units="mm", dpi=600)

PlotMarkers <- ggplot(Markers, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=0.95, color="#e63946", linetype="dashed") +
    labs(title="Marker genotyping frequencies",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 516,598") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotMarkers,
    file="FigureB.png", path="./0201/",
    width=400, height=150, units="mm", dpi=600)

# Import batch and phenotype information
Batch <- read.delim("./ABCD_release3.0_.batch_info.txt")
Pheno <- read.delim("./OCD.pheno", header=F)

# Update column names for loaded dataframes
colnames(Batch) <- c("FID", "IID", "Plate", "Batch")
colnames(Pheno) <- c("IID", "LifetimeOCD")

# Update samples with batch and pheno information
SamplesPlus <- Samples %>%
    left_join(Batch, by=c("FID"="FID", "IID"="IID")) %>%
    left_join(Pheno, by=c("IID"="IID"))

# Tabulate results
table(SamplesPlus$Batch, SamplesPlus$LifetimeOCD)
table(SamplesPlus$Batch)

# Plot batchwise sample genotyping frequency
PlotSamplesBatch <- ggplot(SamplesPlus, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=0.95, color="#e63946", linetype="dashed") +
    labs(title="Sample genotyping frequencies, by batches",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 11,099") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000")) +
    facet_wrap(vars(Batch), 
        ncol=2, nrow=6,
        scales="free")

ggsave(plot=PlotSamplesBatch,
    file="FigureC.png", path="./0201/",
    width=400, height=250, units="mm", dpi=600)

# Summarize genotyping frequencies by batch
tapply(SamplesPlus$GF, SamplesPlus$Batch, summary)

# Tabulate OCD cases in bad batches (WB batches, indices 2 and 4)
table(SamplesPlus$Batch, SamplesPlus$LifetimeOCD)[c(2,4),]

# Extract individuals from whole blood batches and mark them for removal
SamplesPlus %>%
    filter(Batch %in% c("BATCH_1_WB", "BATCH_2_WB")) %>%
    select(c(FID, IID)) %>%
    write.table(file="./0201/BatchRemove.txt",
        row.names=F, col.names=F, quote=F,
        sep="\t", eol="\n")

# Exit out of R and don't save workspace
quit()
n
```

## 02.02. Remove Problematic IDs 
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Create a holding folder
mkdir -p ${WD}/0202

# Remove problematic batches
plink --bfile ${WD}/ABCD_release_3.0_QCed \
    --remove ${WD}/0201/BatchRemove.txt \
    --missing \
    --make-bed \
    --out ${WD}/0202/ABCD

# Unload plink, load R and launch R module
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Load fam file
Fam <- read.delim("./0202/ABCD.fam", header=F, sep = "")

# Check if all rows in FID column have consistend form according to ABCD study design
nrow(Fam) # 10,892
Fam %>%
    filter(str_detect(V1, "^AB")) %>%
    nrow() # 10,892

# Check if all rows in IID column have consistent form according to ABCD study design
nrow(Fam) # 10,892
Fam %>%
    filter(str_detect(V2, "^NDAR_")) %>%
    nrow() # 10,890

# Identify problematic IIDs
Fam %>%
    filter(!str_detect(V2, "^NDAR_"))

# Save the problematic IDs and their corrections
Fam %>%
    filter(V2 %in% c("`NDAR_INVF3FYXH1G", "NDARINVPWLFYWTX")) %>%
    select(c(V1, V2)) %>%
    mutate(V1.2=V1,
        VV2.2=c("NDAR_INVF3FYXH1G", "NDAR_INVPWLFYWTX")) %>%
    write.table(file="./0202/SampleRename.txt",
        row.names=F, col.names=F, quote=F,
        sep="\t", eol="\n")

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Remove problematic batches
plink --bfile ${WD}/0202/ABCD \
    --update-ids ${WD}/0202/SampleRename.txt \
    --make-bed \
    --out ${WD}/0202/ABCD2
```

## 02.03. Check Sex
```bash
# Create a holding folder
mkdir -p ${WD}/0203

# Unload plink, load R and launch R module
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Load fam file
Fam <- read.delim("./0202/ABCD2.fam", header=F, sep = "")

# Load in a pheno file and extract sex information
KSAD1P <- read_delim("./abcd_ksad01.txt", delim = "\t", escape_double = FALSE, col_types = "c",
    trim_ws = TRUE)[-1,]

KSAD1P <- KSAD1P %>%
    select(c(subjectkey, sex)) %>%
    distinct()

# Merge sex information into Fam file
Fam <- Fam %>%
    left_join(KSAD1P, by=c("V2"="subjectkey")) %>%
    mutate(V5=sex) %>%
    select(-sex) %>%
    mutate(V5=recode(.$V5, "M"="1", "F"="2", .missing="0")) %>%
    select(c(V1, V2, V5)) %>%
    write.table(file="./0203/SexUpdate.txt",
        row.names=F, col.names=F, quote=F,
        sep="\t", eol="\n")

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Update FAM file with sex info
plink --bfile ${WD}/0202/ABCD2 \
    --update-sex ${WD}/0203/SexUpdate.txt \
    --make-bed \
    --out ${WD}/0203/ABCD

# Impute sex
plink --bfile ${WD}/0203/ABCD \
    --check-sex \
    --out ${WD}/0203/ABCD

# Flag IDs with problematic sex information
awk '$5 ~ /PROBLEM/' ${WD}/0203/ABCD.sexcheck > ${WD}/0203/ABCD.flag.sexcheck
awk '{$3="PROBLEM"}1' ${WD}/0203/ABCD.nosex > ${WD}/0203/ABCD.flag.nosex

join -j 2 -a 1 -e NULL ${WD}/0203/ABCD.flag.sexcheck \
    ${WD}/0203/ABCD.flag.nosex > ${WD}/0203/ABCD.flag.sexprob

# Create an exclusion file based on sex problematic IDs
awk 'BEGIN {OFS="\t"; ORS="\n"}
    {print $2, $1}' ${WD}/0203/ABCD.flag.sexprob > \
    ${WD}/0203/SexOut.txt

# Remove IDs with problematic sex information
plink --bfile ${WD}/0203/ABCD \
    --remove ${WD}/0203/SexOut.txt \
    --check-sex \
    --make-bed \
    --out ${WD}/0203/ABCD2
```

## 02.04. Remove non-autosomes and low call frequency markers, and samples

```bash
# Create a holding folder
mkdir -p ${WD}/0204

# Remove nonautosomes and calculate genotyping frequencies
plink --bfile ${WD}/0203/ABCD2 \
    --chr 1-22 \
    --make-bed \
    --out ${WD}/0204/ABCD94

plink --bfile ${WD}/0204/ABCD94 \
    --missing \
    --out ${WD}/0204/ABCD94

# Serial filtering
for i in {94..97};
do
    j=$((i+1))
    echo "Start: $i, End: $j"

    cutoff="0.0$((100-j))"

    # SNP
    plink --bfile ${WD}/0204/ABCD${i} \
        --geno ${cutoff} \
        --make-bed \
        --out ${WD}/0204/ABCD_SNP${j};

    # Sample
    plink --bfile ${WD}/0204/ABCD_SNP${j} \
        --mind ${cutoff} \
        --make-bed \
        --out ${WD}/0204/ABCD${j};
done

plink --bfile ${WD}/0204/ABCD98 \
    --missing \
    --out ${WD}/0204/ABCD98

# Unload plink, load R and launch R module
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Load missingness files
Samples <- read.csv("./0204/ABCD98.imiss", sep = "")
Markers <- read.csv("./0204/ABCD98.lmiss", sep = "")

# Update datasets with genotyping rates
Samples <- Samples %>%
    mutate(GF = 1 - F_MISS)
Markers <- Markers %>%
    mutate(GF = 1 - F_MISS)

# Plot a histogram of genotyping frequency in all markers and samples
PlotSamples <- ggplot(Samples, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=c(0.95, 0.98), color="#e63946", linetype="dashed") +
    labs(title="Sample genotyping frequencies, post-QC",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 10,607") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureD.png", path="./0204/",
    width=400, height=150, units="mm", dpi=600)

PlotMarkers <- ggplot(Markers, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=c(0.95, 0.98), color="#e63946", linetype="dashed") +
    labs(title="Marker genotyping frequencies, post-QC",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 480,427") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotMarkers,
    file="FigureE.png", path="./0204/",
    width=400, height=150, units="mm", dpi=600)

# Import batch and phenotype information
Batch <- read.delim("./ABCD_release3.0_.batch_info_corrected.txt")
Pheno <- read.delim("./OCD.pheno", header=F)

# Update column names for loaded dataframes
colnames(Batch) <- c("FID", "IID", "Plate", "Batch")
colnames(Pheno) <- c("IID", "LifetimeOCD")

# Update samples with batch and pheno information
SamplesPlus <- Samples %>%
    left_join(Batch, by=c("FID"="FID", "IID"="IID")) %>%
    left_join(Pheno, by=c("IID"="IID"))

# Tabulate results
table(SamplesPlus$Batch, SamplesPlus$LifetimeOCD)
table(SamplesPlus$Batch)

# Plot batchwise sample genotyping frequency
PlotSamplesBatch <- ggplot(SamplesPlus, aes(x = GF)) +
    geom_histogram(breaks=seq(0.89, 1.0, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=c(0.95, 0.98), color="#e63946", linetype="dashed") +
    labs(title="Sample genotyping frequencies, by batches, post-QC",
        x = "Genotyping Rate",
        y = "Frequency",
        caption = "N = 10,607") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000")) +
    facet_wrap(vars(Batch), 
        ncol=2, nrow=4,
        scales="free")

ggsave(plot=PlotSamplesBatch,
    file="FigureF.png", path="./0204/",
    width=400, height=200, units="mm", dpi=600)

# Summarize genotyping frequencies by batch
tapply(SamplesPlus$GF, SamplesPlus$Batch, summary)

# Exit out of R and don't save workspace
quit()
n
```

# Migrate into parent directory
cd ${WD}

# Load in Plink module
module load plink/1.90b3.39

# Create a holding folder
mkdir -p ${WD}/0301
```
## 03.01. Prune SNPs
```bash
# Create a list of SNPs to be pruned
plink --bfile ${WD}/0204/ABCD98 \
      --indep-pairwise 200 50 0.15 \
      --out ${WD}/0301/ABCD

# Prune the SNPs
plink --bfile ${WD}/0204/ABCD98 \
      --extract ${WD}/0301/ABCD.prune.in \
      --make-bed \
      --out ${WD}/0301/ABCD.Pruned
```

## 03.02. IBD QC

```bash
# Create a holding folder
mkdir -p ${WD}/0302

# Run PLINK IBD estimation
plink --bfile ${WD}/0301/ABCD.Pruned \
      --genome \
      --min 0.15 \
      --out ${WD}/0302/ABCD.Pruned

# Get missingness report for selections
plink --bfile ${WD}/0301/ABCD.Pruned \
    --missing \
    --out ${WD}/0302/ABCD.Pruned

# Unload Plink, load and open R modules
module unload plink/1.90b3.39
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import datasets
Genome <- read.csv("./0302/ABCD.Pruned.genome", sep="", header = TRUE)
CallFreq <- read.csv("./0302/ABCD.Pruned.imiss", sep="")

# Identify contaminated samples
Freq <- rbind(setNames(select(Genome, FID1, IID1), c("FID", "IID")), 
    setNames(select(Genome, FID2, IID2), c("FID", "IID"))) %>%
    group_by(IID, FID) %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    filter(n > 50)

# Save the list of files for removal
Freq %>%
    select(c(FID, IID)) %>%
    write.table(file="./0302/GenoOut1.txt", sep="\t", eol="\n", 
        row.names=F, col.names=F, quote=F)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Remove contaminated individuals
plink --bfile ${WD}/0301/ABCD.Pruned \
    --remove ${WD}/0302/GenoOut1.txt \
    --make-bed \
    --out ${WD}/0302/ABCD2

# Run PLINK IBD estimation on new bed files
plink --bfile ${WD}/0302/ABCD2 \
      --genome \
      --min 0.15 \
      --out ${WD}/0302/ABCD2

# Unload Plink, load and open R modules
module unload plink/1.90b3.39
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import datasets
Genome <- read.csv("./0302/ABCD2.genome", sep="", header = TRUE)
CallFreq <- read.csv("./0302/ABCD.Pruned.imiss", sep="")

# Identify contaminated samples
Freq <- rbind(setNames(select(Genome, FID1, IID1), c("FID", "IID")), 
    setNames(select(Genome, FID2, IID2), c("FID", "IID"))) %>%
    group_by(IID, FID) %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    filter(n > 10)

# Save the list of files for removal
Freq %>%
    select(c(FID, IID)) %>%
    write.table(file="./0302/GenoOut2.txt", sep="\t", eol="\n", 
        row.names=F, col.names=F, quote=F)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Remove contaminated individuals
plink --bfile ${WD}/0302/ABCD2 \
    --remove ${WD}/0302/GenoOut2.txt \
    --make-bed \
    --out ${WD}/0302/ABCD3

# Run PLINK IBD estimation on new bed files
plink --bfile ${WD}/0302/ABCD3 \
      --genome \
      --min 0.15 \
      --out ${WD}/0302/ABCD3

# Unload Plink, load and open R modules
module unload plink/1.90b3.39
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import datasets
Genome <- read.csv("./0302/ABCD3.genome", sep="", header = TRUE)
CallFreq <- read.csv("./0302/ABCD.Pruned.imiss", sep="")
Pheno <- read.delim("./OCD.pheno", header=F)

# Identify contaminated samples
Freq <- rbind(setNames(select(Genome, FID1, IID1), c("FID", "IID")), 
    setNames(select(Genome, FID2, IID2), c("FID", "IID"))) %>%
    group_by(IID) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

# Prepare dataframe for analysis
Genome <- Genome %>%
    select(c("FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2", "PI_HAT", "PHE")) %>%
    left_join(CallFreq, by=c("FID1"="FID", "IID1"="IID")) %>%
    select(-c("MISS_PHENO", "N_MISS", "N_GENO")) %>%
    rename(MISS1=F_MISS) %>%
    left_join(CallFreq, by=c("FID2"="FID", "IID2"="IID")) %>%
    select(-c("MISS_PHENO", "N_MISS", "N_GENO")) %>%
    rename(MISS2=F_MISS) %>%
    left_join(Pheno, by=c("IID1"="V1")) %>%
    rename(PHE1=V2) %>%
    left_join(Pheno, by=c("IID2"="V1")) %>%
    rename(PHE2=V2) %>%
    left_join(Freq, by=c("IID1"="IID")) %>%
    rename(N1=n) %>%
    left_join(Freq, by=c("IID2"="IID")) %>%
    rename(N2=n)

# Change NA to -9
Genome$PHE1[is.na(Genome$PHE1)] <- -9
Genome$PHE2[is.na(Genome$PHE2)] <- -9

# Merge PHE columns for analysis
Genome <- Genome %>%
    mutate(PHE=paste(PHE1,PHE2, sep=""))

# Define a storage dataframe for sample removal
Out <- data.frame(FID = as.character(NULL),
    IID = as.character(NULL))

# A: Extract -9-9, 00, 11, and 22 related samples
GenA <- Genome %>%
    filter(PHE %in% c("-9-9", "00", "11", "22")) %>%
    mutate(Remove_N = case_when(
            N2 == N1 ~ NA_real_,
            N2 > N1 ~ 2,
            N2 < N1 ~ 1),
        Remove_MISS = case_when(
            MISS2 == MISS1 ~ NA_real_,
            MISS2 > MISS1 ~ 2,
            MISS2 < MISS1 ~ 1)) %>%
    mutate(REMOVE = case_when(
        Remove_N == 1 ~ 1,
        Remove_N == 2 ~ 2,
        is.na(Remove_N) & Remove_MISS == 1 ~ 1,
        is.na(Remove_N) & Remove_MISS == 2 ~ 2,
        TRUE ~ 1))

for(i in 1:nrow(GenA)){
    if (GenA$REMOVE[i] == 1) {
        OUT <- GenA[i, c(1, 2)]
    } else if (GenA$REMOVE[i] == 2) {
        OUT <- GenA[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

# B: Extract -90, 0-9, -91, 1-9, -92, 2-9 samples
GenB <- Genome %>%
    filter(PHE %in% c("-90", "0-9", "-91", "1-9", "-91",  "1-9", "-92", "2-9")) %>%
    mutate(REMOVE = case_when(
            PHE1 == -9 ~ 1,
            PHE2 == -9 ~ 2))

for(i in 1:nrow(GenB)){
    if (GenB$REMOVE[i] == 1) {
        OUT <- GenB[i, c(1, 2)]
    } else if (GenB$REMOVE[i] == 2) {
        OUT <- GenB[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

# C: Extract 01, 10, 02, 20
GenC <- Genome %>%
    filter(PHE %in% c("01", "10", "02", "20")) %>%
    mutate(REMOVE = case_when(
            PHE1 == 0 ~ 1,
            PHE2 == 0 ~ 2))

for(i in 1:nrow(GenC)){
    if (GenC$REMOVE[i] == 1) {
        OUT <- GenC[i, c(1, 2)]
    } else if (GenC$REMOVE[i] == 2) {
        OUT <- GenC[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

# D: Extract 12, 21
GenD <- Genome %>%
    filter(PHE %in% c("21", "12")) %>%
    mutate(REMOVE = case_when(
            PHE1 == 1 ~ 1,
            PHE2 == 1 ~ 2))

for(i in 1:nrow(GenD)){
    if (GenD$REMOVE[i] == 1) {
        OUT <- GenD[i, c(1, 2)]
    } else if (GenD$REMOVE[i] == 2) {
        OUT <- GenD[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

Out <- unique(Out)

# Save the list of samples to remove
write.table(Out, file="./0302/GenoOut3.txt", sep="\t", eol="\n", 
    row.names=F, col.names=F, quote=F)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Remove contaminated individuals
plink --bfile ${WD}/0302/ABCD3 \
    --remove ${WD}/0302/GenoOut3.txt \
    --make-bed \
    --out ${WD}/0302/ABCD4

# Run PLINK IBD estimation on new bed files
plink --bfile ${WD}/0302/ABCD4 \
    --genome \
    --min 0.15 \
    --out ${WD}/0302/ABCD4

# Unload Plink, load and open R modules
module unload plink/1.90b3.39
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Import datasets
Genome <- read.csv("./0302/ABCD4.genome", sep="", header = TRUE)
CallFreq <- read.csv("./0302/ABCD.Pruned.imiss", sep="")
Pheno <- read.delim("./OCD.pheno", header=F)

# Import datasets
Genome <- read.csv("./0302/ABCD3.genome", sep="", header = TRUE)
CallFreq <- read.csv("./0302/ABCD.Pruned.imiss", sep="")
Pheno <- read.delim("./OCD.pheno", header=F)

# Identify contaminated samples
Freq <- rbind(setNames(select(Genome, FID1, IID1), c("FID", "IID")), 
    setNames(select(Genome, FID2, IID2), c("FID", "IID"))) %>%
    group_by(IID) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

# Prepare dataframe for analysis
Genome <- Genome %>%
    select(c("FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2", "PI_HAT", "PHE")) %>%
    left_join(CallFreq, by=c("FID1"="FID", "IID1"="IID")) %>%
    select(-c("MISS_PHENO", "N_MISS", "N_GENO")) %>%
    rename(MISS1=F_MISS) %>%
    left_join(CallFreq, by=c("FID2"="FID", "IID2"="IID")) %>%
    select(-c("MISS_PHENO", "N_MISS", "N_GENO")) %>%
    rename(MISS2=F_MISS) %>%
    left_join(Pheno, by=c("IID1"="V1")) %>%
    rename(PHE1=V2) %>%
    left_join(Pheno, by=c("IID2"="V1")) %>%
    rename(PHE2=V2) %>%
    left_join(Freq, by=c("IID1"="IID")) %>%
    rename(N1=n) %>%
    left_join(Freq, by=c("IID2"="IID")) %>%
    rename(N2=n)

# Change NA to -9
Genome$PHE1[is.na(Genome$PHE1)] <- -9
Genome$PHE2[is.na(Genome$PHE2)] <- -9

# Merge PHE columns for analysis
Genome <- Genome %>%
    mutate(PHE=paste(PHE1,PHE2, sep=""))

# Define a storage dataframe for sample removal
Out <- data.frame(FID = as.character(NULL),
    IID = as.character(NULL))

# Pull out Individuals who are matched with a lot of people
Pull <- Freq %>% 
    filter(n > 2) %>% 
    pull(IID)

Pull <- CallFreq %>%
    filter(IID %in% Pull) %>%
    select(c(FID, IID))

Out <- rbind(Out, Pull)

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

# A: Extract -9-9, 00, 11, and 22 related samples
GenA <- Genome %>%
    filter(PHE %in% c("-9-9", "00", "11", "22")) %>%
    mutate(Remove_N = case_when(
            N2 == N1 ~ NA_real_,
            N2 > N1 ~ 2,
            N2 < N1 ~ 1),
        Remove_MISS = case_when(
            MISS2 == MISS1 ~ NA_real_,
            MISS2 > MISS1 ~ 2,
            MISS2 < MISS1 ~ 1)) %>%
    mutate(REMOVE = case_when(
        Remove_N == 1 ~ 1,
        Remove_N == 2 ~ 2,
        is.na(Remove_N) & Remove_MISS == 1 ~ 1,
        is.na(Remove_N) & Remove_MISS == 2 ~ 2,
        TRUE ~ 1))

for(i in 1:nrow(GenA)){
    if (GenA$REMOVE[i] == 1) {
        OUT <- GenA[i, c(1, 2)]
    } else if (GenA$REMOVE[i] == 2) {
        OUT <- GenA[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

# C: Extract 01, 10, 02, 20
GenC <- Genome %>%
    filter(PHE %in% c("01", "10", "02", "20")) %>%
    mutate(REMOVE = case_when(
            PHE1 == 0 ~ 1,
            PHE2 == 0 ~ 2))

for(i in 1:nrow(GenC)){
    if (GenC$REMOVE[i] == 1) {
        OUT <- GenC[i, c(1, 2)]
    } else if (GenC$REMOVE[i] == 2) {
        OUT <- GenC[i, c(3, 4)]
    }
    colnames(OUT) <- c("FID", "IID")
    Out <- rbind(Out, OUT)
}

Genome <- Genome %>%
    filter(!(IID1 %in% Out$IID)) %>%
    filter(!(IID2 %in% Out$IID))

Out <- unique(Out)

# Save the list of samples to remove
write.table(Out, file="./0302/GenoOut4.txt", sep="\t", eol="\n", 
    row.names=F, col.names=F, quote=F)

# Exit out of R and don't save workspace
quit()
n
```

## 03.03. Keep non-related samples

```bash
# Unload R module and load in Plink module
module unload R
module load plink/1.90b3.39

# Remove contaminated individuals
plink --bfile ${WD}/0302/ABCD4 \
    --remove ${WD}/0302/GenoOut4.txt \
    --make-bed \
    --out ${WD}/0302/ABCD5

# Run PLINK IBD estimation on new bed files
plink --bfile ${WD}/0302/ABCD5 \
    --genome \
    --min 0.15 \
    --out ${WD}/0302/ABCD5

# Create a new directory
mkdir ${WD}/0303

# Keep the passing samples in non-pruned files
awk 'NR>=1 {OFS="\t"; ORS="\n"}
    {print $1, $2}' ${WD}/0302/ABCD5.fam > ${WD}/0303/ABCDFinalIn.txt

plink --bfile ${WD}/0204/ABCD98 \
      --keep ${WD}/0303/ABCDFinalIn.txt \
      --make-bed \
      --out ${WD}/0303/ABCD
```

## 04.01 Hardy-Weinberg Equilibrium QC

```bash
# Make a holding directory
mkdir -p ${WD}/0401

# Calculate Hardy-Weinberg Equilibrium Statistics
plink --bfile ${WD}/0303/ABCD \
      --hardy \
      --out ${WD}/0401/ABCD

# Pull out strongly deviating SNPs
awk '{if ($9 <0.00001) print $0}' \
    ${WD}/0401/ABCD.hwe > \
    ${WD}/0401/ABCD_Diseq.hwe

# Load and open R, unload Plink
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Plot HWE
# Import HWE reports
HWE <- read.csv("./0401/ABCD.hwe", sep="", header = TRUE)
HWE_Diseq <- read.csv("./0401/ABCD_Diseq.hwe", sep="", header = FALSE)
colnames(HWE_Diseq) <- colnames(HWE)

# Plot all HWE
PlotSamples <- ggplot(HWE, aes(x = P)) +
    geom_histogram(breaks=seq(1, 0, -0.01), color="#1d3557", fill="#457b9d") +
    #geom_vline(xintercept=c(0.95, 0.98), color="#e63946", linetype="dashed") +
    labs(title="Hardy-Weinberg Equilibrium, All",
        x = "HWE P-Value",
        y = "Frequency",
        caption = "N = 8,718") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureG.png", path="./0401/",
    width=400, height=150, units="mm", dpi=600)

# Plot all HWD
PlotSamples <- ggplot(HWE_Diseq, aes(x = -log10(P+10^(-100)))) +
    geom_histogram(breaks=seq(6, 100, 1), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=c(7), color="#e63946", linetype="dashed") +
    labs(title="Hardy-Weinberg Equilibrium, P < e-6",
        x = "HWE -log10 P-Value",
        y = "Frequency",
        caption = "N = 8,718") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureH.png", path="./0401/",
    width=400, height=150, units="mm", dpi=600)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Load in Plink module, unload R
module unload R
module load plink/1.90b3.39

# Filter out HWE < 10-7
plink --bfile ${WD}/0303/ABCD \
    --hwe 1e-7 \
    --hwe-all \
    --make-bed \
    --out ${WD}/0401/ABCD
```

## 04.02 Heterozygosity
```bash
# Make directory
mkdir ${WD}/0402

# Eclude high inversion region and prune
## 6 25500000    33500000    HLA
## 8 8135000 12000000    Inversion8
## 17    40900000    45000000    Inversion17
plink --bfile ${WD}/0401/ABCD \
    --exclude ${WD}/0402/InversionExclude.txt \
    --range \
    --indep-pairwise 50 5 0.2 \
    --out ${WD}/0402/ABCD_Pruned

plink --bfile ${WD}/0401/ABCD \
    --extract ${WD}/0402/ABCD_Pruned.prune.in \
    --het \
    --out ${WD}/0402/ABCD_Het

# Load and open R, unload Plink
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Plot Heterozygosity
# Import Het reports
Het <- read.csv("./0402/ABCD_Het.het", sep="", header = TRUE)

Het <- Het %>%
    mutate(HET_RATE = (N.NM. - O.HOM.)/N.NM.)

# Find min and max rate cutoffs
HetMin <- mean(Het$HET_RATE) - 4*sd(Het$HET_RATE)
HetMax <- mean(Het$HET_RATE) + 4*sd(Het$HET_RATE)

# Plot all Het
PlotSamples <- ggplot(Het, aes(x = HET_RATE)) +
    geom_histogram(breaks=seq(0.17, 0.24, 0.0005), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=c(HetMin, HetMax), color="#e63946", linetype="dashed") +
    labs(title="Heterozygosity",
        x = "Heterozygosity Rate",
        y = "Frequency",
        caption = "N = 8,718") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureI.png", path="./0402/",
    width=400, height=150, units="mm", dpi=600)

# Mark samples for exclusions
Het %>%
    filter(HET_RATE < HetMin | HET_RATE > HetMax) %>%
    select(c(FID, IID)) %>%
    write.table(file="./0402/HetOut.txt", sep="\t", eol="\n", 
    row.names=F, col.names=F, quote=F)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Load in Plink module
module unload R
module load plink/1.90b3.39

# Remove het
plink --bfile ${WD}/0401/ABCD \
    --remove ${WD}/0402/HetOut.txt \
    --make-bed \
    --out  ${WD}/0402/ABCD
```

## 04.03 Minor Allele Frequency QC

```bash
# Make a holding directory
mkdir -p ${WD}/0403

# Calculate and plot sample MAF
plink --bfile ${WD}/0402/ABCD \
    --freq \
    --out ${WD}/0403/ABCD

# Load and open R, unload Plink
module unload plink
module load R
R
```
```r
# Set working directory
setwd("...")

# Load necessary libraries
library(tidyverse)

# Plot Heterozygosity
# Import Het reports
MAF <- read.csv("./0403/ABCD.frq", sep="", header = TRUE)

# Plot all Het
PlotSamples <- ggplot(MAF, aes(x = MAF)) +
    geom_histogram(breaks=seq(0, 0.5, 0.001), color="#1d3557", fill="#457b9d") +
    geom_vline(xintercept=seq(0.01, 0.05, 0.01), color="#e63946", linetype="dashed") +
    labs(title="Minor Allele Frequency",
        x = "MAF",
        y = "Frequency",
        caption = "N = 8,713") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(),
        legend.position="none",
        axis.text.y=element_text(angle=90, hjust=0.3),
        strip.text.x=element_blank(),
        plot.margin=margin(20,7,7,7,"pt"),
        panel.border=element_blank(),
        axis.line=element_line(color="#000000"))

ggsave(plot=PlotSamples,
    file="FigureJ.png", path="./0403/",
    width=400, height=150, units="mm", dpi=600)

# Exit out of R and don't save workspace
quit()
n
```
```bash
# Load in Plink module
module unload R
module load plink/1.90b3.39

# Remove MAF < 0.05
plink --bfile ${WD}/0402/ABCD \
    --maf 0.05 \
    --make-bed \
    --out ${WD}/0403/ABCD

# N Variants: 294757
# N Individuals: 8713
```

## 05.00 Download and Prepare Reference Files
```bash
# Create a holding folders
mkdir -p ${RF}/SHAPEIT4MAPS ${RF}/1000Gp3hg37

###############################
## Get 1000Genomes Reference ##
###############################
cd ${RF}/1000Gp3hg37

for chr in {1..22};
do
    wget  http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --output-document=chr${chr}.vcf.gz;
done

# Load BCFTOOLS module bcftools/1.13 and filter files for SNPs and MAF
module load bcftools/1.13

for chr in {1..22};
do
    bcftools view \
    --include 'VT="SNP" && EAS_AF>=0.001 && EUR_AF>=0.001 && AFR_AF>=0.001 && AMR_AF>=0.001 && SAS_AF>=0.001' \
    --output-type z \
    --output chr${chr}.vcf.gz \
    --threads 5 \
    chr${chr}.vcf.gz;

    bcftools index \
    chr${chr}.vcf.gz;
done

###################
## Get MAP Files ##
###################
cd ${RF}/SHAPEIT4MAPS

wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz |
    --output-document=genetic_maps.b37.tar.gz

tar xvf genetic_maps.b37.tar.gz
rm genetic_maps.b37.tar.gz
```

## 05.01 Change Plink files to VCF files

```bash
# Set up directory
mkdir -p ${WD}/0501

# Load in Plink modul
module load plink/1.90b3.39

# Filter out SNPs which have problematic format or are indels
awk '($5=="A" || $5=="T" || $5=="G" || $5=="C") && \
    ($6=="A" || $6=="T" || $6=="G" || $6=="C") \
    {print $2}' ${WD}/0403/ABCD.bim > \
    ${WD}/0501/KeepSNPs.txt

# Convert plink files to VCF format and split into individual chromosomes
for chr in {1..22};
do
plink --bfile ${WD}/0403/ABCD \
    --chr ${chr} \
    --extract ${WD}/0501/KeepSNPs.txt \
    --recode vcf-iid bgz \
    --threads 5 \
    --out ${WD}/0501/ABCD.chr${chr};
done

# Load in BCFtools
module load bcftools/1.13

# Convert vcf.gz to bcf.gz, index files, remove vcf.gz files
for chr in {1..22};
do
    bcftools index \
    ${WD}/0501/ABCD.chr${chr}.vcf.gz;
done
```

## 05.02 Phase data with SHAPEIT4
```bash
# Make a holding directory
mkdir -p ${WD}/0502

# Run SHAPEIT4 on the data, but first exit out of dev environment
cd ${WD}/0502
sbatch ${SLURM}/PhaseITwithSHAPEIT.sh
```

### PhaseITwithSHAPEIT.sh
```bash
#!/bin/bash
#SBATCH --job-name=PhaseITwithShapeIT4                # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=njofrica@ufl.edu                # Where to send mail
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=12                           # Number of CPU cores per task
#SBATCH --mem=25gb                                   # Total memory limit
#SBATCH --time=24:00:00                             # Time limit hrs:min:sec
#SBATCH --output=.../SLURM_Scripts/LogFiles/SLURM_%j-%x.log   # Standard output and error log

date;hostname;pwd

export OMP_NUM_THREADS=12

#----------------------------------------------------------------------
# Load necessary modules
#----------------------------------------------------------------------
module load shapeit4/4.2.2

#----------------------------------------------------------------------
# Assign paths to variable names go to WorkDir
#----------------------------------------------------------------------
PD="..."
WD="..."
SLURM="..."
LOGS="..."
RF="..."

mkdir -p ${WD}/0502
cd ${WD}/0502

#----------------------------------------------------------------------
# Phase data with SHAPEIT4
#----------------------------------------------------------------------
for i in {1..22};
do
shapeit4 --input ${WD}/0501/ABCD.chr${i}.vcf.gz \
          --map ${RF}/SHAPEIT4MAPS/chr${i}.b37.gmap.gz \
          --region ${i} \
          --reference ${RF}/1000Gp3hg37/chr${i}.vcf.gz \
          --use-PS 0.0001 \
          --thread 12 \
          --log ${WD}/0502/ABCD.chr${i}.phased.log \
          --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
          --pbwt-depth 8 \
          --output ${WD}/0502/ABCD.chr${i}.phased.vcf
done

echo "Finished."
```

## 05.03 Impue data with IMPUTE5

```bash
# Make a holding directory
mkdir -p ${WD}/0503

# Copy IMPUTE5 executables to working directory, until module is fixed
cp /home/njofrica/impute5_v1.1.5/impute5_1.1.5_static 0503/impute5_1.1.5_static
cp /home/njofrica/impute5_v1.1.5/imp5Converter_1.1.5_static 0503/imp5Converter_1.1.5_static
cp /home/njofrica/impute5_v1.1.5/imp5Chunker_1.1.5_static 0503/imp5Chunker_1.1.5_static

# Run SHAPEIT4 on the data, but first exit out of dev environment
cd ${WD}/0502
sbatch ${SLURM}/PhaseITwithSHAPEIT.sh
```

### PhaseITwithSHAPEIT.sh
```
#!/bin/bash
#SBATCH --job-name=ImputeITwithIMPUTE5              # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=njofrica@ufl.edu                # Where to send mail
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=12                           # Number of CPU cores per task
#SBATCH --mem=30gb                                   # Total memory limit
#SBATCH --time=7-00:00:00                             # Time limit hrs:min:sec
#SBATCH --output=.../SLURM_Scripts/LogFiles/SLURM_%j-%x.log   # Standard output and error log

date;hostname;pwd

export OMP_NUM_THREADS=12

#----------------------------------------------------------------------
# Load necessary modules
#----------------------------------------------------------------------
module load bcftools/1.13

#----------------------------------------------------------------------
# Pull in IMPUTE5 software
#----------------------------------------------------------------------
CHUNKER="/imp5Chunker_1.1.5_static"
CONVERTER="/imp5Converter_1.1.5_static"
IMPUTER="/impute5_1.1.5_static"

#----------------------------------------------------------------------
# Assign paths to variable names go to WorkDir
#----------------------------------------------------------------------
PD="..."
WD="..."
SLURM="..."
LOGS="..."
RF="..."

mkdir -p ${WD}/0503
cd ${WD}/0503

#----------------------------------------------------------------------
# Impute data with IMPUTE5
#----------------------------------------------------------------------
for i in {21..22};
do
    #------------------------------------------------------------------
    # Convert and index ShapeIT4 Phased Data to vcf.gz format
    #------------------------------------------------------------------
    bcftools convert \
        --output ${WD}/0502/ABCD.chr${i}.phased.vcf.gz \
        --output-type z \
        --threads 12 \
        ${WD}/0502/ABCD.chr${i}.phased.vcf;

    bcftools index \
        --threads 12 \
        --output ${WD}/0502/ABCD.chr${i}.phased.vcf.gz.csi \
        ${WD}/0502/ABCD.chr${i}.phased.vcf.gz;

    #------------------------------------------------------------------
    # Create positioning and index file for reference
    #------------------------------------------------------------------
    bcftools convert \
        --output ${RF}/1000Gp3hg37/chr${i}.vcf.gz \
        --output-type z \
        --threads 12 \
        ${RF}/1000Gp3hg37/chr${i}.vcf.gz;

    bcftools view \
        --threads 12 \
        --drop-genotypes \
        --output-type z \
        --output ${RF}/1000Gp3hg37/chr${i}.forChunking.vcf.gz \
        ${RF}/1000Gp3hg37/chr${i}.vcf.gz;

    bcftools index \
        --threads 12 \
        --output ${RF}/1000Gp3hg37/chr${i}.forChunking.vcf.gz.csi \
        ${RF}/1000Gp3hg37/chr${i}.forChunking.vcf.gz;

    bcftools index \
        --threads 12 \
        --output ${RF}/1000Gp3hg37/chr${i}.vcf.gz.csi \
        ${RF}/1000Gp3hg37/chr${i}.vcf.gz;

    #------------------------------------------------------------------
    # Chunk the chromosome
    #------------------------------------------------------------------
    .${CHUNKER} \
        --h ${RF}/1000Gp3hg37/chr${i}.forChunking.vcf.gz \
        --g ${WD}/0502/ABCD.chr${i}.phased.vcf.gz \
        --r ${i} \
        --o ${WD}/0503/coordinates_chr${i}.txt;

    #------------------------------------------------------------------
    # Convert reference to imp5 format
    #------------------------------------------------------------------
    .${CONVERTER} \
        --h ${RF}/1000Gp3hg37/chr${i}.vcf.gz \
        --r ${i} \
        --threads 12 \
        --o ${RF}/1000Gp3hg37/chr${i}.imp5;

    #------------------------------------------------------------------
    # Impute chunks
    #------------------------------------------------------------------
    while read line; do
        ID=$(echo $line|awk '{print $1}')
        BUFFER=$(echo $line|awk '{print $3}')
        CHUNK=$(echo $line|awk '{print $4}')
        
        .${IMPUTER} \
            --h ${RF}/1000Gp3hg37/chr${i}.vcf.gz \
            --m ${RF}/SHAPEIT4MAPS/chr${i}.b37.gmap.gz \
            --g ${WD}/0502/ABCD.chr${i}.phased.vcf.gz \
            --r ${CHUNK} \
            --buffer-region ${BUFFER} \
            --threads 12 \
            --pbwt-depth 8 \
            --o ${WD}/0503/ABCD.chr${i}.chunk${ID}.imputed.vcf.gz;
    done < ${WD}/0503/coordinates_chr${i}.txt
done

echo "Finished."
```

## 05.04 Stitch and filter imputed data 
```bash
# Make a holding directory
mkdir -p ${WD}/0504

# Run the script
sbatch SLURM_Scripts/FilterImputed.sh
```

### FilterImputed.sh
```bash
#!/bin/bash
#SBATCH --job-name=FilterImputed              # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=njofrica@ufl.edu                # Where to send mail
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=12                           # Number of CPU cores per task
#SBATCH --mem=30gb                                   # Total memory limit
#SBATCH --time=7-00:00:00                             # Time limit hrs:min:sec
#SBATCH --output=.../LogFiles/SLURM_%j-%x.log   # Standard output and error log

date;hostname;pwd

export OMP_NUM_THREADS=12

#----------------------------------------------------------------------
# Load necessary modules
#----------------------------------------------------------------------
module load bcftools/1.13

#----------------------------------------------------------------------
# Assign paths to variable names go to WorkDir
#----------------------------------------------------------------------
PD="..."
WD="..."
SLURM="..."
LOGS="..."
RF="..."

mkdir -p ${WD}/0504 ${WD}/0504/TEMP
cd ${WD}

#----------------------------------------------------------------------
# Process Imputed Files
#----------------------------------------------------------------------
for chr in {1..22};
do
    #------------------------------------------------------------------
    # List files to stitch 
    #------------------------------------------------------------------
    ls ${WD}/0503/ABCD.chr${chr}.chunk* > \
    ${WD}/0504/ImputedChunks_chr${chr}.txt;

    #------------------------------------------------------------------
    # Concatonate files 
    #------------------------------------------------------------------
    bcftools concat \
        --naive \
        --file-list ${WD}/0504/ImputedChunks_chr${chr}.txt \
        --output-type z \
        --threads 12 \
        --output ${WD}/0504/ABCD.chr${chr}.imputed.vcf.gz;

    #------------------------------------------------------------------
    # Exclude low MAF and low INFO SNPs
    #------------------------------------------------------------------
    bcftools view \
        --include "INFO>0.8 & AF>0.01" \
        --output-type z \
        --threads 12 \
        --output ${WD}/0504/ABCD.chr${chr}.imputed.filtered.vcf.gz \
        ${WD}/0504/ABCD.chr${chr}.imputed.vcf.gz;

    #------------------------------------------------------------------
    # Sort files and output to bcf.gz
    #------------------------------------------------------------------
    bcftools sort \
        --output-type b \
        --output ${WD}/0504/ABCD.chr${chr}.imputed.filtered.sorted.bcf.gz \
        --temp-dir ${WD}/0504/TEMP \
        ${WD}/0504/ABCD.chr${chr}.imputed.filtered.vcf.gz;

    #------------------------------------------------------------------
    # Index files 
    #------------------------------------------------------------------
    bcftools index \
        ${WD}/0504/ABCD.chr${chr}.imputed.filtered.sorted.bcf.gz   

    #------------------------------------------------------------------
    # List files for archive 
    #------------------------------------------------------------------
    cat ${WD}/0504/ImputedChunks_chr${chr}.txt >> ${WD}/0504/ImputedChunksList.txt
    echo "${WD}/0504/ABCD.chr${chr}.imputed.vcf.gz" >> ${WD}/0504/ImputedStichedList.txt
    echo "${WD}/0504/ABCD.chr${chr}.imputed.filtered.vcf.gz" >> ${WD}/0504/ImputedFilteredList.txt;
done

#----------------------------------------------------------------------
# Archive and transfer to Orange
#----------------------------------------------------------------------
tar -cvf ${WD}/0504/ImputedChunks.tar -T ${WD}/0504/ImputedChunksList.txt
tar -cvf ${WD}/0504/ImputedStiched.tar -T ${WD}/0504/ImputedStichedList.txt
tar -cvf ${WD}/0504/ImputedFiltered.tar -T ${WD}/0504/ImputedFilteredList.txt

mv ${WD}/0504/ImputedChunks.tar /orange/carolmathews/njofrica/ABCD/ImputedChunks.tar
mv ${WD}/0504/ImputedStiched.tar /orange/carolmathews/njofrica/ABCD/ImputedStiched.tar
mv ${WD}/0504/ImputedFiltered.tar /orange/carolmathews/njofrica/ABCD/ImputedFiltered.tar

#----------------------------------------------------------------------
# Delete excess files
#----------------------------------------------------------------------
rm ${WD}/0504/ImputedChunksList.txt
rm ${WD}/0504/ImputedStichedList.txt
rm ${WD}/0504/ImputedFilteredList.txt
rm ${WD}/0504/ImputedChunks_chr*
rm ${WD}/0504/*imputed.filtered.vcf.gz
rm ${WD}/0504/*.imputed.vcf.gz
rm ${WD}/0503/*.imputed.vcf.gz
```

### Preprocess data
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
cd ${WD}/CSI

awk 'NR == 1 || $9 < 0.005' ${WD}/1002/ocd_aug2017 | cut -f2 | awk 'NR>1' >  ${WD}/CSI/OCD_SigHit

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
    --include "INFO>0.9 & AF>=0.005 & AF<=0.995" \
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
```
