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

abcd_abcls01 <- read_delim("./ABCD_Data/abcd_abcls01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE)[-1,]

abcd_bpmt01 <- read_delim("./ABCD_Data/abcd_bpmt01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE)[-1,]

cct01 <- read_delim("./ABCD_Data/cct01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE)[-1,]

abcd_cna01 <- read_delim("./ABCD_Data/abcd_cna01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999"))[-1,]

crpbi01 <- read_delim("./ABCD_Data/crpbi01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888"))[-1,]

abcd_cb01 <- read_delim("./ABCD_Data/abcd_cb01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

dhx01 <- read_delim("./ABCD_Data/dhx01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

abcd_eatqp01 <- read_delim("./ABCD_Data/abcd_eatqp01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

fhxp102 <- read_delim("./ABCD_Data/fhxp102.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

fhxp201 <- read_delim("./ABCD_Data/fhxp201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

abcd_gdss01 <- read_delim("./ABCD_Data/abcd_gdss01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

abcd_hsss01 <- read_delim("./ABCD_Data/abcd_hsss01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

lmtp201 <- read_delim("./ABCD_Data/lmtp201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

abcd_lpds01 <- read_delim("./ABCD_Data/abcd_lpds01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777"))[-1,]

abcd_lssmh01 <- read_delim("./ABCD_Data/abcd_lssmh01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA", "999", "888", "777", "-1"))[-1,]

abcd_lsssa01 <- read_delim("./ABCD_Data/abcd_lsssa01.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip101 <- read_delim("./ABCD_Data/abcd_drsip101.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip201 <- read_delim("./ABCD_Data/abcd_drsip201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip301 <- read_delim("./ABCD_Data/abcd_drsip301.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip401 <- read_delim("./ABCD_Data/abcd_drsip401.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip501 <- read_delim("./ABCD_Data/abcd_drsip501.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip601 <- read_delim("./ABCD_Data/abcd_drsip601.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_drsip701 <- read_delim("./ABCD_Data/abcd_drsip701.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_mrisdp10201 <- read_delim("./ABCD_Data/abcd_mrisdp10201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_mrisdp20201 <- read_delim("./ABCD_Data/abcd_mrisdp20201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_mrisdp30201 <- read_delim("./ABCD_Data/abcd_mrisdp30201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_smrip10201 <- read_delim("./ABCD_Data/abcd_smrip10201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_smrip20201 <- read_delim("./ABCD_Data/abcd_smrip20201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_smrip30201 <- read_delim("./ABCD_Data/abcd_smrip30201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_ddtidp101 <- read_delim("./ABCD_Data/abcd_ddtidp101.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_ddtidp201 <- read_delim("./ABCD_Data/abcd_ddtidp201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_ddtifp101 <- read_delim("./ABCD_Data/abcd_ddtifp101.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_ddtifp201 <- read_delim("./ABCD_Data/abcd_ddtifp201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_dmdtifp101 <- read_delim("./ABCD_Data/abcd_dmdtifp101.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_dmdtifp201 <- read_delim("./ABCD_Data/abcd_dmdtifp201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_dti_p101 <- read_delim("./ABCD_Data/abcd_dti_p101.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]

abcd_dti_p201 <- read_delim("./ABCD_Data/abcd_dti_p201.txt", delim = "\t", escape_double = FALSE, col_types = "c", trim_ws = TRUE, na = c("", "NA"))[-1,]


# RxVars
RxALL <- "Concerta|Citalopram|Zoloft|Zofran|Amoxicillin|Azithromycin|Prozac|Buspirone|Diphenhydramine|Vitamin|Synthroid|levothyroxine|Lexapro|Celexa|Bupropion|Progesterone|Wellbutrin|Folic Acid|Insulin|Metformin|Estrogens|Clomiphene|Advair|Albuterol|Claritin-D|Heparin|Abilify|Sertraline|Lipitor|Glipizide|Methyldopa|Amlodipine|Ambien|Tylenol|Percocet|Flovent|Levothyroxine|Thyroxine|Norvasc|Symbicort|Metoprolol|Terbutaline|Acetaminophen|Zyrtec|Fluticasone|Hydrocodone|Effexor|Valtrex|Truvada|Lamictal|Paxil|Tranylcypromine|Atarax|Rocaltrol|Phenergan|Humalog|Flexeril|NovoLog|Methylprednisolone|Azathioprine|Topamax|Fluoxetine|PYRIDOXINE|Antibiotic|Cetirizine|Singulair|Propylthiouracil|Prilosec|Sumatriptan|ferrous sulfate|Lisinopril|Nausea Control|Oxycontin|PreNexa|Nexium|Levo-T|insulin|Aspirin|Lortab|Omeprazole|Levoxyl|Labetalol|Keflex|ProAir|Acyclovir|thyroid|Lovenox|Ventolin|Xanax|Lithium|Morphine|Plaquenil|Methadone|Claritin|Remicade|Glucophage|Protonix|Paroxetine|Vicodin|Promethazine|Nifedipine|Lantus|Bactrim|Prednisone|Famotidine|Escitalopram|Carbamazepine|Rhinocort|Pentasa|Glucagon|Cymbalta|oxcarbazepine|Flagyl|Hydrochlorothiazide|ORTHO-NOVUM|Allegra|Sulfasalazine|Imitrex|Prevacid|Ranitidine|Maxalt|Colace|Nitrofurantoin|Thyroid|Zantac|Levothroid|Oxycodone|Divalproex|Calcium|Iron|quetiapine|Hyoscyamine|Tums|Lasix|Sotalol|Imuran|Golytely|Keppra|lansoprazole|Atenolol|Etanercept|1151133 Pill|284743 284743|3407 3407|861206 861206|11103 Vaginal Cream|5093 5093|865575 865575|10369 10369|6915 6915|Xyzal|Reglan|Ondansetron|Losartan|Asthmahaler|Metoclopramide|Proventil|Lyrica|218090 218090|YASMIN|Pepcid|Fish Oils|Trazodone|Seroquel|Bisoprolol|152699 152699|Clindamycin|Dovonex|Trileptal|Zoladex|Penicillin|Adderall|Fioricet|Norco|T125 Hormone|butalbital|Prometrium|Glyburide|Methimazole|prednisolone|Cephalexin|Z-PAK|224913 224913|203423 203423|Sudafed|Suboxone|Warfarin|Prezista|Magnesium Sulfate|Levonorgestrel|mesalamine|Nasonex|Propranolol|Lialda|7417 7417|Diovan|tolterodine|Vistaril|Phenobarbital|pantoprazole|Sulfate Ion|Clozapine|Ibuprofen|Depakote|Nasal Inhalant|Midrin|glatiramer|Cytomel|Viramune|Pulmicor|Vancomycin|Benadryl|heparin|cyclobenzaprine|Neurontin|Dexamethasone|Valium|Diflucan|Levemir|Kaletra|Zovirax|153668 153668|Soma|Detrol|Flonase|Nifedical|Fentanyl|Aciphex|Celebrex|1151133.0 Pill|Diabeta|montelukast|204527 tetanus|Folgard|Unisom|Ultram|Folate|Loratadine|7417.0 7417|1649574.0 Injection|24941.0 24941|4821.0 4821|Furosemide|Folplex|Augmentin|Omega-3|Lorazepam|Betamethasone|1649574 Injection|Zithromax|1536581 1536581|Gonadorelin|Humira|Follistim|Advil|Clomicalm|Contrave|Mirena|Klonopin|ORTHO-CEPT|Gonadotropin|Atripla|Follistim|ORTHO TRI-CYCLEN|Glucovance|Benicar|Tamoxifen|telithromycin|Enalapril|Baclofen|Codeine|Clobetasol|Adipex-P|NuvaRing|Estradiol|Evra|Doxycycline|Rifampin|LO/OVRAL|Tamiflu|Doxepin|Temazepam|Cyproheptadine|Nortriptyline|5093.0 5093|Antihistamine|Ziac|Ativan|rifapentine|Gonal|Follicle|Vyvanse|ORTHO-CYCLEN|202421 202421|Portia|predniSONE|Levora|Xulane|Ketorolac|Tetracycline|ESTROSTEP|Amitriptyline|Femara|Retin-A|Norvir|Naproxen|methotrexate|LaMICtal|Natelle|Ibuprohm|1310235 1310235|Bromocriptine|Copaxone|Lotrel|metoprolol|Travatan|Mircette|ezetimibe|Asacol|Sulfamethoxazole|Methylphenidate|Geodon|lamotrigine|Calcitriol|sennosides|Depo-Provera|Tri-Sprintec|BC Headache|Toprol|Colazal|Lodine|Cortisone Pill|Clonazepam|Lo Loestrin|1-Propanol|1-propoxy-2-propanol|gabapentin|seasonique|Oxymorphone|Roxicet|Lupron|tramadol|Beclomethasone|227528.0 227528|Gardasil|227528.0 227528|Compazine|Lidocaine|cefprozil|Natafort"

RxPsych <- "Concerta|Citalopram|Zoloft|Prozac|Buspirone|Lexapro|Celexa|Bupropion|Wellbutrin|Abilify|Sertraline|Ambien|Effexor|Paxil|Tranylcypromine|Fluoxetine|Xanax|Lithium|Paroxetine|Escitalopram|Cymbalta|Maxalt|quetiapine|Trazodone|Seroquel|Adderall|Phenobarbital|Clozapine|Depakote|Valium|Lorazepam|Clomicalm|Contrave|Klonopin|Adipex-P|Doxepin|Temazepam|Nortriptyline|Ativan|Vyvanse|Amitriptyline|Methylphenidate|Geodon|Clonazepam|Compazine"

RxGI <- "Zofran|Phenergan|Prilosec|Nausea Control|Nexium|Omeprazole|Protonix|Promethazine|Famotidine|Prevacid|Ranitidine|Colace|Zantac|Hyoscyamine|Tums|Golytely|lansoprazole|Reglan|Ondansetron|Metoclopramide|Pepcid|mesalamine|Lialda|Vistaril|pantoprazole|Pulmicort|Aciphex|sennosides|Asacol|Colazal|Compazine"

RxAntibiotic <- "Amoxicillin|Azithromycin|Antibiotic|Keflex|Bactrim|Flagyl|Nitrofurantoin|Clindamycin|Penicillin|Cephalexin|Z-PAK|Vancomycin|Augmentin|Zithromax|telithromycin|Doxycycline|Rifampin|rifapentine|Tetracycline|Sulfamethoxazole|cefprozil"

RxAntihistamine <- "Diphenhydramine|Claritin-D|Zyrtec|Fluticasone|Atarax|Phenergan|Methylprednisolone|Cetirizine|Claritin|Promethazine|Prednisone|Famotidine|Rhinocort|Allegra|Ranitidine|Zantac|Xyzal|Pepcid|prednisolone|Sudafed|Nasonex|Vistaril|Benadryl|montelukast|Unisom|Loratadine|Betamethasone|Cyproheptadine|Antihistamine|predniSONE|Cortisone Pill"

RxSupplments <- "Vitamin|Folic Acid|Rocaltrol|PYRIDOXINE|ferrous sulfate|PreNexa|Calcium|Iron|Fish Oils|Dovonex|Magnesium Sulfate|Sulfate Ion|Folgard|Folate|Folplex|Omega-3|Natelle|Calcitriol|Natafort"

RxThyroid <- "Synthroid|levothyroxine|Levothyroxine|Thyroxine|Propylthiouracil|Levo-T|Levoxyl|thyroid|Thyroid|Levothroid|T125 Hormone|Methimazole|Cytomel"

RxSexHor <- "Progesterone|Estrogens|Clomiphene|ORTHO-NOVUM|YASMIN|Zoladex|Prometrium|Levonorgestrel|Gonadorelin|Follistim|ORTHO-CEPT|Gonadotropin|Follistim|ORTHO TRI-CYCLEN|Tamoxifen|NuvaRing|Estradiol|Evra|LO/OVRAL|Gonal|Follicle|ORTHO-CYCLEN|Portia|Levora|Xulane|ESTROSTEP|Femara|Bromocriptine|Mircette|Depo-Provera|Tri-Sprintec|Lo Loestrin|seasonique|Lupron"

RxDiab <- "Insulin|Metformin|Glipizide|Humalog|NovoLog|insulin|Glucophage|Lantus|Glyburide|Levemir|Diabeta|Glucovance"

RxAsthma <- "Advair|Albuterol|Flovent|Symbicort|Terbutaline|Singulair|ProAir|Ventolin|Asthmahaler|Proventil|Pulmicort|Nasal Inhalant|Flonase|montelukast|Cortisone Pill|Beclomethasone"

RxAnticoag <- "Heparin|Lovenox|Warfarin|heparin"

RxHypertens <- "Lipitor|Methyldopa|Amlodipine|Norvasc|Metoprolol|Lisinopril|Labetalol|Nifedipine|Glucagon|Hydrochlorothiazide|Lasix|Sotalol|Atenolol|Losartan|Bisoprolol|Propranolol|Diovan|tolterodine|Detrol|Nifedical|Furosemide|Benicar|Enalapril|Ziac|Lotrel|metoprolol|ezetimibe|Toprol"

RxPain <- "Tylenol|Percocet|Acetaminophen|Hydrocodone|Flexeril|Sumatriptan|Oxycontin|Aspirin|Lortab|Morphine|Methadone|Vicodin|Pentasa|Imitrex|Oxycodone|Lyrica|Fioricet|Norco|butalbital|Suboxone|Ibuprofen|Midrin|cyclobenzaprine|Dexamethasone|Soma|Fentanyl|Celebrex|Ultram|Advil|Baclofen|Codeine|Ketorolac|Naproxen|Ibuprohm|BC Headache|Lodine|Oxymorphone|Roxicet|tramadol|Lidocaine"

RxAntivir <- "Valtrex|Truvada|Acyclovir|Prezista|Viramune|Kaletra|Zovirax|Atripla|Tamiflu|Norvir"

RxAnticonv <- "Lamictal|Topamax|Carbamazepine|oxcarbazepine|Divalproex|Keppra|Lyrica|Trileptal|Phenobarbital|Depakote|Neurontin|Lorazepam|Klonopin|Ativan|LaMICtal|lamotrigine|Clonazepam|gabapentin"

RxImmunosup <- "Azathioprine|Plaquenil|Remicade|Imuran|Etanercept|glatiramer|Humira|methotrexate|Copaxone"

RxAntirheum <- "Sulfasalazine|Etanercept"

RxAntifung <- "Diflucan"

# Filter datasets 
abcd_abcls01 <- abcd_abcls01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, ends_with("_r"))) %>%
    mutate_at(., vars(contains("abcl")), as.numeric) 

abcd_bpmt01 <- abcd_bpmt01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("bpmt_q"))) %>%
    mutate_at(., vars(contains("bpmt_q")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("bpmt_q")), mean, na.rm=T)

cct01 <- cct01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, cash_choice_task)) %>%
    mutate_at(., vars(contains("cash_choice_task")), as.numeric) 

abcd_cna01 <- abcd_cna01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, ends_with("_p"))) %>%
    mutate_at(., vars(contains("cna_")), as.numeric) 

crpbi01 <- crpbi01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, ends_with("_y"))) %>%
    select(-crpbi_caregiver1_y) %>%
    mutate_at(., vars(contains("_y")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("_y")), mean, na.rm=T)

abcd_cb01 <- abcd_cb01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("cybb_phenx"))) %>%
    mutate_at(., vars(contains("cybb_phenx")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("cybb_phenx")), mean, na.rm=T)

dhx01 <- dhx01 %>%
    filter(visit %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    rowwise() %>%
    mutate(BirthWeight=as.numeric(birth_weight_lbs)*0.453592+as.numeric(birth_weight_lbs)*0.0283495,
        MomAgeBirth=as.numeric(devhx_3_p),
        DadAgeBirth=as.numeric(devhx_4_p),
        WasTwin=as.numeric(devhx_5_p),
        WasPlanned=as.numeric(devhx_6_p),
        FoundPregWeek=as.numeric(devhx_7_p),
        BeforeKnowPregMedAny=as.numeric(ifelse(
            grepl(RxALL, devhx_8_rxnorm_med1) | 
            grepl(RxALL, devhx_8_rxnorm_med2) | 
            grepl(RxALL, devhx_8_rxnorm_med3), 1, 0)),
        BeforePsychotropic=as.numeric(ifelse(
            grepl(RxPsych, devhx_8_rxnorm_med1) | 
            grepl(RxPsych, devhx_8_rxnorm_med2) | 
            grepl(RxPsych, devhx_8_rxnorm_med3), 1, 0)),
        BeforeGI=as.numeric(ifelse(
            grepl(RxGI, devhx_8_rxnorm_med1) | 
            grepl(RxGI, devhx_8_rxnorm_med2) | 
            grepl(RxGI, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAntibiotic=as.numeric(ifelse(
            grepl(RxAntibiotic, devhx_8_rxnorm_med1) | 
            grepl(RxAntibiotic, devhx_8_rxnorm_med2) | 
            grepl(RxAntibiotic, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAntihistamine=as.numeric(ifelse(
            grepl(RxAntihistamine, devhx_8_rxnorm_med1) | 
            grepl(RxAntihistamine, devhx_8_rxnorm_med2) | 
            grepl(RxAntihistamine, devhx_8_rxnorm_med3), 1, 0)),
        BeforeSupplement=as.numeric(ifelse(
            grepl(RxSupplments, devhx_8_rxnorm_med1) | 
            grepl(RxSupplments, devhx_8_rxnorm_med2) | 
            grepl(RxSupplments, devhx_8_rxnorm_med3), 1, 0)),
        BeforeThyroid=as.numeric(ifelse(
            grepl(RxThyroid, devhx_8_rxnorm_med1) | 
            grepl(RxThyroid, devhx_8_rxnorm_med2) | 
            grepl(RxThyroid, devhx_8_rxnorm_med3), 1, 0)),
        BeforeSexHormone=as.numeric(ifelse(
            grepl(RxSexHor, devhx_8_rxnorm_med1) | 
            grepl(RxSexHor, devhx_8_rxnorm_med2) | 
            grepl(RxSexHor, devhx_8_rxnorm_med3), 1, 0)),
        BeforeDiabetic=as.numeric(ifelse(
            grepl(RxDiab, devhx_8_rxnorm_med1) | 
            grepl(RxDiab, devhx_8_rxnorm_med2) | 
            grepl(RxDiab, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAsthma=as.numeric(ifelse(
            grepl(RxAsthma, devhx_8_rxnorm_med1) | 
            grepl(RxAsthma, devhx_8_rxnorm_med2) | 
            grepl(RxAsthma, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAnticoagulant=as.numeric(ifelse(
            grepl(RxAnticoag, devhx_8_rxnorm_med1) | 
            grepl(RxAnticoag, devhx_8_rxnorm_med2) | 
            grepl(RxAnticoag, devhx_8_rxnorm_med3), 1, 0)),
        BeforeHypertensive=as.numeric(ifelse(
            grepl(RxHypertens, devhx_8_rxnorm_med1) | 
            grepl(RxHypertens, devhx_8_rxnorm_med2) | 
            grepl(RxHypertens, devhx_8_rxnorm_med3), 1, 0)),
        BeforePain=as.numeric(ifelse(
            grepl(RxPain, devhx_8_rxnorm_med1) | 
            grepl(RxPain, devhx_8_rxnorm_med2) | 
            grepl(RxPain, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAntiviral=as.numeric(ifelse(
            grepl(RxAntivir, devhx_8_rxnorm_med1) | 
            grepl(RxAntivir, devhx_8_rxnorm_med2) | 
            grepl(RxAntivir, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAnticonvulsant=as.numeric(ifelse(
            grepl(RxAnticonv, devhx_8_rxnorm_med1) | 
            grepl(RxAnticonv, devhx_8_rxnorm_med2) | 
            grepl(RxAnticonv, devhx_8_rxnorm_med3), 1, 0)),
        BeforeImmunosuppressant=as.numeric(ifelse(
            grepl(RxImmunosup, devhx_8_rxnorm_med1) | 
            grepl(RxImmunosup, devhx_8_rxnorm_med2) | 
            grepl(RxImmunosup, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAntirheumatic=as.numeric(ifelse(
            grepl(RxAntirheum, devhx_8_rxnorm_med1) | 
            grepl(RxAntirheum, devhx_8_rxnorm_med2) | 
            grepl(RxAntirheum, devhx_8_rxnorm_med3), 1, 0)),
        BeforeAntifungal=as.numeric(ifelse(
            grepl(RxAntifung, devhx_8_rxnorm_med1) | 
            grepl(RxAntifung, devhx_8_rxnorm_med2) | 
            grepl(RxAntifung, devhx_8_rxnorm_med3), 1, 0)),
        AfterKnowPregMedAny=as.numeric(ifelse(
            grepl(RxALL, devhx_9_med1_rxnorm) | 
            grepl(RxALL, devhx_9_med2_rxnorm) | 
            grepl(RxALL, devhx_9_med3_rxnorm) | 
            grepl(RxALL, devhx_9_med4_rxnorm) | 
            grepl(RxALL, devhx_9_med5_rxnorm), 1, 0)),
        AfterPsychotropic=as.numeric(ifelse(
            grepl(RxPsych, devhx_9_med1_rxnorm) | 
            grepl(RxPsych, devhx_9_med2_rxnorm) | 
            grepl(RxPsych, devhx_9_med3_rxnorm) | 
            grepl(RxPsych, devhx_9_med4_rxnorm) | 
            grepl(RxPsych, devhx_9_med5_rxnorm), 1, 0)),
        AfterGI=as.numeric(ifelse(
            grepl(RxGI, devhx_9_med1_rxnorm) | 
            grepl(RxGI, devhx_9_med2_rxnorm) | 
            grepl(RxGI, devhx_9_med3_rxnorm) | 
            grepl(RxGI, devhx_9_med4_rxnorm) | 
            grepl(RxGI, devhx_9_med5_rxnorm), 1, 0)),
        AfterAntibiotic=as.numeric(ifelse(
            grepl(RxAntibiotic, devhx_9_med1_rxnorm) | 
            grepl(RxAntibiotic, devhx_9_med2_rxnorm) | 
            grepl(RxAntibiotic, devhx_9_med3_rxnorm) | 
            grepl(RxAntibiotic, devhx_9_med4_rxnorm) | 
            grepl(RxAntibiotic, devhx_9_med5_rxnorm), 1, 0)),
        AfterAntihistamine=as.numeric(ifelse(
            grepl(RxAntihistamine, devhx_9_med1_rxnorm) | 
            grepl(RxAntihistamine, devhx_9_med2_rxnorm) | 
            grepl(RxAntihistamine, devhx_9_med3_rxnorm) | 
            grepl(RxAntihistamine, devhx_9_med4_rxnorm) | 
            grepl(RxAntihistamine, devhx_9_med5_rxnorm), 1, 0)),
        AfterSupplement=as.numeric(ifelse(
            grepl(RxSupplments, devhx_9_med1_rxnorm) | 
            grepl(RxSupplments, devhx_9_med2_rxnorm) | 
            grepl(RxSupplments, devhx_9_med3_rxnorm) | 
            grepl(RxSupplments, devhx_9_med4_rxnorm) | 
            grepl(RxSupplments, devhx_9_med5_rxnorm), 1, 0)),
        AfterThyroid=as.numeric(ifelse(
            grepl(RxThyroid, devhx_9_med1_rxnorm) | 
            grepl(RxThyroid, devhx_9_med2_rxnorm) | 
            grepl(RxThyroid, devhx_9_med3_rxnorm) | 
            grepl(RxThyroid, devhx_9_med4_rxnorm) | 
            grepl(RxThyroid, devhx_9_med5_rxnorm), 1, 0)),
        AfterSexHormone=as.numeric(ifelse(
            grepl(RxSexHor, devhx_9_med1_rxnorm) | 
            grepl(RxSexHor, devhx_9_med2_rxnorm) | 
            grepl(RxSexHor, devhx_9_med3_rxnorm) | 
            grepl(RxSexHor, devhx_9_med4_rxnorm) | 
            grepl(RxSexHor, devhx_9_med5_rxnorm), 1, 0)),
        AfterDiabetic=as.numeric(ifelse(
            grepl(RxDiab, devhx_9_med1_rxnorm) | 
            grepl(RxDiab, devhx_9_med2_rxnorm) | 
            grepl(RxDiab, devhx_9_med3_rxnorm) | 
            grepl(RxDiab, devhx_9_med4_rxnorm) | 
            grepl(RxDiab, devhx_9_med5_rxnorm), 1, 0)),
        AfterAsthma=as.numeric(ifelse(
            grepl(RxAsthma, devhx_9_med1_rxnorm) | 
            grepl(RxAsthma, devhx_9_med2_rxnorm) | 
            grepl(RxAsthma, devhx_9_med3_rxnorm) | 
            grepl(RxAsthma, devhx_9_med4_rxnorm) | 
            grepl(RxAsthma, devhx_9_med5_rxnorm), 1, 0)),
        AfterAnticoagulant=as.numeric(ifelse(
            grepl(RxAnticoag, devhx_9_med1_rxnorm) | 
            grepl(RxAnticoag, devhx_9_med2_rxnorm) | 
            grepl(RxAnticoag, devhx_9_med3_rxnorm) | 
            grepl(RxAnticoag, devhx_9_med4_rxnorm) | 
            grepl(RxAnticoag, devhx_9_med5_rxnorm), 1, 0)),
        AfterHypertensive=as.numeric(ifelse(
            grepl(RxHypertens, devhx_9_med1_rxnorm) | 
            grepl(RxHypertens, devhx_9_med2_rxnorm) | 
            grepl(RxHypertens, devhx_9_med3_rxnorm) | 
            grepl(RxHypertens, devhx_9_med4_rxnorm) | 
            grepl(RxHypertens, devhx_9_med5_rxnorm), 1, 0)),
        AfterPain=as.numeric(ifelse(
            grepl(RxPain, devhx_9_med1_rxnorm) | 
            grepl(RxPain, devhx_9_med2_rxnorm) | 
            grepl(RxPain, devhx_9_med3_rxnorm) | 
            grepl(RxPain, devhx_9_med4_rxnorm) | 
            grepl(RxPain, devhx_9_med5_rxnorm), 1, 0)),
        AfterAntiviral=as.numeric(ifelse(
            grepl(RxAntivir, devhx_9_med1_rxnorm) | 
            grepl(RxAntivir, devhx_9_med2_rxnorm) | 
            grepl(RxAntivir, devhx_9_med3_rxnorm) | 
            grepl(RxAntivir, devhx_9_med4_rxnorm) | 
            grepl(RxAntivir, devhx_9_med5_rxnorm), 1, 0)),
        AfterAnticonvulsant=as.numeric(ifelse(
            grepl(RxAnticonv, devhx_9_med1_rxnorm) | 
            grepl(RxAnticonv, devhx_9_med2_rxnorm) | 
            grepl(RxAnticonv, devhx_9_med3_rxnorm) | 
            grepl(RxAnticonv, devhx_9_med4_rxnorm) | 
            grepl(RxAnticonv, devhx_9_med5_rxnorm), 1, 0)),
        AfterImmunosuppressant=as.numeric(ifelse(
            grepl(RxImmunosup, devhx_9_med1_rxnorm) | 
            grepl(RxImmunosup, devhx_9_med2_rxnorm) | 
            grepl(RxImmunosup, devhx_9_med3_rxnorm) | 
            grepl(RxImmunosup, devhx_9_med4_rxnorm) | 
            grepl(RxImmunosup, devhx_9_med5_rxnorm), 1, 0)),
        AfterAntirheumatic=as.numeric(ifelse(
            grepl(RxAntirheum, devhx_9_med1_rxnorm) | 
            grepl(RxAntirheum, devhx_9_med2_rxnorm) | 
            grepl(RxAntirheum, devhx_9_med3_rxnorm) | 
            grepl(RxAntirheum, devhx_9_med4_rxnorm) | 
            grepl(RxAntirheum, devhx_9_med5_rxnorm), 1, 0)),
        AfterAntifungal=as.numeric(ifelse(
            grepl(RxAntifung, devhx_9_med1_rxnorm) | 
            grepl(RxAntifung, devhx_9_med2_rxnorm) | 
            grepl(RxAntifung, devhx_9_med3_rxnorm) | 
            grepl(RxAntifung, devhx_9_med4_rxnorm) | 
            grepl(RxAntifung, devhx_9_med5_rxnorm), 1, 0)),
        BeforeCigsDay=as.numeric(devhx_8_cigs_per_day),
        BeforeAlcoholAvgWeek=as.numeric(devhx_8_alchohol_avg),
        BeforeMarijuanaDay=as.numeric(devhx_8_marijuana_amt),
        BeforeCrackCocaineDat=as.numeric(devhx_8_coc_crack_amt),
        BeforeHeroineMorphine=as.numeric(devhx_8_her_morph_amt),
        BeforeOxycontin=as.numeric(devhx_8_oxycont_amt),
        AfterCigsDay=as.numeric(devhx_9_cigs_per_day),
        AfterAlcoholAvgWeek=as.numeric(devhx_9_alchohol_avg),
        AfterMarijuanaDay=as.numeric(devhx_9_marijuana_amt),
        AfterCrackCocaineDat=as.numeric(devhx_9_coc_crack_amt),
        AfterHeroineMorphine=as.numeric(devhx_9_her_morph_amt),
        AfterOxycontin=as.numeric(devhx_9_oxycont_amt),
        Caffeine=as.numeric(devhx_caffeine_amt),
        Nausea=as.numeric(devhx_10a3_p),
        Bleeding=as.numeric(devhx_10b3_p),
        Eclampsia=as.numeric(devhx_10c3_p),
        GallBladder=as.numeric(devhx_10d3_p),
        Proteinuria=as.numeric(devhx_10e3_p),
        Rubella=as.numeric(devhx_10f3_p),
        Anemia=as.numeric(devhx_10g3_p),
        UTI=as.numeric(devhx_10h3_p),
        Diabetes=as.numeric(devhx_10i3_p),
        Hypertension=as.numeric(devhx_10j3_p),
        Previa=as.numeric(devhx_10k3_p),
        Injury=as.numeric(devhx_10l3_p),
        OtherConditioons=as.numeric(devhx_10m3_p),
        PrematureWeeks=as.numeric(devhx_12_p),
        Ceasarian=as.numeric(devhx_13_3_p),
        ComplicationsBlue=as.numeric(devhx_14a3_p),
        ComplicationsHeartbeat=as.numeric(devhx_14b3_p),
        ComplicationsBreathe=as.numeric(devhx_14c3_p),
        ComplicationsConvulsionas=as.numeric(devhx_14d3_p),
        ComplicationsJaundice=as.numeric(devhx_14e3_p),
        ComplicationsOxygen=as.numeric(devhx_14f3_p),
        ComplicationsTransfusion=as.numeric(devhx_14g3_p),
        ComplicationsRh=as.numeric(devhx_14h3_p),
        IncubatorDays=as.numeric(devhx_15),
        FeverDays=as.numeric(devhx_16_p),
        InfectionDays=as.numeric(devhx_17_p),
        BreastDays=as.numeric(devhx_18_p),
        RolloverAge=as.numeric(devhx_19a_p),
        SitAge=as.numeric(devhx_19b_p),
        WalkAge=as.numeric(devhx_19c_p),
        WordAge=as.numeric(devhx_19d_p),
        MotorDevlopment=as.numeric(devhx_20_p),
        SpeechDevelopmeent=as.numeric(devhx_21_p),
        BedWetting=as.numeric(devhx_22_3_p),
        ContinenceAge=as.numeric(devhx_23b_p),
        .keep="unused") %>%
    select(-c(visit, sex, dataset_id, dhx01_id, starts_with(c("devhx_", "accult_", "src_", "interview_", "birth_", "collection_"))))

abcd_eatqp01 <- abcd_eatqp01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, ends_with("_p"))) %>%
    mutate_at(., vars(ends_with("_p")), as.numeric)

fhxp102 <- fhxp102 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(-c(contains(c("_dk_", "_999", "interview_", "_id", "_title")), sex, eventname, famhx_select_language___1, famhx_1)) %>%
    mutate_at(., vars(contains("_")), as.numeric)
fhxp102$famhx_4_p[fhxp102$famhx_4_p == 7] <- NA
fhxp102$fam_history_5_yes_no[fhxp102$fam_history_5_yes_no == 7] <- NA
fhxp102$fam_history_6_yes_no[fhxp102$fam_history_6_yes_no == 7] <- NA

fhxp201 <- fhxp201 %>%
    filter(visit %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(-c(contains(c("_dk_", "_999", "interview_", "_id", "_title")), sex, visit)) %>%
    mutate_at(., vars(contains("_")), as.numeric)
fhxp201$fam_history_7_yes_no[fhxp201$fam_history_7_yes_no == 7] <- NA
fhxp201$fam_history_8_yes_no[fhxp201$fam_history_8_yes_no == 7] <- NA
fhxp201$fam_history_9_yes_no[fhxp201$fam_history_9_yes_no == 7] <- NA
fhxp201$fam_history_10_yes_no[fhxp201$fam_history_10_yes_no == 7] <- NA
fhxp201$fam_history_11_yes_no[fhxp201$fam_history_11_yes_no == 7] <- NA
fhxp201$fam_history_12_yes_no[fhxp201$fam_history_12_yes_no == 7] <- NA
fhxp201$fam_history_13_yes_no[fhxp201$fam_history_13_yes_no == 7] <- NA

abcd_gdss01 <- abcd_gdss01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(contains(c("gdt_scr_values_", "gdt_scr_parameters_")), subjectkey, gdt_scr_script_elapsedtime)) %>%
    mutate_at(., vars(contains("gdt_")), as.numeric)

abcd_hsss01 <- abcd_hsss01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, hormone_scr_dhea_mean, hormone_scr_hse_mean, hormone_scr_ert_mean)) %>%
    mutate_at(., vars(contains("_mean")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("_mean")), mean, na.rm=T)

lmtp201 <- lmtp201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("lmt_scr_"))) %>%
    mutate_at(., vars(contains("lmt_scr_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("lmt_scr_")), mean, na.rm=T)

abcd_lpds01 <- abcd_lpds01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("demo_"))) %>%
    mutate(demog_biomother = as.numeric(ifelse(demo_prim_l=="1", 1, 0)),
    demog_biofather = as.numeric(ifelse(demo_prim_l=="2", 1, 0)),
    demog_adoptive = as.numeric(ifelse(demo_prim_l=="3", 1, 0)),
    demog_custodial = as.numeric(ifelse(demo_prim_l=="4", 1, 0)),
    demog_cisgender = as.numeric(ifelse(demo_gender_id_v2_l %in% c("1", "2"), 1, 0)),
    demog_transgender = as.numeric(ifelse(demo_gender_id_v2_l %in% c("3", "4", "5"), 1, 0)),
    demog_nateng = as.numeric(ifelse(demo_nat_lang_l == "58", 1, 0)),
    demog_englprim = as.numeric(demo_nat_lang_2_l),
    demog_duallang = as.numeric(demo_dual_lang_v2_l),
    demog_duallangyr = as.numeric(demo_dual_lang_years_p___1) + as.numeric(demo_dual_lang_years_p___2) + as.numeric(demo_dual_lang_years_p___1) + as.numeric(demo_dual_lang_years_p___1) + as.numeric(demo_dual_lang_years_p___3) + as.numeric(demo_dual_lang_years_p___4) + as.numeric(demo_dual_lang_years_p___5) + as.numeric(demo_dual_lang_years_p___6) + as.numeric(demo_dual_lang_years_p___7) + as.numeric(demo_dual_lang_years_p___8) + as.numeric(demo_dual_lang_years_p___9) + as.numeric(demo_dual_lang_years_p___10),
    demog_christian = as.numeric(ifelse(demo_relig_v2_l %in% c("1", "2", "3", "4", "5", "6", "7", "11", "12", "13"), 1, 0)),
    demog_parage = as.numeric(demo_prnt_age_v2_l), 
    demog_parcis = as.numeric(ifelse(demo_prnt_gender_id_v2_l %in% c("1", "2"), 1, 0)),
    demog_parwhite = as.numeric(demo_prnt_race_a_v2_l___10),
    demog_parblack = as.numeric(demo_prnt_race_a_v2_l___11),
    demog_parnatam = as.numeric(as.numeric(demo_prnt_race_a_v2_l___12) | as.numeric(demo_prnt_race_a_v2_l___13) | as.numeric(demo_prnt_race_a_v2_l___14)),
    demog_paras = as.numeric(as.numeric(demo_prnt_race_a_v2_l___18) | as.numeric(demo_prnt_race_a_v2_l___19) | as.numeric(demo_prnt_race_a_v2_l___20) | as.numeric(demo_prnt_race_a_v2_l___21) | as.numeric(demo_prnt_race_a_v2_l___22) | as.numeric(demo_prnt_race_a_v2_l___23) | as.numeric(demo_prnt_race_a_v2_l___24)),
    demog_parhisp = as.numeric(demo_prnt_ethn_v2_l),
    demog_parnateng = as.numeric(ifelse(demo_prnt_nat_lang_l == "58", 1, 0)),
    demog_parmar = as.numeric(ifelse(demo_prnt_marital_v2_l == "1", 1, 0)),
    demog_edu = as.numeric(demo_prnt_ed_v2_l), 
    demog_work = as.numeric(ifelse(demo_prnt_empl_v2_l == "1", 1, 0)),
    demog_income = as.numeric(demo_comb_income_v2_l),
    demog_food = as.numeric(demo_fam_exp1_v2_l), 
    demog_phone = as.numeric(demo_fam_exp2_v2_l), 
    demog_rent = as.numeric(demo_fam_exp3_v2_l), 
    demog_evict = as.numeric(demo_fam_exp4_v2_l), 
    demog_offserv = as.numeric(demo_fam_exp5_v2_l), 
    demog_doc = as.numeric(demo_fam_exp6_v2_l), 
    demog_dent = as.numeric(demo_fam_exp7_v2_l), 
    demog_reltime = as.numeric(demo_yrs_1_l), 
    demog_relimp = as.numeric(demo_yrs_2_l),
    demog_relalc = as.numeric(demo_yrs_2a_l),
    demog_reldrug = as.numeric(demo_yrs_2b_l),
    demog_insured = as.numeric(as.numeric(demo_med_insur_a_p) | as.numeric(demo_med_insur_b_p) | as.numeric(demo_med_insur_c_p) | as.numeric(demo_med_insur_d_p) | as.numeric(demo_med_insur_e_p)),
    .keep="unused") %>%
    select(c(subjectkey, starts_with("demog_"))) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("demog_")), mean, na.rm=T)

abcd_lssmh01 <- abcd_lssmh01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("medhx_"))) %>%
    mutate_at(., vars(contains("medhx_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("medhx_")), mean, na.rm=T)

abcd_lsssa01 <- abcd_lsssa01 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("sai_ss_"))) %>%
    mutate_at(., vars(contains("sai_ss_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("sai_ss_")), mean, na.rm=T)

QUE <- abcd_abcls01 %>%
    full_join(crpbi01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_bpmt01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_cb01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_eatqp01, by=c("subjectkey"="subjectkey")) %>%
    full_join(cct01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_gdss01, by=c("subjectkey"="subjectkey")) %>%
    full_join(lmtp201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_cna01, by=c("subjectkey"="subjectkey")) %>%
    full_join(dhx01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_hsss01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_lssmh01, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_lsssa01, by=c("subjectkey"="subjectkey")) %>%
    full_join(fhxp102, by=c("subjectkey"="subjectkey")) %>%
    full_join(fhxp201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_lpds01, by=c("subjectkey"="subjectkey")) 

saveRDS(QUE, file="./CSI/Pitstop_QUE.rds")

# dMRI
abcd_drsip101 <- abcd_drsip101 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    select(-dmri_rsi_visitid) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip201 <- abcd_drsip201 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip301 <- abcd_drsip301 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip401 <- abcd_drsip401 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip501 <- abcd_drsip501 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip601 <- abcd_drsip601 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_drsip701 <- abcd_drsip701 %>% 
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_ddtidp101 <- abcd_ddtidp101 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("ddtidp_"))) %>%
    mutate_at(., vars(contains("ddtidp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("ddtidp_")), mean, na.rm=T)

abcd_ddtidp201 <- abcd_ddtidp201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("ddtidp_"))) %>%
    mutate_at(., vars(contains("ddtidp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("ddtidp_")), mean, na.rm=T)

abcd_ddtifp101 <- abcd_ddtifp101 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("ddtifp_"))) %>%
    mutate_at(., vars(contains("ddtifp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("ddtifp_")), mean, na.rm=T)

abcd_ddtifp201 <- abcd_ddtifp201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("ddtifp_"))) %>%
    mutate_at(., vars(contains("ddtifp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("ddtifp_")), mean, na.rm=T)

abcd_dmdtifp101 <- abcd_dmdtifp101 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmdtifp1_"))) %>%
    mutate_at(., vars(contains("dmdtifp1_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmdtifp1_")), mean, na.rm=T)

abcd_dmdtifp201 <- abcd_dmdtifp201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmdtifp1_"))) %>%
    mutate_at(., vars(contains("dmdtifp1_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmdtifp1_")), mean, na.rm=T)

abcd_dti_p101 <- abcd_dti_p101 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    select(-dmri_dti_visitid) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

abcd_dti_p201 <- abcd_dti_p201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("dmri_"))) %>%
    mutate_at(., vars(contains("dmri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("dmri_")), mean, na.rm=T)

dMRI <- abcd_drsip101 %>%
    full_join(abcd_drsip201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_drsip301, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_drsip401, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_drsip501, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_drsip601, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_drsip701, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_ddtidp101, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_ddtidp201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_ddtifp101, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_ddtifp201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_dmdtifp101, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_dmdtifp201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_dti_p101, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_dti_p201, by=c("subjectkey"="subjectkey"))

saveRDS(dMRI, file="./CSI/Pitstop_dMRI.rds")

# sMRI 
abcd_mrisdp10201 <- abcd_mrisdp10201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("mrisdp_"))) %>%
    mutate_at(., vars(contains("mrisdp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("mrisdp_")), mean, na.rm=T)

abcd_mrisdp20201 <- abcd_mrisdp20201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("mrisdp_"))) %>%
    mutate_at(., vars(contains("mrisdp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("mrisdp_")), mean, na.rm=T)

abcd_mrisdp30201 <- abcd_mrisdp30201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("mrisdp_"))) %>%
    mutate_at(., vars(contains("mrisdp_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("mrisdp_")), mean, na.rm=T)

abcd_smrip10201 <- abcd_smrip10201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("smri_"))) %>%
    select(-smri_visitid) %>%
    mutate_at(., vars(contains("smri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("smri_")), mean, na.rm=T)

abcd_smrip20201 <- abcd_smrip20201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("smri_"))) %>%
    mutate_at(., vars(contains("smri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("smri_")), mean, na.rm=T)

abcd_smrip30201 <- abcd_smrip30201 %>%
    filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1",        "2_year_follow_up_y_arm_1")) %>%
    select(c(subjectkey, starts_with("smri_"))) %>%
    mutate_at(., vars(contains("smri_")), as.numeric) %>%
    group_by(subjectkey) %>%
    summarise_at(., vars(contains("smri_")), mean, na.rm=T)

sMRI <- abcd_mrisdp10201 %>%
    full_join(abcd_mrisdp20201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_mrisdp30201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_smrip10201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_smrip20201, by=c("subjectkey"="subjectkey")) %>%
    full_join(abcd_smrip30201, by=c("subjectkey"="subjectkey"))

saveRDS(sMRI, file="./CSI/Pitstop_sMRI.rds")

# Preprocess genetic dataset
GEN <- GEN %>%
    select(-V2)

colnames(GEN) <- c("SampleID", SAM$V1)

GENt <- as_tibble(cbind(SampleID = names(GEN), t(GEN)))
colnames(GENt) <- GENt[1,]
GENt <- GENt[-1, ] 

saveRDS(GENt, file="./CSI/Pitstop_GENt.rds")

# Keep only samples with information on both
KEEP <- intersect(intersect(intersect(intersect(
    GENt$SampleID, OCD$V1), 
    dMRI$subjectkey),
    sMRI$subjectkey),
    QUE$subjectkey)

OCD <- OCD %>%
    filter(V1 %in% KEEP)

GENt <- GENt %>%
    filter(SampleID %in% KEEP)

dMRI <- dMRI %>%
    filter(subjectkey %in% KEEP)

sMRI <- sMRI %>%
    filter(subjectkey %in% KEEP)

QUE <- QUE %>%
    filter(subjectkey %in% KEEP)

########### PAUSE LINE

# Change class to numeric where applicable and save
GENt <- GENt %>%
    mutate_at(vars(rs1806446:rs111926263), as.numeric)
GENt[,-1] <- apply(GENt[,-1], 2, function(x) {(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))})
saveRDS(GENt, file="./CSI/Genetic.rds")

# Update OCD colnames and save phenotype data
colnames(OCD) <- c("SampleID", "OCD")
OCD$OCD[OCD$OCD == 1] <- NA
OCD$OCD[OCD$OCD == 2] <- 1
saveRDS(OCD, file="./CSI/OCD.rds")

# Update QUE colname and save que data
QUE <- QUE %>%
    rename(SampleID = subjectkey) %>%
    rename_at(vars(-SampleID),function(x) paste0("VAR_", x))
QUE <- QUE[colSums(!is.na(QUE)) > 0.2]
QUE[,-1] <- apply(QUE[,-1], 2, function(x) {(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))})
for(i in 2:ncol(QUE)){
  QUE[is.na(QUE[,i]), i] <- colMeans(QUE[,i], na.rm = TRUE)
}
sum(is.na(QUE))
saveRDS(QUE, file="./CSI/QUE.rds")

# Update and save dMRI data
dMRI <- dMRI %>%
    rename(SampleID = subjectkey) %>%
    rename_at(vars(-SampleID),function(x) paste0("VAR_", x))
dMRI <- dMRI[colSums(!is.na(dMRI)) > 0]
dMRI[,-1] <- apply(dMRI[,-1], 2, function(x) {(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))})
for(i in 2:ncol(dMRI)){
  dMRI[is.na(dMRI[,i]), i] <- colMeans(dMRI[,i], na.rm = TRUE)
}
sum(is.na(dMRI))
saveRDS(dMRI, file="./CSI/dMRI.rds")

# Update and save sMRI data
sMRI <- sMRI %>%
    rename(SampleID = subjectkey) %>%
    rename_at(vars(-SampleID),function(x) paste0("VAR_", x))
sMRI <- sMRI[colSums(!is.na(sMRI)) > 0.2]
sMRI[,-1] <- apply(sMRI[,-1], 2, function(x) {(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))})
for(i in 2:ncol(sMRI)){
  sMRI[is.na(sMRI[,i]), i] <- colMeans(sMRI[,i], na.rm = TRUE)
}
sum(is.na(sMRI))
saveRDS(sMRI, file="./CSI/sMRI.rds")

# Create Batches
BAT <- OCD %>%
    mutate(Experimental = ifelse(is.na(OCD), 1, 0))

BAT_Case <- OCD %>%
    filter(OCD==1) %>%
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
    mutate_at(vars(Experimental:Train_76), as.numeric) %>%
    rowwise() %>% 
    mutate(ALL = sum(c_across(all_of(names(BAT)[3:81]))),
        TRAIN = sum(c_across(all_of(names(BAT)[6:81]))))
saveRDS(BAT, file="./CSI/Batch.rds")
```
