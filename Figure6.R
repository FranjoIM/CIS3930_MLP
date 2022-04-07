library(tidyverse)

# Kernel-wise analysis
DFA <- data.frame(
  KERNEL = rep(c("Total", "Phenotype", "sMRI", "dMRI", "rsfMRI", "tsfMRI", "Genotype"),2),
  VALUE = c(0.89, 0.90, 0.58, 0.53, 0.56, 0.49, 0.50, 0.93, 0.93, 0.63, 0.59, 0.60, 0.48, 0.50),
  TYPE = c(rep("Accuracy", 7), rep("Precision", 7)))

DFA$KERNEL <- factor(DFA$KERNEL, levels=c("Total", "Phenotype", "sMRI", "dMRI", "rsfMRI", "tsfMRI", "Genotype"))

DFB <- data.frame(
  KERNEL = rep(c("Total", "Phenotype", "sMRI", "dMRI", "rsfMRI", "tsfMRI", "Genotype"),2),
  VALUE = c(0.94, 0.94, 0.78, 0.86, 0.76, 0.78, 0.28, 0.84, 0.86, 0.38, 0.20, 0.36, 0.20, 0.72),
  TYPE = c(rep("Sensitivity", 7), rep("Specificity", 7)))

DFB$KERNEL <- factor(DFB$KERNEL, levels=c("Total", "Phenotype", "sMRI", "dMRI", "rsfMRI", "tsfMRI", "Genotype"))

ggplot(DFA, aes(x=KERNEL, y=VALUE, color=TYPE, fill=TYPE)) +
  geom_bar(stat="identity", position=position_dodge2(width = 0.1, preserve = "single", padding = 0.2)) +
  ylim(0,1) +
  geom_hline(yintercept = 0.5, color = "darkred", size = 1, linetype = "dashed") +
  theme_bw() +
  labs(fill="Measure",
       color="Measure",
       x=NULL, y=NULL) +
  theme(legend.position = "top")

ggplot(DFB, aes(x=KERNEL, y=VALUE, color=TYPE, fill=TYPE)) +
  geom_bar(stat="identity", position=position_dodge2(width = 0.1, preserve = "single", padding = 0.2)) +
  ylim(0,1) +
  geom_hline(yintercept = 0.5, color = "darkred", size = 1, linetype = "dashed") +
  theme_bw() +
  labs(fill="Measure",
       color="Measure",
       x=NULL, y=NULL) +
  theme(legend.position = "top")
