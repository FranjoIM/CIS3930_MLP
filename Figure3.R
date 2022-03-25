# Figure 3

library(tidyverse)

# False positive Dx
DFA <- data.frame(
  KERNEL = c("Total", "Genetic", "dMRI", "sMRI", "Questionnaire"),
  VALUE = c(0.54, 0.61, 0.54, 0.61, 0.45, 0.53, 0.58, 0.53, 0.59, 0.46),
  TYPE = c(rep("Accuracy", 5), rep("Precision", 5)))

DFA$KERNEL <- factor(DFA$KERNEL, levels=c("Total", "Genetic", "dMRI", "sMRI", "Questionnaire"))

DFB <- data.frame(
  KERNEL = c("Total", "Genetic", "dMRI", "sMRI", "Questionnaire"),
  VALUE = c(0.76, 0.84, 0.64, 0.74, 0.56, 0.32, 0.38, 0.44, 0.48, 0.34),
  TYPE = c(rep("Sensitivity", 5), rep("Specificity", 5)))

DFB$KERNEL <- factor(DFB$KERNEL, levels=c("Total", "Genetic", "dMRI", "sMRI", "Questionnaire"))

ggplot(DFA, aes(x=KERNEL, y=VALUE, color=TYPE, fill=TYPE)) +
  geom_bar(stat="identity", position=position_dodge2(width = 0.1, preserve = "single", padding = 0.2)) +
  ylim(0,1) +
  theme_bw() +
  labs(fill="Measure",
       color="Measure",
       x=NULL, y=NULL) +
  theme(legend.position = "bottom")

ggplot(DFB, aes(x=KERNEL, y=VALUE, color=TYPE, fill=TYPE)) +
  geom_bar(stat="identity", position=position_dodge2(width = 0.1, preserve = "single", padding = 0.2)) +
  ylim(0,1) +
  theme_bw() +
  labs(fill="Measure",
       color="Measure",
       x=NULL, y=NULL) +
  theme(legend.position = "bottom")
