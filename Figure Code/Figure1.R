# Figure 1

library(tidyverse)

# False positive Dx
DF <- data.frame(
  NCASES             = c( 280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280),
  NCONTROLS          = c(1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120),
  MIS_OCD            = c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.12, 0.14, 0.15))

DF6 <- DF %>%
  mutate(AFREQ_CASE = 0.3, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF5 <- DF %>%
  mutate(AFREQ_CASE = 0.275, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF4 <- DF %>%
  mutate(AFREQ_CASE = 0.25, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF3 <- DF %>%
  mutate(AFREQ_CASE = 0.225, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF2 <- DF %>%
  mutate(AFREQ_CASE = 0.2, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF1 <- DF %>%
  mutate(AFREQ_CASE = 0.175, AFREQ_CONT = 0.1) %>%
  mutate(CASE_As = (1 - AFREQ_CASE) * NCASES * (1-MIS_OCD) + (1 - AFREQ_CONT) * NCASES * MIS_OCD,
         CASE_Bs = AFREQ_CASE * NCASES * (1-MIS_OCD) + AFREQ_CONT * NCASES * MIS_OCD,
         CONT_As = (1 - AFREQ_CONT) * NCONTROLS,
         CONT_Bs = AFREQ_CONT * NCONTROLS) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CASE_As, CASE_Bs, CONT_As, CONT_Bs), nrow=2, byrow=T))$p.value))

DF_MISSOCD <- rbind(DF1, DF2, DF3, DF4, DF5, DF6)

ggplot(DF_MISSOCD, aes(x=MIS_OCD, y=P, color=as.factor(AFREQ_CASE))) + # 520 by 500
  geom_line() +
  theme_bw() +
  xlab("OCD Misdiagnosis Rate (false positive)") +
  ylab("-log10 (P)") +
  theme(legend.position = "none")

# False negatives Dx
CF <- data.frame(
  NCASES             = c( 280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280,  280),
  NCONTROLS          = c(1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120),
  MIS_OCD            = c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.12, 0.14, 0.15))

CF6 <- CF %>%
  mutate(AFREQ_CA = 0.3, AFREQ_CO = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF5 <- CF %>%
  mutate(AFREQ_CA = 0.275, AFREQ_CO = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF4 <- CF %>%
  mutate(AFREQ_CA = 0.25, AFREQ_CO = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF3 <- CF %>%
  mutate(AFREQ_CA = 0.225, AFREQ_CO = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF2 <- CF %>%
  mutate(AFREQ_CA = 0.2, AFREQ_CO = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF1 <- CF %>%
  mutate(AFREQ_CO = 0.175, AFREQ_CA = 0.1) %>%
  mutate(CONT_A = (1 - AFREQ_CO) * NCONTROLS * (1-MIS_OCD) + (1 - AFREQ_CA) * NCONTROLS * MIS_OCD,
         CONT_B = AFREQ_CO * NCONTROLS * (1-MIS_OCD) + AFREQ_CA * NCONTROLS * MIS_OCD,
         CASE_A = (1 - AFREQ_CA) * NCASES,
         CASE_B = AFREQ_CA * NCASES) %>%
  rowwise() %>%
  mutate(P = -log10(fisher.test(matrix(c(CONT_A, CONT_B, CASE_A, CASE_B), nrow=2, byrow=T))$p.value))

CF_MISSOCD <- rbind(CF1, CF2, CF3, CF4, CF5, CF6)

ggplot(CF_MISSOCD, aes(x=MIS_OCD, y=P, color=as.factor(AFREQ_CA))) + # 520 by 500
  geom_line() +
  theme_bw() +
  xlab("OCD Misdiagnosis Rate (false negative)") +
  ylab("-log10 (P)") +
  theme(legend.position = "none")
