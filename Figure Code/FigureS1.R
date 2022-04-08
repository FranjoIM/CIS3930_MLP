# Load Libraries
library(tidyverse)

# Import data and change to appropriate shape
FMM <- read_csv("~/FinalModelMetrics.csv") %>%
  pivot_longer(cols=AUC:Precision,
               names_to = "WeighingMetric",
               values_to = "Value")

# Define colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7", "#0072B2")

# Plot
ggplot(FMM, aes(y=Value, fill=Model, x=WeighingMetric)) +
  geom_col(position="dodge") +
  coord_cartesian(ylim=c(0.50,0.95)) +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(. ~ Metric, scales = "free")
