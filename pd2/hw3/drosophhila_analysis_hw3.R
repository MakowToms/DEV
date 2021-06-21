library(tidyverse)

df <- read_csv("drosophhila_net_node_table.csv")
ggplot(df, aes(x = EdgeCount)) +
  geom_histogram() +
  scale_x_log10() +
  scale_y_log10() +
  ylab("node count") +
  xlab("node degree") +
  ggtitle("Log node degree histogram")

df %>%
  mutate(BetweennessCentrality = BetweennessCentrality / max(BetweennessCentrality),
         ClosenessCentrality = ClosenessCentrality / max(ClosenessCentrality),
         MeanCentrality = (BetweennessCentrality + ClosenessCentrality) / 2) %>%
  arrange(-MeanCentrality) %>%
  select(name, SYMBOL_BD, MeanCentrality)
