library("data.table")

cor_df <- read.csv("/home/estagaman/benchmarking_project/test_simulated_data/correlation_results/cor_values_default.csv")

cor_dt <- as.data.table(cor_df)

cor_melt <- melt(cor_dt, id.vars = "X", variable.name = "method", value.name = "cor")

cor_melt_df <- as.data.frame(cor_melt)

cor_melt_df$X <- factor(cor_melt_df$X, levels = unique(cor_melt_df$X))
cor_melt_df$method <- factor(cor_melt_df$method, levels = c("Slingshot", "PAGA", "monocle"))

# Now plot
heatmap <- ggplot(cor_melt_df, aes(x = X, y = method, fill = cor)) +
  geom_tile() +
  labs(title = "Correlation with True Pseudotime", x = "Dropout Rate", y = "Method") +
  scale_fill_gradient(low = "white", high = "#07075c")

ggsave("/home/estagaman/benchmarking_project/test_simulated_data/correlation_results/heatmap_default.png", heatmap, width = 10, height = 5)
