

wrap_plot_t2 <- function(data) {
  ggplot(data, aes(x = as.factor(tau_squared_values), group=1)) +
    geom_line(aes( y = bias.tau2.unadj.ML, color = "Unadjusted ML"), linewidth = 0.6) +
    geom_line(aes(y = bias.tau2.unadj.REML, color = "Unadjusted REML"), linewidth = 0.6) +

    geom_line(aes( y = bias.tau2.adj.ML, color = "Copas ML"), linewidth = 0.6) +
    geom_line(aes( y = bias.tau2.adj.REML, color = "Copas REML"), linewidth = 0.6) +

    geom_line(aes(y = bias.tau2.adj.neg.exp.ML, color = "Neg Exp ML"), linewidth = 0.6) +
    geom_line(aes( y = bias.tau2.adj.neg.exp.REML, color = "Neg Exp REML"), linewidth = 0.6) +

    geom_line(aes( y = bias.tau2.adj.sigmoid.ML, color = "Sigmoid ML"), linewidth = 0.6) +
    geom_line(aes(y = bias.tau2.adj.sigmoid.REML, color = "Sigmoid REML"), linewidth = 0.6) +

    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.3) +
    labs(
      x = expression(tau^2),
      y = expression(E(hat(tau^2)) - tau^2),
      title = bquote("K" == .(data$k_values) ~ "," ~ mu == .(as.numeric(data$mu_values)))
    ) +
    scale_x_discrete(labels = c("0.01", "0.04", "0.37"))+

    ylim(-1, 1) +
    #ylim(min(data$bias.unadj_ML,
    #				 data$bias.unadj_REML,
    #				 data$bias.adj_ML,
    #				 data$bias.adj_REML,
    #				 data$bias.adj.neg.exp_ML,
    #				 data$bias.adj.neg.exp_REML,
    #				 data$bias.adj.sigmoid_ML,
    #				 data$bias.adj.sigmoid_REML),

    #		 max(data$bias.unadj_ML,
    #		 		data$bias.unadj_REML,
  #		 		data$bias.adj_ML,
  #		 		data$bias.adj_REML,
  #		 		data$bias.adj.neg.exp_ML,
  #		 		data$bias.adj.neg.exp_REML,
  #		 		data$bias.adj.sigmoid_ML,
  #		 		data$bias.adj.sigmoid_REML)) +
  theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Filter results_df based on specific values of mu
filtered_results <- parameters %>%
  filter(mu_values %in% c(0, 0.1, 0.3, 0.5))

# Generate plots
plots_bias_t2 <- filtered_results %>%
  group_by(mu_values, k_values) %>%
  do(plot = wrap_plot_t2(.)) %>%
  pull(plot)


# Create a separate plot for the legend
legend_plot_t2 <- ggplot() +
  geom_line(aes(x = 0, y = 0, color = "Unadjusted ML")) +
  geom_line(aes(x = 0, y = 0, color = "Unadjusted REML")) +

  geom_line(aes(x = 0, y = 0, color = "Copas ML")) +
  geom_line(aes(x = 0, y = 0, color = "Copas REML")) +

  geom_line(aes(x = 0, y = 0, color = "Neg Exp ML")) +
  geom_line(aes(x = 0, y = 0, color = "Neg Exp REML")) +

  geom_line(aes(x = 0, y = 0, color = "Sigmoid ML")) +
  geom_line(aes(x = 0, y = 0, color = "Sigmoid REML")) +

  labs(color = "Heterogeneity Estimate") +
  theme_void()+
  theme(
    legend.text = element_text(size = 5),
    legend.key.size = unit(3, "mm") # Adjust the size of plot title
  )+
  guides(color = guide_legend(title = NULL))



# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_bias_t2, nrow = 4, ncol = 3),
  legend_plot_t2,
  nrow = 1,
  widths = c(5, 0.6)
)


