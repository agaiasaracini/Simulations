
#Function to wrap ggplot
wrap_plot_bias2 <- function(data) {
  ggplot(data, aes(x = as.factor(tau_squared_values), group=1)) +
    geom_line(aes(y = bias.adj, color = "Copas HR"), linewidth = 0.6) +
    geom_line(aes(y = bias.adj.neg.exp, color = "Neg Exp HR+LR"), linewidth = 0.6) +
    geom_line(aes(y = bias.adj.wrong, color = "Copas HR+LR"), linewidth = 0.6) +
    geom_line(aes(y = bias.unadj, color = "Unadjusted"), linewidth = 0.6) +
    geom_line(aes(y = bias.adj.sigmoid, color = "Sigmoid HR+LR"), linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = expression(tau^2),
      y = expression(E(hat(mu)) - mu),
      title = bquote("K" == .(data$k_values) ~ "," ~ mu == .(data$mu_values))
    ) +
    #ylim(-0.5, 0.5) +
    ylim(-0.8, 0.8) +
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
  filter(mu_values %in% c(0, 0.3, 0.5))

# Generate plots
plots_bias2 <- filtered_results %>%
  group_by(mu_values, k_values) %>%
  do(plot = wrap_plot_bias2(.)) %>%
  pull(plot)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_bias2, nrow = 3, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6)
)

