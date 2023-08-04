

wrap_plot_mse <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = mse.adj_hetfixed, color = "Copas HR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.neg.exp_hetfixed, color = "Neg Exp HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.wrong_hetfixed, color = "Copas HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.unadj_hetfixed, color = "Unadjusted"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.sigmoid_hetfixed, color = "Sigmoid HR+LR"), linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = expression(mu),
      y = expression(E(hat(mu)) - mu),
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared, 2)))
    ) +
    ylim(0, 0.35) +
    #ylim(0, max(data$mse.unadj)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Generate plots
plots_mse_hetfixed <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot_2hetfixed = wrap_plot_mse(.)) %>%
  pull(plot_2hetfixed)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_mse_hetfixed, nrow = 3, ncol = 3),
  legend_plot,
  nrow = 1,
  widths = c(5, 0.5)
)
