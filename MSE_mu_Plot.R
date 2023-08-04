

wrap_plot_mse <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = mse.adj, color = "Copas HR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.neg.exp, color = "Neg Exp HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.wrong, color = "Copas HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.unadj, color = "Unadjusted"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = mse.adj.sigmoid, color = "Sigmoid HR+LR"), linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = expression(mu),
      y = expression(E((hat(mu) - mu)^2)),
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared,1)))
    ) +
    #ylim(-0.5, 0.5) +
    ylim(0, 0.6) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Generate plots
plots_mse <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_mse(.)) %>%
  pull(plot)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_mse, nrow = 3, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6)
)
