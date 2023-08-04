

wrap_plot_width <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = width.adj, color = "Copas HR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = width.adj.neg.exp, color = "Neg Exp HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = width.adj.wrong, color = "Copas HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = width.unadj, color = "Unadjusted"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = width.adj.sigmoid, color = "Sigmoid HR+LR"), linewidth = 0.6) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.3) +
    labs(
      x = expression(mu),
      y = "Widths of CI",
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared,1)))
    ) +
    ylim(0.1, 1.6) +
    # ylim(min(data$width.unadj,
    #				 data$width.adj,
    #				 data$width.adj.neg.exp,
    #				 data$width.adj.wrong,
    #				 data$width.adj.sigmoid),
    #		 max(data$width.unadj,
    #				 data$width.adj,
    #				 data$width.adj.neg.exp,
    #				 data$width.adj.wrong,
    #				 data$width.adj.sigmoid)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Generate plots
plots_width <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_width(.)) %>%
  pull(plot)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_width, nrow = 3, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.65)
)
