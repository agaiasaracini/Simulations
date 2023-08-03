

# Function to wrap ggplot
wrap_plot_cov <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = cov.adj, color = "Copas HR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = cov.adj.neg.exp, color = "Neg Exp HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = cov.adj.wrong, color = "Copas HR+LR"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = cov.unadj, color = "Unadjusted"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = cov.adj.sigmoid, color = "Sigmoid HR+LR"), linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = expression(mu),
      y = expression(E(hat(mu)) - mu),
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared,1)))
    ) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}

# Create a separate plot for the legend
legend_plot_bias <- ggplot() +
  geom_line(aes(x = 0, y = 0, color =  "Unadjusted")) +
  geom_line(aes(x = 0, y = 0, color =  "Copas HR")) +
  geom_line(aes(x = 0, y = 0, color =  "Neg Exp HR+LR")) +
  geom_line(aes(x = 0, y = 0, color = "Copas HR+LR"))+
  geom_line(aes(x = 0, y = 0, color = "Sigmoid HR+LR"))+

  theme_void()+
  theme(
    legend.text = element_text(size = 4),
    legend.key.size = unit(2, "mm"), # Adjust the size of plot title
    axis.title.x = element_blank(),             # Remove x-axis label
    axis.title.y = element_blank()
  )+
  guides(color = guide_legend(title = NULL))


# Generate plots
plots_cov <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot = wrap_plot_cov(.)) %>%
  pull(plot)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_cov, nrow = 3, ncol = 3),
  legend_plot_bias,
  nrow = 1,
  widths = c(5, 0.6)
)
