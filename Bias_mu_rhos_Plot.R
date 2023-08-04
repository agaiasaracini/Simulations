
# Function to wrap ggplot
wrap_plot_rhos <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = bias.adj.neg.exp, color = "Neg Exp Rho=7"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = bias.adj.neg.exp3, color = "Neg Exp Rho=3"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = bias.adj.neg.exp30, color = "Neg Exp Rho=30"), linewidth = 0.6) +
    geom_line(aes(x = mu_values, y = bias.adj, color = "Copas HR"), linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      x = expression(mu),
      y = expression(E(hat(mu)) - mu),
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared,1)))
    ) +
    ylim(-0.5, 0.5) +
    # ylim(min(data$bias.adj.neg.exp,
    #				 data$bias.adj.neg.exp3,
    #				 data$bias.adj.neg.exp.rho30) - 0.01,
    #		 max(data$bias.adj.neg.exp,
    #		 		data$bias.adj.neg.exp3,
    #		 		data$bias.adj.neg.exp.rho30) + 0.01) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}


# Create a separate plot for the legend
legend_plot_rhos <- ggplot() +
  geom_line(aes(x = 0, y = 0, color =  "Neg Exp Rho=3")) +
  geom_line(aes(x = 0, y = 0, color =  "Neg Exp Rho=7")) +
  geom_line(aes(x = 0, y = 0, color =  "Neg Exp Rho=30")) +
  geom_line(aes(x = 0, y = 0, color =  "Copas HR")) +
  theme_void()+
  theme(
    legend.text = element_text(size = 4),
    legend.key.size = unit(2, "mm"), # Adjust the size of plot title
    axis.title.x = element_blank(),             # Remove x-axis label
    axis.title.y = element_blank()
  )+
  guides(color = guide_legend(title = NULL))


# Generate plots
plots_rhos <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot_rho = wrap_plot_rhos(.)) %>%
  pull(plot_rho)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_rhos, nrow = 3, ncol = 3),
  legend_plot_rhos,
  nrow = 1,
  widths = c(5, 0.6)
)

