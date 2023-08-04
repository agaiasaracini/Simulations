
# Function to wrap ggplot
wrap_plot_descr <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = perc_HR.av, color = "% HR"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_HR.av, color = "% HR")) +
    geom_line(aes(x = mu_values, y = perc_LR.av, color = "% LR"), linewidth = 0.6)+
    geom_point(aes(x = mu_values, y = perc_LR.av, color = "% LR"))+
    geom_line(aes(x = mu_values, y = perc_REP.av, color = "% Rep"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_REP.av, color = "% Rep")) +
    labs(
      x = expression(mu),
      y = "Average % of LR/HR/Rep",
      title = bquote("K" == .(data$k_values) ~ "," ~ I^2 == .(round(data$I_squared,1)))
    ) +
    #ylim(-0.5, 0.5) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7)
    )
}


legend_plot_descr <- ggplot() +
  geom_line(aes(x = 0, y = 0, color =  "% HR")) +
  geom_line(aes(x = 0, y = 0, color =  "% LR")) +
  geom_line(aes(x = 0, y = 0, color = "% Rep"))+

  theme_void()+
  theme(
    legend.text = element_text(size = 4),
    legend.key.size = unit(2, "mm"), # Adjust the size of plot title
    axis.title.x = element_blank(),             # Remove x-axis label
    axis.title.y = element_blank()
  )+
  guides(color = guide_legend(title = NULL))


# Generate plots
plots_descr <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot_descr = wrap_plot_descr(.)) %>%
  pull(plot_descr)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_descr, nrow = 3, ncol = 3),
  legend_plot_descr,
  nrow = 1,
  widths = c(5, 0.5)
)

