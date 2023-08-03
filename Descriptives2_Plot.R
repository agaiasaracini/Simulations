

# Function to wrap ggplot
wrap_plot_descr_sign <- function(data) {
  ggplot(data) +
    geom_line(aes(x = mu_values, y = perc_REP_nonsign.av, color = "% Rep Not Sign"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_REP_nonsign.av, color = "% Rep Not Sign")) +

    geom_line(aes(x = mu_values, y = perc_HR_sign.av, color = "% HR Sign"), linewidth = 0.6)+
    geom_point(aes(x = mu_values, y = perc_HR_sign.av, color = "% HR Sign"))+

    geom_line(aes(x = mu_values, y = perc_HR_nonsign.av, color = "% HR Not Sign"), linewidth = 0.6)+
    geom_point(aes(x = mu_values, y = perc_HR_nonsign.av, color = "% HR Not Sign"))+

    geom_line(aes(x = mu_values, y = perc_LR_sign.av, color = "% LR Sign"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_LR_sign.av, color = "% LR Sign")) +

    geom_line(aes(x = mu_values, y = perc_LR_nonsign.av, color = "% LR Not Sign"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_LR_nonsign.av, color = "% LR Not Sign")) +

    geom_line(aes(x = mu_values, y = perc_REP_sign.av, color = "% Rep Sign"), linewidth = 0.6) +
    geom_point(aes(x = mu_values, y = perc_REP_sign.av, color = "% Rep Sign")) +
    labs(
      x = expression(mu),
      y = "Average % of One-sided Significant",
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


legend_plot_descr_sign <- ggplot() +
  geom_line(aes(x = 0, y = 0, color =  "% Rep Not Sign")) +
  geom_line(aes(x = 0, y = 0, color =  "% HR Sign")) +
  geom_line(aes(x = 0, y = 0, color =  "% HR Not Sign")) +
  geom_line(aes(x = 0, y = 0, color = "% LR Sign"))+
  geom_line(aes(x = 0, y = 0, color = "% LR Not Sign"))+
  geom_line(aes(x = 0, y = 0, color = "% Rep Sign"))+

  theme_void()+
  theme(
    legend.text = element_text(size = 4),
    legend.key.size = unit(2, "mm"), # Adjust the size of plot title
    axis.title.x = element_blank(),             # Remove x-axis label
    axis.title.y = element_blank()
  )+
  guides(color = guide_legend(title = NULL))


# Generate plots
plots_descr_sign <- parameters %>%
  group_by(tau_squared_values, k_values) %>%
  do(plot_descr = wrap_plot_descr_sign(.)) %>%
  pull(plot_descr)

# Arrange the plots and legend
grid.arrange(
  arrangeGrob(grobs = plots_descr_sign, nrow = 3, ncol = 3),
  legend_plot_descr_sign,
  nrow = 1,
  widths = c(5, 0.6)
)

