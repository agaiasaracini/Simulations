

library(foreach)
library(doParallel)
library(dplyr)
library(doRNG)

source("reORBgen.function.R")
source("simulate_meta_data.function.R")
source("reORBgen_hetfixed.function.R")


# Set the number of cores to use
num_cores <- 2
registerDoParallel(cores = num_cores)

# Define number of data sets to generate
n_datasets <- 2

#Parameter values
parameters <- expand.grid(
  gamma = 1.5,
  mu_values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
  k_values = c(5, 15, 30),
  tau_squared_values = c(0.01, 0.04, 0.37),

  stringsAsFactors = FALSE
)


seed <- 123



#For each parameter combinatin
for (i in seq_len(nrow(parameters))) {

  # get parameters (i.e. see grid above)
  pars <- parameters[i, ]
  K <- pars$k_values
  gamma <- pars$gamma
  tau_squared <- pars$tau_squared_values
  mu <- pars$mu_values

  result_df <- foreach(
    j = seq_len(n_datasets),
    .options.RNG = seed  # this is the seed
  ) %dorng% {

    #Calculations should go here: for every dataset, what I calculate
    #Create a dataset
    meta_data_miss <- simulate_meta_data(n_studies = K,
                                         mu = mu,
                                         tau_squared = tau_squared,
                                         n_treatment = 50,
                                         n_control = 50,
                                         gamma = gamma)


    #Descriptives: in each dataset
    perc_HR <- sum(meta_data_miss$y == "high")/K
    perc_LR <- sum(meta_data_miss$y == "low")/K
    perc_REP <- sum(meta_data_miss$y != "low" & meta_data_miss$y != "high")/K

    perc_sign <- sum(meta_data_miss$p_value <= 0.05)/K

    perc_HR_sign <- sum(meta_data_miss$p_value <= 0.05 & meta_data_miss$y == "high")/K
    perc_HR_nonsign <- sum(meta_data_miss$p_value >= 0.05 & meta_data_miss$y == "high")/K
    perc_LR_sign <- sum(meta_data_miss$p_value <= 0.05 & meta_data_miss$y == "low")/K
    perc_LR_nonsign <- sum(meta_data_miss$p_value >= 0.05 & meta_data_miss$y == "low")/K
    perc_REP_sign <- sum(meta_data_miss$p_value <= 0.05 & meta_data_miss$y != "low" & meta_data_miss$y != "high")/K
    perc_REP_nonsign <- sum(meta_data_miss$p_value >= 0.05 & meta_data_miss$y != "low" & meta_data_miss$y != "high")/K


    # Apply functions to dataset
    result <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE, #threshold 1.96
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Copas.oneside",
      opt_method = "Nelder-Mead",
      low.risk = FALSE
    )

    result_hetfixed <- reORBgen_hetfixed(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = 0.1,
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE, #threshold 1.96
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Copas.oneside",
      opt_method = "Brent",
      low.risk = FALSE,
      tau_squared_fixed = tau_squared
    )

    result.wrong <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Copas.oneside",
      opt_method = "Nelder-Mead",
      low.risk = TRUE
    )

    result.wrong_hetfixed <- reORBgen_hetfixed(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = 0.1,
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Copas.oneside",
      opt_method = "Brent",
      low.risk = TRUE,
      tau_squared_fixed = tau_squared
    )

    result.neg.exp <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Neg.exp.piecewise",
      weight_param = 7,
      opt_method = "L-BFGS-B",
      low.risk = TRUE
    )

    result.neg.exp_hetfixed <- reORBgen_hetfixed(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = 0.1,
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Neg.exp.piecewise",
      weight_param = 7,
      opt_method = "Brent",
      low.risk = TRUE,
      tau_squared_fixed = tau_squared
    )

    result.neg.exp3 <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Neg.exp.piecewise",
      weight_param = 3,
      opt_method = "L-BFGS-B",
      low.risk = TRUE
    )

    result.neg.exp30 <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Neg.exp.piecewise",
      weight_param = 30,
      opt_method = "L-BFGS-B",
      low.risk = TRUE
    )

    result.sigmoid <- reORBgen(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = c(0.1, 0.01),
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Sigmoid.cont",
      weight_param = 7,
      opt_method = "L-BFGS-B",
      low.risk = TRUE
    )

    result.sigmoid_hetfixed <- reORBgen_hetfixed(
      y = meta_data_miss$y,
      s = meta_data_miss$s,
      n1 = meta_data_miss$n_treatment,
      n2 = meta_data_miss$n_control,
      outcome = "benefit",
      init_param = 0.1,
      alpha_ben = 0.05,
      alpha_ben_one.sided = FALSE,
      true.SE = NULL,
      LR.CI = TRUE,
      selection.benefit = "Sigmoid.cont",
      weight_param = 7,
      opt_method = "Brent",
      low.risk = TRUE,
      tau_squared_fixed = tau_squared
    )





    #Put these things in a dataframe
    data.frame(

      perc_HR = perc_HR,
      perc_LR = perc_LR,
      perc_REP = perc_REP,

      perc_sign = perc_sign,

      perc_HR_sign = perc_HR_sign,
      perc_HR_nonsign = perc_HR_nonsign,
      perc_LR_sign = perc_LR_sign,
      perc_LR_nonsign = perc_LR_nonsign,
      perc_REP_sign = perc_REP_sign,
      perc_REP_nonsign = perc_REP_nonsign,


      mu_unadj = result$mu_unadjusted,
      LR_low_unadj = result$LR_mu_unadjusted_low,
      LR_up_unadj = result$LR_mu_unadjusted_up,
      LR_unadj_yes = (mu >= result$LR_mu_unadjusted_low && mu <= result$LR_mu_unadjusted_up),
      width_unadj = abs(result$LR_mu_unadjusted_up - result$LR_mu_unadjusted_low),
      tau2_ML_unadj = result$tau_squared_unadjusted,
      tau2_REML_unadj = result$tau_squared_unadjusted_REML,

      mu_unadj_hetfixed = result_hetfixed$mu_unadjusted,
      LR_low_unadj_hetfixed = result_hetfixed$LR_mu_unadjusted_low,
      LR_up_unadj_hetfixed = result_hetfixed$LR_mu_unadjusted_up,
      LR_unadj_yes_hetfixed = (mu >= result_hetfixed$LR_mu_unadjusted_low && mu <= result_hetfixed$LR_mu_unadjusted_up),
      width_unadj_hetfixed = abs(result_hetfixed$LR_mu_unadjusted_up - result_hetfixed$LR_mu_unadjusted_low),

      mu_adj = result$mu_adjusted_benefit,
      LR_low_adj = result$LR_mu_adjusted_low,
      LR_up_adj = result$LR_mu_adjusted_up,
      LR_adj_yes = (mu >= result$LR_mu_adjusted_low && mu <= result$LR_mu_adjusted_up),
      width_adj = abs(result$LR_mu_adjusted_up - result$LR_mu_adjusted_low),
      tau2_ML_adj = result$tau_squared_adjusted,
      tau2_REML_adj = result$tau_squared_adjusted_REML,

      mu_adj_hetfixed = result_hetfixed$mu_adjusted_benefit,
      LR_low_adj_hetfixed = result_hetfixed$LR_mu_adjusted_low,
      LR_up_adj_hetfixed = result_hetfixed$LR_mu_adjusted_up,
      LR_adj_yes_hetfixed = (mu >= result_hetfixed$LR_mu_adjusted_low && mu <= result_hetfixed$LR_mu_adjusted_up),
      width_adj_hetfixed = abs(result_hetfixed$LR_mu_adjusted_up - result_hetfixed$LR_mu_adjusted_low),


      mu_adj_neg.exp = result.neg.exp$mu_adjusted_benefit,
      LR_low_adj_neg.exp = result.neg.exp$LR_mu_adjusted_low,
      LR_up_adj_neg.exp = result.neg.exp$LR_mu_adjusted_up,
      LR_adj_neg.exp_yes = (mu >= result.neg.exp$LR_mu_adjusted_low && mu <= result.neg.exp$LR_mu_adjusted_up),
      width_adj_neg.exp = abs(result.neg.exp$LR_mu_adjusted_up - result.neg.exp$LR_mu_adjusted_low),
      tau2_ML_adj_neg.exp = result.neg.exp$tau_squared_adjusted,
      tau2_REML_adj_neg.exp = result.neg.exp$tau_squared_adjusted_REML,


      mu_adj_neg.exp_hetfixed = result.neg.exp_hetfixed$mu_adjusted_benefit,
      LR_low_adj_neg.exp_hetfixed = result.neg.exp_hetfixed$LR_mu_adjusted_low,
      LR_up_adj_neg.exp_hetfixed = result.neg.exp_hetfixed$LR_mu_adjusted_up,
      LR_adj_neg.exp_yes_hetfixed = (mu >= result.neg.exp_hetfixed$LR_mu_adjusted_low && mu <= result.neg.exp_hetfixed$LR_mu_adjusted_up),
      width_adj_neg.exp_hetfixed = abs(result.neg.exp_hetfixed$LR_mu_adjusted_up - result.neg.exp_hetfixed$LR_mu_adjusted_low),

      mu_adj_neg.exp3 = result.neg.exp3$mu_adjusted_benefit,
      LR_low_adj_neg.exp3 = result.neg.exp3$LR_mu_adjusted_low,
      LR_up_adj_neg.exp3 = result.neg.exp3$LR_mu_adjusted_up,
      LR_adj_neg.exp_yes3 = (mu >= result.neg.exp3$LR_mu_adjusted_low && mu <= result.neg.exp3$LR_mu_adjusted_up),
      width_adj_neg.exp3 = abs(result.neg.exp3$LR_mu_adjusted_up - result.neg.exp3$LR_mu_adjusted_low),
      tau2_ML_adj_neg.exp3 = result.neg.exp3$tau_squared_adjusted,
      tau2_REML_adj_neg.exp3 = result.neg.exp3$tau_squared_adjusted_REML,

      mu_adj_neg.exp30 = result.neg.exp30$mu_adjusted_benefit,
      LR_low_adj_neg.exp30 = result.neg.exp30$LR_mu_adjusted_low,
      LR_up_adj_neg.exp30 = result.neg.exp30$LR_mu_adjusted_up,
      LR_adj_neg.exp_yes30 = (mu >= result.neg.exp30$LR_mu_adjusted_low && mu <= result.neg.exp30$LR_mu_adjusted_up),
      width_adj_neg.exp30 = abs(result.neg.exp30$LR_mu_adjusted_up - result.neg.exp30$LR_mu_adjusted_low),
      tau2_ML_adj_neg.exp30 = result.neg.exp30$tau_squared_adjusted,
      tau2_REML_adj_neg.exp30 = result.neg.exp30$tau_squared_adjusted_REML,

      mu_adj_sigmoid = result.sigmoid$mu_adjusted_benefit,
      LR_low_adj_sigmoid = result.sigmoid$LR_mu_adjusted_low,
      LR_up_adj_sigmoid = result.sigmoid$LR_mu_adjusted_up,
      LR_adj_sigmoid = (mu >= result.sigmoid$LR_mu_adjusted_low && mu <= result.sigmoid$LR_mu_adjusted_up),
      width_adj_sigmoid = abs(result.sigmoid$LR_mu_adjusted_up - result.sigmoid$LR_mu_adjusted_low),
      tau2_ML_adj_sigmoid = result.sigmoid$tau_squared_adjusted,
      tau2_REML_adj_sigmoid = result.sigmoid$tau_squared_adjusted_REML,

      mu_adj_sigmoid_hetfixed = result.sigmoid_hetfixed$mu_adjusted_benefit,
      LR_low_adj_sigmoid_hetfixed = result.sigmoid_hetfixed$LR_mu_adjusted_low,
      LR_up_adj_sigmoid_hetfixed = result.sigmoid_hetfixed$LR_mu_adjusted_up,
      LR_adj_sigmoid_hetfixed = (mu >= result.sigmoid_hetfixed$LR_mu_adjusted_low && mu <= result.sigmoid_hetfixed$LR_mu_adjusted_up),
      width_adj_sigmoid_hetfixed = abs(result.sigmoid_hetfixed$LR_mu_adjusted_up - result.sigmoid_hetfixed$LR_mu_adjusted_low),

      mu_adj_wrong = result.wrong$mu_adjusted_benefit,
      LR_low_adj_wrong = result.wrong$LR_mu_adjusted_low,
      LR_up_adj_wrong = result.wrong$LR_mu_adjusted_up,
      LR_adj_wrong_yes = (mu >= result.wrong$LR_mu_adjusted_low && mu <= result.wrong$LR_mu_adjusted_up),
      width_adj_wrong = abs(result.wrong$LR_mu_adjusted_up - result.wrong$LR_mu_adjusted_low),
      tau2_ML_adj_wrong = result.wrong$tau_squared_adjusted,
      tau2_REML_adj_wrong = result.wrong$tau_squared_adjusted_REML,

      mu_adj_wrong_hetfixed = result.wrong_hetfixed$mu_adjusted_benefit,
      LR_low_adj_wrong_hetfixed = result.wrong_hetfixed$LR_mu_adjusted_low,
      LR_up_adj_wrong_hetfixed = result.wrong_hetfixed$LR_mu_adjusted_up,
      LR_adj_wrong_yes_hetfixed = (mu >= result.wrong_hetfixed$LR_mu_adjusted_low && mu <= result.wrong_hetfixed$LR_mu_adjusted_up),
      width_adj_wrong_hetfixed = abs(result.wrong_hetfixed$LR_mu_adjusted_up - result.wrong_hetfixed$LR_mu_adjusted_low),


      av_sigma_sq_true = mean(meta_data_miss$s_true^2)


    )

  }

  # save the seed as an attribute to the resulting object
  attr(result_df, "seed") <- seed
  #result_df


  # Calculate true sigma squared across all simulations
  parameters[i, "av_sigma_squared_true"] <- mean(unlist(lapply(result_df, function(df) df$av_sigma_sq_true)))
  #True I squared which ofc depends on which tau squared we have
  parameters[i, "I_squared"] <- mean(tau_squared / (tau_squared + unlist(lapply(result_df, function(df) df$av_sigma_sq_true))))

  #Mean descriptives
  parameters[i, "perc_HR.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_HR)))
  parameters[i, "perc_LR.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_LR)))
  parameters[i, "perc_REP.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_REP)))
  parameters[i, "perc_sign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_sign)))
  parameters[i, "perc_HR_sign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_HR_sign)))
  parameters[i, "perc_HR_nonsign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_HR_nonsign)))
  parameters[i, "perc_LR_sign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_LR_sign)))
  parameters[i, "perc_LR_nonsign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_LR_nonsign)))
  parameters[i, "perc_REP_sign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_REP_sign)))
  parameters[i, "perc_REP_nonsign.av"] <- mean(unlist(lapply(result_df, function(df) df$perc_REP_nonsign)))

  # Calculate summary statistics unadjusted
  parameters[i, "bias.unadj"] <- mean(unlist(lapply(result_df, function(df) df$mu_unadj))) - mu
  parameters[i, "cov.unadj"] <- sum(unlist(lapply(result_df, function(df) df$LR_unadj_yes == 1))) / n_datasets
  parameters[i, "width.unadj"] <- mean(unlist(lapply(result_df, function(df) df$width_unadj)))
  parameters[i, "mse.unadj"] <- mean(unlist(lapply(result_df, function(df) (df$mu_unadj - mu) ^ 2)))
  parameters[i, "bias.tau2.unadj.ML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_unadj))) - tau_squared
  parameters[i, "bias.tau2.unadj.REML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_unadj))) - tau_squared

  #Calculate summary statistics unadjusted het fixed
  parameters[i, "bias.unadj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$mu_unadj_hetfixed))) - mu
  parameters[i, "cov.unadj_hetfixed"] <- sum(unlist(lapply(result_df, function(df) df$LR_unadj_yes_hetfixed == 1))) / n_datasets
  parameters[i, "width.unadj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$width_unadj_hetfixed)))
  parameters[i, "mse.unadj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) (df$mu_unadj_hetfixed - mu) ^ 2)))

  # Calculate summary statistics adjusted Copas HR
  parameters[i, "bias.adj"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj))) - mu
  parameters[i, "cov.adj"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_yes == 1))) / n_datasets
  parameters[i, "width.adj"] <- mean(unlist(lapply(result_df, function(df) df$width_adj)))
  parameters[i, "mse.adj"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.ML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj))) - tau_squared
  parameters[i, "bias.tau2.adj.REML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj))) - tau_squared

  # Calculate summary statistics adjusted Copas HR het fixed
  parameters[i, "bias.adj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_hetfixed))) - mu
  parameters[i, "cov.adj_hetfixed"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_yes_hetfixed == 1))) / n_datasets
  parameters[i, "width.adj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_hetfixed)))
  parameters[i, "mse.adj_hetfixed"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_hetfixed - mu) ^ 2)))


  # Calculate summary statistics adjusted Neg Exp rho=7
  parameters[i, "bias.adj.neg.exp"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_neg.exp))) - mu
  parameters[i, "cov.adj.neg.exp"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_neg.exp_yes == 1))) / n_datasets
  parameters[i, "width.adj.neg.exp"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_neg.exp)))
  parameters[i, "mse.adj.neg.exp"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_neg.exp - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.neg.exp.ML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj_neg.exp))) - tau_squared
  parameters[i, "bias.tau2.adj.neg.exp.REML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj_neg.exp))) - tau_squared


  # Calculate summary statistics adjusted Neg Exp rho=7, het fixed
  parameters[i, "bias.adj.neg.exp_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_neg.exp_hetfixed))) - mu
  parameters[i, "cov.adj.neg.exp_hetfixed"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_neg.exp_yes_hetfixed == 1))) / n_datasets
  parameters[i, "width.adj.neg.exp_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_neg.exp_hetfixed)))
  parameters[i, "mse.adj.neg.exp_hetfixed"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_neg.exp_hetfixed - mu) ^ 2)))

  # Calculate summary statistics adjusted Neg Exp rho=3
  parameters[i, "bias.adj.neg.exp3"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_neg.exp3))) - mu
  parameters[i, "cov.adj.neg.exp3"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_neg.exp_yes3 == 1))) / n_datasets
  parameters[i, "width.adj.neg.exp3"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_neg.exp3)))
  parameters[i, "mse.adj.neg.exp3"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_neg.exp3 - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.neg.exp.ML3"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj_neg.exp3))) - tau_squared
  parameters[i, "bias.tau2.adj.neg.exp.REML3"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj_neg.exp3))) - tau_squared



  # Calculate summary statistics adjusted Neg Exp rho=30
  parameters[i, "bias.adj.neg.exp30"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_neg.exp30))) - mu
  parameters[i, "cov.adj.neg.exp30"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_neg.exp_yes30 == 1))) / n_datasets
  parameters[i, "width.adj.neg.exp30"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_neg.exp30)))
  parameters[i, "mse.adj.neg.exp30"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_neg.exp30 - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.neg.exp.ML30"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj_neg.exp30))) - tau_squared
  parameters[i, "bias.tau2.adj.neg.exp.REML30"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj_neg.exp30))) - tau_squared


  # Calculate summary statistics adjusted Sigmoid
  parameters[i, "bias.adj.sigmoid"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_sigmoid))) - mu
  parameters[i, "cov.adj.sigmoid"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_sigmoid == 1))) / n_datasets
  parameters[i, "width.adj.sigmoid"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_sigmoid)))
  parameters[i, "mse.adj.sigmoid"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_sigmoid - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.sigmoid.ML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj_sigmoid))) - tau_squared
  parameters[i, "bias.tau2.adj.sigmoid.REML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj_sigmoid))) - tau_squared

  # Calculate summary statistics adjusted Sigmoid het fixed
  parameters[i, "bias.adj.sigmoid_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_sigmoid_hetfixed))) - mu
  parameters[i, "cov.adj.sigmoid_hetfixed"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_sigmoid_hetfixed == 1))) / n_datasets
  parameters[i, "width.adj.sigmoid_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_sigmoid_hetfixed)))
  parameters[i, "mse.adj.sigmoid_hetfixed"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_sigmoid_hetfixed - mu) ^ 2)))


  # Calculate summary statistics adjusted WRONG
  parameters[i, "bias.adj.wrong"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_wrong))) - mu
  parameters[i, "cov.adj.wrong"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_wrong == 1))) / n_datasets
  parameters[i, "width.adj.wrong"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_wrong)))
  parameters[i, "mse.adj.wrong"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_wrong - mu) ^ 2)))
  parameters[i, "bias.tau2.adj.wrong.ML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_ML_adj_wrong))) - tau_squared
  parameters[i, "bias.tau2.adj.wrong.REML"] <- mean(unlist(lapply(result_df, function(df) df$tau2_REML_adj_wrong))) - tau_squared

  # Calculate summary statistics adjusted WRONG het fixed
  parameters[i, "bias.adj.wrong_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$mu_adj_wrong_hetfixed))) - mu
  parameters[i, "cov.adj.wrong_hetfixed"] <- sum(unlist(lapply(result_df, function(df) df$LR_adj_wrong_hetfixed == 1))) / n_datasets
  parameters[i, "width.adj.wrong_hetfixed"] <- mean(unlist(lapply(result_df, function(df) df$width_adj_wrong_hetfixed)))
  parameters[i, "mse.adj.wrong_hetfixed"] <- mean(unlist(lapply(result_df, function(df) (df$mu_adj_wrong_hetfixed - mu) ^ 2)))



}

#relative file path
write.csv(parameters, "~/Simulations_new.csv")
