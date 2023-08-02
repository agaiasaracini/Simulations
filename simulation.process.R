

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
n_datasets <- 10

# Define parameters
gamma <- 1.5
tau_squared_values <- c(0.009, 0.039, 0.36)
k_values <- c(5, 15, 30)
mu_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)



#### Simulation ########
results_df <- foreach(tau_squared = tau_squared_values, .combine = rbind) %:%
  foreach(K = k_values, .combine = rbind) %:%
  foreach(mu = mu_values, .combine = rbind) %dopar% {

    # Create empty dataframe to store output
    random_mus_df <- foreach(i = 1:n_datasets, .combine = rbind) %dopar% {
      # Generate meta data with missing values

      set.seed(123+i)

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


      # Apply function to datasets
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

      # Apply function to datasets
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


    # Calculate true sigma squared across all simulations
    av_sigma_squared_true <- mean(random_mus_df$av_sigma_sq_true)
    #True I squared which ofc depends on which tau squared we have
    I_squared <- mean(tau_squared / (tau_squared + av_sigma_squared_true))

    #Mean descriptives
    perc_HR.av = mean(random_mus_df$perc_HR)
    perc_LR.av = mean(random_mus_df$perc_LR)
    perc_REP.av = mean(random_mus_df$perc_REP)
    perc_sign.av = mean(random_mus_df$perc_sign)
    perc_HR_sign.av = mean(random_mus_df$perc_HR_sign)
    perc_HR_nonsign.av = mean(random_mus_df$perc_HR_nonsign)
    perc_LR_sign.av = mean(random_mus_df$perc_LR_sign)
    perc_LR_nonsign.av = mean(random_mus_df$perc_LR_nonsign)
    perc_REP_sign.av = mean(random_mus_df$perc_REP_sign)
    perc_REP_nonsign.av = mean(random_mus_df$perc_REP_nonsign)


    # Calculate biases
    bias.unadj <- mean(random_mus_df$mu_unadj) - mu
    cov.unadj <- sum(random_mus_df$LR_unadj_yes == 1) / n_datasets
    width.unadj <- mean(random_mus_df$width_unadj)
    mse.unadj <- mean((random_mus_df$mu_unadj - mu) ^ 2)
    bias.tau2.unadj.ML <- mean(random_mus_df$tau2_ML_unadj) - tau_squared
    bias.tau2.unadj.REML <- mean(random_mus_df$tau2_REML_unadj) - tau_squared

    bias.unadj_hetfixed <- mean(random_mus_df$mu_unadj_hetfixed) - mu
    cov.unadj_hetfixed <- sum(random_mus_df$LR_unadj_yes_hetfixed == 1) / n_datasets
    width.unadj_hetfixed <- mean(random_mus_df$width_unadj_hetfixed)
    mse.unadj_hetfixed <- mean((random_mus_df$mu_unadj_hetfixed - mu) ^ 2)


    bias.adj <- mean(random_mus_df$mu_adj) - mu
    cov.adj <- sum(random_mus_df$LR_adj_yes == 1) / n_datasets
    width.adj <- mean(random_mus_df$width_adj)
    mse.adj <- mean((random_mus_df$mu_adj - mu) ^ 2)
    bias.tau2.adj.ML <- mean(random_mus_df$tau2_ML_adj) - tau_squared
    bias.tau2.adj.REML <- mean(random_mus_df$tau2_REML_adj) - tau_squared

    bias.adj_hetfixed <- mean(random_mus_df$mu_adj_hetfixed) - mu
    cov.adj_hetfixed <- sum(random_mus_df$LR_adj_yes_hetfixed == 1) / n_datasets
    width.adj_hetfixed <- mean(random_mus_df$width_adj_hetfixed)
    mse.adj_hetfixed <- mean((random_mus_df$mu_adj_hetfixed - mu) ^ 2)


    bias.adj.neg.exp <- mean(random_mus_df$mu_adj_neg.exp) - mu
    cov.adj.neg.exp <- sum(random_mus_df$LR_adj_neg.exp_yes == 1) / n_datasets
    width.adj.neg.exp <- mean(random_mus_df$width_adj_neg.exp)
    mse.adj.neg.exp <- mean((random_mus_df$mu_adj_neg.exp - mu) ^ 2)
    bias.tau2.adj.neg.exp.ML <- mean(random_mus_df$tau2_ML_adj_neg.exp) - tau_squared
    bias.tau2.adj.neg.exp.REML <- mean(random_mus_df$tau2_REML_adj_neg.exp) - tau_squared

    bias.adj.neg.exp_hetfixed <- mean(random_mus_df$mu_adj_neg.exp_hetfixed) - mu
    cov.adj.neg.exp_hetfixed <- sum(random_mus_df$LR_adj_neg.exp_yes_hetfixed == 1) / n_datasets
    width.adj.neg.exp_hetfixed <- mean(random_mus_df$width_adj_neg.exp_hetfixed)
    mse.adj.neg.exp_hetfixed <- mean((random_mus_df$mu_adj_neg.exp_hetfixed - mu) ^ 2)


    bias.adj.neg.exp3 <- mean(random_mus_df$mu_adj_neg.exp3) - mu
    cov.adj.neg.exp3 <- sum(random_mus_df$LR_adj_neg.exp_yes3 == 1) / n_datasets
    width.adj.neg.exp3 <- mean(random_mus_df$width_adj_neg.exp3)
    mse.adj.neg.exp3 <- mean((random_mus_df$mu_adj_neg.exp3 - mu) ^ 2)
    bias.tau2.adj.neg.exp.ML3 <- mean(random_mus_df$tau2_ML_adj_neg.exp3) - tau_squared
    bias.tau2.adj.neg.exp.REML3 <- mean(random_mus_df$tau2_REML_adj_neg.exp3) - tau_squared


    bias.adj.neg.exp30 <- mean(random_mus_df$mu_adj_neg.exp30) - mu
    cov.adj.neg.exp30 <- sum(random_mus_df$LR_adj_neg.exp_yes30 == 1) / n_datasets
    width.adj.neg.exp30 <- mean(random_mus_df$width_adj_neg.exp30)
    mse.adj.neg.exp30 <- mean((random_mus_df$mu_adj_neg.exp30 - mu) ^ 2)
    bias.tau2.adj.neg.exp.ML30 <- mean(random_mus_df$tau2_ML_adj_neg.exp30) - tau_squared
    bias.tau2.adj.neg.exp.REML30 <- mean(random_mus_df$tau2_REML_adj_neg.exp30) - tau_squared


    bias.adj.sigmoid <- mean(random_mus_df$mu_adj_sigmoid) - mu
    cov.adj.sigmoid <- sum(random_mus_df$LR_adj_sigmoid == 1) / n_datasets
    width.adj.sigmoid <- mean(random_mus_df$width_adj_sigmoid)
    mse.adj.sigmoid <- mean((random_mus_df$mu_adj_sigmoid - mu) ^ 2)
    bias.tau2.adj.sigmoid.ML <- mean(random_mus_df$tau2_ML_adj_sigmoid) - tau_squared
    bias.tau2.adj.sigmoid.REML <- mean(random_mus_df$tau2_REML_adj_sigmoid) - tau_squared

    bias.adj.sigmoid_hetfixed <- mean(random_mus_df$mu_adj_sigmoid_hetfixed) - mu
    cov.adj.sigmoid_hetfixed <- sum(random_mus_df$LR_adj_sigmoid_hetfixed == 1) / n_datasets
    width.adj.sigmoid_hetfixed <- mean(random_mus_df$width_adj_sigmoid_hetfixed)
    mse.adj.sigmoid_hetfixed <- mean((random_mus_df$mu_adj_sigmoid_hetfixed - mu) ^ 2)


    bias.adj.wrong <- mean(random_mus_df$mu_adj_wrong) - mu
    cov.adj.wrong <- sum(random_mus_df$LR_adj_wrong_yes == 1) / n_datasets
    width.adj.wrong <- mean(random_mus_df$width_adj_wrong)
    mse.adj.wrong <- mean((random_mus_df$mu_adj_wrong - mu) ^ 2)
    bias.tau2.adj.wrong.ML <- mean(random_mus_df$tau2_ML_adj_wrong) - tau_squared
    bias.tau2.adj.wrong.REML <- mean(random_mus_df$tau2_REML_adj_wrong) - tau_squared

    bias.adj.wrong_hetfixed <- mean(random_mus_df$mu_adj_wrong_hetfixed) - mu
    cov.adj.wrong_hetfixed <- sum(random_mus_df$LR_adj_wrong_yes_hetfixed == 1) / n_datasets
    width.adj.wrong_hetfixed <- mean(random_mus_df$width_adj_wrong_hetfixed)
    mse.adj.wrong_hetfixed <- mean((random_mus_df$mu_adj_wrong_hetfixed - mu) ^ 2)



    # Return the results
    data.frame(
      tau_squared = tau_squared,
      I_squared = I_squared,
      K = K,
      mu = mu,

      perc_HR.av = perc_HR.av,
      perc_LR.av = perc_LR.av,
      perc_REP.av = perc_REP.av,
      perc_sign.av = perc_sign.av,
      perc_HR_sign.av = perc_HR_sign.av,
      perc_HR_nonsign.av = perc_HR_nonsign.av,
      perc_LR_sign.av = perc_LR_sign.av,
      perc_LR_nonsign.av = perc_LR_nonsign.av,
      perc_REP_sign.av = perc_REP_sign.av,
      perc_REP_nonsign.av = perc_REP_nonsign.av,

      bias.unadj = bias.unadj,
      cov.unadj = cov.unadj,
      width.unadj = width.unadj,
      mse.unadj = mse.unadj,
      bias.tau2.unadj.ML = bias.tau2.unadj.ML,
      bias.tau2.unadj.REML = bias.tau2.unadj.REML,

      bias.unadj_hetfixed = bias.unadj_hetfixed,
      cov.unadj_hetfixed = cov.unadj_hetfixed,
      width.unadj_hetfixed = width.unadj_hetfixed,
      mse.unadj_hetfixed = mse.unadj_hetfixed,


      bias.adj = bias.adj,
      cov.adj = cov.adj,
      width.adj = width.adj,
      mse.adj = mse.adj,
      bias.tau2.adj.ML = bias.tau2.adj.ML,
      bias.tau2.adj.REML = bias.tau2.adj.REML,

      bias.adj_hetfixed = bias.adj_hetfixed,
      cov.adj_hetfixed = cov.adj_hetfixed,
      width.adj_hetfixed = width.adj_hetfixed,
      mse.adj_hetfixed = mse.adj_hetfixed,

      bias.adj.neg.exp = bias.adj.neg.exp,
      cov.adj.neg.exp = cov.adj.neg.exp,
      width.adj.neg.exp = width.adj.neg.exp,
      mse.adj.neg.exp = mse.adj.neg.exp,
      bias.tau2.adj.neg.exp.ML = bias.tau2.adj.neg.exp.ML,
      bias.tau2.adj.neg.exp.REML = bias.tau2.adj.neg.exp.REML,

      bias.adj.neg.exp_hetfixed = bias.adj.neg.exp_hetfixed,
      cov.adj.neg.exp_hetfixed = cov.adj.neg.exp_hetfixed,
      width.adj.neg.exp_hetfixed = width.adj.neg.exp_hetfixed,
      mse.adj.neg.exp_hetfixed = mse.adj.neg.exp_hetfixed,


      bias.adj.neg.exp3 = bias.adj.neg.exp3,
      cov.adj.neg.exp3 = cov.adj.neg.exp3,
      width.adj.neg.exp3 = width.adj.neg.exp3,
      mse.adj.neg.exp3 = mse.adj.neg.exp3,
      bias.tau2.adj.neg.exp.ML3 = bias.tau2.adj.neg.exp.ML3,
      bias.tau2.adj.neg.exp.REML3 = bias.tau2.adj.neg.exp.REML3,

      bias.adj.neg.exp30 = bias.adj.neg.exp30,
      cov.adj.neg.exp30 = cov.adj.neg.exp30,
      width.adj.neg.exp30 = width.adj.neg.exp30,
      mse.adj.neg.exp30 = mse.adj.neg.exp30,
      bias.tau2.adj.neg.exp.ML30 = bias.tau2.adj.neg.exp.ML30,
      bias.tau2.adj.neg.exp.REML30 = bias.tau2.adj.neg.exp.REML30,

      bias.adj.sigmoid = bias.adj.sigmoid,
      cov.adj.sigmoid = cov.adj.sigmoid,
      width.adj.sigmoid = width.adj.sigmoid,
      mse.adj.sigmoid = mse.adj.sigmoid,
      bias.tau2.adj.sigmoid.ML = bias.tau2.adj.sigmoid.ML,
      bias.tau2.adj.sigmoid.REML = bias.tau2.adj.sigmoid.REML,

      bias.adj.sigmoid_hetfixed = bias.adj.sigmoid_hetfixed,
      cov.adj.sigmoid_hetfixed = cov.adj.sigmoid_hetfixed,
      width.adj.sigmoid_hetfixed = width.adj.sigmoid_hetfixed,
      mse.adj.sigmoid_hetfixed = mse.adj.sigmoid_hetfixed,


      bias.adj.wrong = bias.adj.wrong,
      cov.adj.wrong = cov.adj.wrong,
      width.adj.wrong = width.adj.wrong,
      mse.adj.wrong = mse.adj.wrong,
      bias.tau2.adj.wrong.ML = bias.tau2.adj.wrong.ML,
      bias.tau2.adj.wrong.REML = bias.tau2.adj.wrong.REML,

      bias.adj.wrong_hetfixed = bias.adj.wrong_hetfixed,
      cov.adj.wrong_hetfixed = cov.adj.wrong_hetfixed,
      width.adj.wrong_hetfixed = width.adj.wrong_hetfixed,
      mse.adj.wrong_hetfixed = mse.adj.wrong_hetfixed


    )


  }



# Stop parallel processing
stopImplicitCluster()

# Combine the results into a single dataframe
results_df <- bind_rows(results_df)

#write.csv(results_df, "/Users/alessandrasaracini/Documents/ETH/Msc THESIS/Simulations_all_parallel_10k.csv")
