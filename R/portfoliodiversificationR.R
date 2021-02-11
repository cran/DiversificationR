#################################################################################################################################
#                                                                                                                               #
#                                                         Functions                                                             #
#                                                                                                                               #
#################################################################################################################################

# Required packages
# library(stats)

#################################################################################################################################
#                                                                                                                               #
#                                                           f_VaR                                                              #
#                                                                                                                               #
#################################################################################################################################

# Function computing the Cornish-Fisher expansion for the Value-at-Risk
# Reference: Cornish, E. A., & Fisher, R. A. (1938). Moments and cumulants in the specification of distributions. Revue de l'Institut International de Statistique / Review of the International Statistical Institute, 5 (4), 307.
#f_VaR <- function(v_input_data = v_return, b_input_var_modified = b_var_modified, input_prob = input_prob)

f_VaR <- function(v_input_data, b_input_var_modified, input_prob)
{
  # Initialization
  v_input_data <- as.vector(v_input_data)
  number_observations <- length(v_input_data)
  mu <- mean(v_input_data)
  sigma <- stats::sd(v_input_data)

  # Compute Kurtosis - Step 1
  K_sum <- sum((v_input_data - mu) ^ 4 / (sigma ^ 2 * (number_observations - 1) / number_observations) ^ 2)
  # Compute Kurtosis - Step 2: moment case
  K <- (1 / number_observations) * K_sum

  # Compute Skewness - Step 1
  S_sum <- sum((v_input_data - mu) ^ 3)

  # Compute Skewness - Step 2: Moment case
  S <- (1 / number_observations) * S_sum / ((sigma ^ 2) ^ (3 / 2))

  # Compute the input_prob quantile of normal distribution
  Zn <- stats::qnorm(input_prob)

    # Compute classic or Cornish-Fisher expansion for the Value-at-Risk
  if(b_input_var_modified == TRUE)
  {
    Z <- Zn + (1 / 6) * (Zn ^ 2 - 1) * S + (1 / 24) * (Zn ^ 3 - 3 * Zn) * (K - 3) - (1 / 36) * (2 * Zn ^ 3 - 5 * Zn) * (S ^ 2)
  }else{
    Z <- Zn
  }

  # Compute final VaR
  result <- mu + Z * sigma
  return(result)
}

#################################################################################################################################
#                                                                                                                               #
#                                                            f_SR                                                               #
#                                                                                                                               #
#################################################################################################################################

# Function computing the Sharpe Ratio or modified Sharpe ratio
# Reference: Bali, Turan G., Stephen J. Brown, and K. Ozgur Demirtas. "Do hedge funds outperform stocks and bonds?." Management Science 59.8 (2013): 1887-1903.
# Reference: Favre, Laurent, and José-Antonio Galeano. "Mean-modified value-at-risk optimization with hedge funds." The journal of alternative investments 5.2 (2002): 21-25.
# Reference: Gregoriou, Greg N., and Jean-Pierre Gueyie. "Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio." The Journal of wealth management 6.3 (2003): 77-83.
# Reference: Sharpe, William F. "The sharpe ratio." Journal of Portfolio Management 21.1 (1994): 49-58.
# Reference: Sharpe, William F. "Mutual fund performance." The Journal of business 39.1 (1966): 119-138.
#f_SR <- function(v_input_data_portfolio = v_portfolio_return, v_input_data_benchmark = v_benchmark_return, c_input_method = c_method, input_prob = prob)

f_SR <- function(v_input_data_portfolio, v_input_data_benchmark, c_input_method, input_prob)
{
  number_observations <- length(v_input_data_portfolio)
  v_D <- numeric(number_observations)
  v_risk_free_rate <- numeric(number_observations)
  v_D <- (v_input_data_portfolio - v_input_data_benchmark)
  D <- mean(v_D)

  # Sharpe (1994) - the ex-post Sharpe ratio
  if (c_input_method == "S"){
    S <- D / stats::sd(v_D)
    result <- S
  # Bali et al. (2013) - the modified Sharpe ratio
  }else if (c_input_method == "BBD"){
    mS <- D / f_VaR(v_input_data = v_D, b_input_var_modified = FALSE, input_prob = input_prob)
    result <- mS
  # Favre and Galeano (2002) ; Gregoriou and Gueyie (2003) - the modified Sharpe ratio
  }else if (c_input_method == "FG-GG"){
    mS <- D / f_VaR(v_input_data = v_D, b_input_var_modified = TRUE, input_prob = input_prob)
    result <- mS
  # Sharpe (1966) - the Reward-to-Variablity Ratio
  }else{
    v_risk_free_rate <- v_input_data_benchmark
    RV <- (v_input_data_portfolio - v_risk_free_rate) / stats::sd(v_input_data_portfolio)
    result <- RV
  }
  return(result)
}

#################################################################################################################################
#                                                                                                                               #
#                                                           f_RSRL                                                              #
#                                                                                                                               #
#################################################################################################################################

# Function computing the RSRL or mRSRL. The purpose of this function is to measure the (under)diversification potential of a given portfolio with its benchmark
# Reference: Calvet, L. E., Campbell, J. Y., & Sodini, P. (2007). Down or out: Assessing the welfare costs of household investment mistakes. Journal of Political Economy, 115(5), 707-747.
# Reference: Candelon, B., Fuerst, F., & Hasse, J-B. (2020). Diversification Potential in Real Estate Portfolios. Cambridge Working Paper.
f_RSRL <- function(v_input_data_portfolio, v_input_data_benchmark, b_input_RSRL_modified, input_prob)
{
  v_data_portfolio <- as.vector(v_input_data_portfolio)
  v_data_benchmark <- as.vector(v_input_data_benchmark)

  # Choice RSRL or mRSRL
  if (b_input_RSRL_modified == TRUE)
  {
    # mRSRL as in Candelon, Fuerst and Hasse (2020)
    mSh <- f_SR(v_input_data_portfolio = v_data_portfolio, v_input_data_benchmark = 0, c_input_method = "FG-GG", input_prob = input_prob)
    mSb <- f_SR(v_input_data_portfolio = v_data_benchmark, v_input_data_benchmark = 0, c_input_method = "FG-GG", input_prob = input_prob)
    mRSRL <- 1 - mSh / mSb
    result <- mRSRL
  }else{
    # RSRL as in Calvet, Campbell and Sodini (2007)
    Sh <- f_SR(v_input_data_portfolio = v_data_portfolio, v_input_data_benchmark = 0, c_input_method = "S", input_prob = input_prob)
    Sb <- f_SR(v_input_data_portfolio = v_data_benchmark, v_input_data_benchmark = 0, c_input_method = "S", input_prob = input_prob)
    RSRL <- 1 - Sh / Sb
    result <- RSRL
  }
  return(result)
}

#################################################################################################################################
#                                                                                                                               #
#                                                   f_diversification_measurement                                               #
#                                                                                                                               #
#################################################################################################################################

# Function computing several diversification measures: Portfolio Diversification Index (PDI), Diversification Ratio (DR), Diversification Delta (DD) and Diversification Delta Star (DD*)
# Reference: Rudin, Alexander M. "A portfolio diversification index." The Journal of Portfolio Management 32.2 (2006): 81-89.
# Reference: Choueifaty, Yves, and Yves Coignard. "Toward maximum diversification." The Journal of Portfolio Management 35.1 (2008): 40-51.
# Reference: Vermorken, Maximilian A., Francesca R. Medda, and Thomas Schroder. "The diversification delta: A higher-moment measure for portfolio diversification." The Journal of Portfolio Management 39.1 (2012): 67-74.
# Reference: Flores, Yuri Salazar, et al. "The diversification delta: A different perspective." The Journal of Portfolio Management 43.4 (2017): 112-124.

# Function to measure portfolio diversification
f_diversification_measurement <- function(v_input_weights, m_input_returns, c_input_method)
{
  # Compute covariance and assets / portfolio variances
  m_returns <- m_input_returns / 100
  m_cov <- as.matrix(stats::cov(as.matrix(m_returns)))
  v_weights <- as.vector(v_input_weights)
  Nb_Assets <- length(v_input_weights)
  sigma_assets <- as.vector(apply(m_returns, 2, stats::sd))
  sigma_port <- sqrt((t(v_input_weights) %*% m_cov %*% v_input_weights))

  if (c_input_method == "Portfolio_Diversification_Index"){

    # Rudin, Alexander M. "A portfolio diversification index." The Journal of Portfolio Management 32.2 (2006): 81-89.
    # Interpretation: 1 < PDI < Nbre of Assets ; PDI = 1 -> maximum concentration ; PDI = Nbre of Assets -> maximum diversification
    pca <- stats::prcomp(m_cov, center = TRUE, scale = TRUE)
    pca_variance <- as.vector(apply(pca$x, 2, stats::var) )
    pca_proportion_of_variance <- as.vector(pca_variance / sum(pca_variance))
    k_Wk <- 0

    for (i in 1:Nb_Assets) {
      k_Wk <- k_Wk + sum(i * pca_proportion_of_variance[i])
    }
    PDI <- 2 * k_Wk - 1
    DM <- as.vector(PDI)


  }else if (c_input_method == "Diversification_Ratio"){

    # Choueifaty, Yves, and Yves Coignard. "Toward maximum diversification." The Journal of Portfolio Management 35.1 (2008): 40-51.
    # Interpretation: 1 < DR < sqrt(Nbre of Assets) ; DR = 1 -> maximum concentration ; DR = sqrt(Nbre of Assets) -> maximum diversification
    DM <- as.vector((t(v_input_weights) %*% (sigma_assets)) / sigma_port )

  }else if(c_input_method == "Diversification_Delta"){

    # Vermorken, Maximilian A., Francesca R. Medda, and Thomas Schroder. "The diversification delta: A higher-moment measure for portfolio diversification." The Journal of Portfolio Management 39.1 (2012): 67-74.
    # Interpretation: should be 0 < DD < 1 but Flores et al. (2017) show that DD can be negative ; should be DD = 0 -> maximum concentration ; DD = 1 -> maximum diversification
    H_X <- as.vector(0.5 * log(2 * pi * exp(1) * sigma_assets * sigma_assets))
    H_wX <- as.vector(0.5 * log(2 * pi * exp(1) * sigma_port * sigma_port))
    num <- exp(t(v_input_weights) %*% H_X) - exp(H_wX)
    den <- exp(v_input_weights %*% H_X)
    DM <- as.vector(num / den)

  }else if(c_input_method == "Diversification_Delta_Star"){

    # Flores, Yuri Salazar, et al. "The diversification delta: A different perspective." The Journal of Portfolio Management 43.4 (2017): 112-124.
    # Interpretation: 0 < DD* < 1 ; DD* = 0 -> maximum concentration ; DD* = 1 -> maximum diversification
    H_X <- as.vector(0.5 * log(2 * pi * exp(1) * sigma_assets * sigma_assets))
    H_wX <- as.vector(0.5 * log(2 * pi * exp(1) * sigma_port * sigma_port))
    num <- t(v_input_weights) %*% exp(H_X) - exp(H_wX)
    den <- t(v_input_weights) %*% exp(H_X)
    DM <- as.vector(num / den)

  }else{
    print("Method invalid")
    DM <- NULL
  }
  return(round(DM,3))
}

#################################################################################################################################
#                                                                                                                               #
#                                                     f_circular_bloc_bootstrap                                                 #
#                                                                                                                               #
#################################################################################################################################

# Function computing a circular block bootstrap
# Reference: Politis, D. N., & Romano, J. P. (1992). A circular block-resampling procedure for stationary data. Exploring the limits of bootstrap, 2635270.

f_circular_bloc_bootstrap <- function(m_input_data_series, input_c, input_b, input_prob)
{

  # Preliminary treatment
  number_series <- length(m_input_data_series[1,])
  length_series <- length(m_input_data_series[,1])
  m_output_series <- matrix(0, nrow = length_series, ncol = number_series)
  m_bootstrapped_series <- matrix(99.99, nrow = length_series, ncol = number_series)
  m_extended_input_series <- matrix(0, nrow = (length_series * input_c), ncol = number_series)

  # "circularization"
  m_extended_input_series <- do.call(rbind, replicate(input_c, m_input_data_series, simplify = FALSE))
  length_extended_series <- length(m_extended_input_series[,1])

  # Block length
  if (input_b == 0)
  {
    # Optimal block size from Hall, Horowitz and Jing (1995, Biometrika)
    optimal_size <- floor(length_series^(1/4))
  }else{
    optimal_size <- input_b
  }

  v_size_block <- c(optimal_size, number_series)
  number_block <- floor(length_series / v_size_block[1])

  # Bootstrap
  for(cpt in 1:number_block)
  {
    # sample
    target_left <- floor(stats::runif(1,1,length_extended_series - v_size_block[1] + 1))
    target_right <- v_size_block[1] + target_left - 1
    m_bootstrapped_series[(v_size_block[1] * cpt - v_size_block[1] + 1):(v_size_block[1] * cpt), ] <- m_extended_input_series[(target_left:target_right), ]
  }

  # Last rows assignations
  if ((length_series %% v_size_block[1]) != 0)
  {
    size_last_block <- length(m_bootstrapped_series[m_bootstrapped_series == 99.99]) / number_series
    last_target_left <- floor(stats::runif(1,1,length_extended_series - size_last_block + 1))
    last_target_right <- size_last_block + last_target_left - 1
    m_bootstrapped_series[((length_series - size_last_block + 1):length_series), ] <- m_extended_input_series[last_target_left:last_target_right,]
    m_output_series <- m_bootstrapped_series
  }else{
    m_output_series <- m_bootstrapped_series
  }

  # Compute bootstrapped RSRL and mRSRL
  bootstrapped_RSRL <- f_RSRL(v_input_data_portfolio = as.vector(m_output_series[,2]), v_input_data_benchmark = as.vector(m_output_series[,1]), b_input_RSRL_modified = FALSE, input_prob = input_prob)
  bootstrapped_mRSRL <- f_RSRL(v_input_data_portfolio = as.vector(m_output_series[,2]), v_input_data_benchmark = as.vector(m_output_series[,1]), b_input_RSRL_modified = TRUE, input_prob = input_prob)

  # Prepare results
  l_result <- list("RSRL" = bootstrapped_RSRL, "bootstapped_series" = m_output_series, "mRSRL" = bootstrapped_mRSRL)
  return(l_result)
}

#################################################################################################################################
#                                                                                                                               #
#                                                             f_test_RSRL                                                       #
#                                                                                                                               #
#################################################################################################################################

# Function estimating both coefficients and significance levels of RSRL and mRSRL. The purpose of this function is to test for diversification potential.
# Reference: Candelon, B., Fuerst, F., & Hasse, J-B. (2020). Diversification Potential in Real Estate Portfolios. Cambridge Working Paper.

f_test_RSRL <- function(v_input_p_r, v_input_b_r, input_c, input_b, input_sim, b_input_s, input_prob)
{

  # Compute bootstrapped RSRL
  v_bootstrapped_RSRL <- numeric(input_sim)
  v_bootstrapped_mRSRL <- numeric(input_sim)
  v_boot_RSRL <- numeric(input_sim)
  v_boot_mRSRL <- numeric(input_sim)

  # Number of simulations has to be >= 1000
  if (input_sim < 1000)
  {
    input_sim <- 1000
  }else{
    input_sim <- as.numeric(input_sim)
  }

  # Compute "real" RSRL
  RSRL <- f_RSRL(v_input_data_portfolio = v_input_p_r, v_input_data_benchmark = v_input_b_r, b_input_RSRL_modified = FALSE, input_prob = input_prob)

  # Compute "real" mRSRL
  mRSRL <- f_RSRL(v_input_data_portfolio = v_input_p_r, v_input_data_benchmark = v_input_b_r, b_input_RSRL_modified = TRUE, input_prob = input_prob)

  m_input_series <- cbind(v_input_b_r, v_input_p_r)

  # Call Circular Bloc Bootstrap in a loop
  for (cpt_simul in 1:input_sim)
  {
    l_res <- f_circular_bloc_bootstrap(m_input_data_series = m_input_series, input_c = input_c, input_b = input_b, input_prob = input_prob)
    v_bootstrapped_RSRL[cpt_simul] <- l_res$RSRL
    v_bootstrapped_mRSRL[cpt_simul] <- l_res$mRSRL
  }

  # Studentized or not ?
  if (b_input_s == TRUE)
  {

    # Studentized RSRL

    v_tau <- (v_bootstrapped_RSRL) / stats::sd(v_bootstrapped_RSRL)
    v_sorted_boot_RSRL <- sort(v_tau)

    # Compute 90, 95 and 99% confidence intervals

    right_side_3 <- v_sorted_boot_RSRL[(950/1000) * input_sim]
    left_side_3 <- v_sorted_boot_RSRL[(50/1000) * input_sim]

    right_side_2 <- v_sorted_boot_RSRL[(975/1000) * input_sim]
    left_side_2 <- v_sorted_boot_RSRL[(25/1000) * input_sim]

    right_side_1 <- v_sorted_boot_RSRL[(995/1000) * input_sim]
    left_side_1 <- v_sorted_boot_RSRL[(5/1000) * input_sim]

    # Compute significance levels (0 is or isn't in the confidence interval)

    if (left_side_1 > 0)
    {
      c_test <- "***"
    }else{
      if (left_side_2 > 0)
      {
        c_test <- "**"
      }else{
        if (left_side_3 > 0)
        {
          c_test <- "*"
        }else{
          c_test <- "."
        }
      }
    }

    # Studentized mRSRL

    v_tau <- (v_bootstrapped_mRSRL) / stats::sd(v_bootstrapped_mRSRL)
    v_sorted_boot_RSRL <- sort(v_tau)

    # Compute 90, 95 and 99% confidence intervals

    right_side_3 <- v_sorted_boot_RSRL[(950/1000) * input_sim]
    left_side_3 <- v_sorted_boot_RSRL[(50/1000) * input_sim]

    right_side_2 <- v_sorted_boot_RSRL[(975/1000) * input_sim]
    left_side_2 <- v_sorted_boot_RSRL[(25/1000) * input_sim]

    right_side_1 <- v_sorted_boot_RSRL[(995/1000) * input_sim]
    left_side_1 <- v_sorted_boot_RSRL[(5/1000) * input_sim]

    # Compute significance levels (0 is or isn't in the confidence interval)

    if (left_side_1 > 0)
    {
      c_test_2 <- "***"
    }else{
      if (left_side_2 > 0)
      {
        c_test_2 <- "**"
      }else{
        if (left_side_3 > 0)
        {
          c_test_2 <- "*"
        }else{
          c_test_2 <- "."
        }
      }
    }

  }else{

    # Non studentized RSRL

    v_boot_RSRL <- v_bootstrapped_RSRL
    v_sorted_boot_RSRL <- sort(v_boot_RSRL)

    # Compute 90, 95 and 99% confidence intervals

    right_side_3 <- v_sorted_boot_RSRL[(950/1000) * input_sim]
    left_side_3 <- v_sorted_boot_RSRL[(50/1000) * input_sim]

    right_side_2 <- v_sorted_boot_RSRL[(975/1000) * input_sim]
    left_side_2 <- v_sorted_boot_RSRL[(25/1000) * input_sim]

    right_side_1 <- v_sorted_boot_RSRL[(995/1000) * input_sim]
    left_side_1 <- v_sorted_boot_RSRL[(5/1000) * input_sim]

    # Compute significance levels (0 is or isn't in the confidence interval)

    if (left_side_1 > 0)
    {
      c_test <- "***"
    }else{
      if (left_side_2 > 0)
      {
        c_test <- "**"
      }else{
        if (left_side_3 > 0)
        {
          c_test <- "*"
        }else{
          c_test <- "."
        }
      }
    }

    # Non studentized mRSRL

    v_boot_mRSRL <- v_bootstrapped_mRSRL
    v_sorted_boot_mRSRL <- sort(v_boot_mRSRL)

    # Compute 90, 95 and 99% confidence intervals

    right_side_3 <- v_sorted_boot_mRSRL[(950/1000) * input_sim]
    left_side_3 <- v_sorted_boot_mRSRL[(50/1000) * input_sim]

    right_side_2 <- v_sorted_boot_mRSRL[(975/1000) * input_sim]
    left_side_2 <- v_sorted_boot_mRSRL[(25/1000) * input_sim]

    right_side_1 <- v_sorted_boot_mRSRL[(995/1000) * input_sim]
    left_side_1 <- v_sorted_boot_mRSRL[(5/1000) * input_sim]

    # Compute significance levels (0 is or isn't in the confidence interval)

    if (left_side_1 > 0)
    {
      c_test_2 <- "***"
    }else{
      if (left_side_2 > 0)
      {
        c_test_2 <- "**"
      }else{
        if (left_side_3 > 0)
        {
          c_test_2 <- "*"
        }else{
          c_test_2 <- "."
        }
      }
    }

  }
  l_result <- list("RSRL" = RSRL, "Signif_level_RSRL" = c_test, "mRSRL" = mRSRL, "Signif_level_mRSRL" = c_test_2)
  return(l_result)
}
