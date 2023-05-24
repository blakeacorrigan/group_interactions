#####################################################################################################
#####################################################################################################
#####################################################################################################
# R Code for running Stan model, plotting and running simulation study for paper 
# "Inferring Stochastic Group Interactions within Structured Populations via Coupled Autoregression"
#####################################################################################################
#####################################################################################################
#####################################################################################################

# Clear the workspace
rm(list=ls())

###########
###########
# Packages
###########
###########

library(rstan); library(tidyverse); library(extraDistr)
library(bayesplot); library(MCMCvis); library(patchwork);
library(cowplot); library(RColorBrewer)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################
####################
# wolf-elk data-set
####################
####################

wolves_final <- read.csv(file = "Yellowstone.csv")

############################
############################
# Wolf-elk time series plot
############################
############################

wolves <- wolves_final %>%
  mutate(pack = as.factor(Pack)) %>%
  select(Year, Pack, Pack.Size, Number.Packs, Population.Size, Elk) %>%
  arrange(Pack)
wolves <- as.data.frame(wolves_final)
co <- count(wolves_final, Pack) 
cou <- co[co$n > 3,]

filter_wolves <- wolves_final %>%
  filter(Pack %in% cou$vars)

filter_wolf <- wolves_final %>%
  ggplot(aes(x = Year, y = Pack.Size))+
  theme_bw()+
  geom_line(aes(color = Pack), show.legend = F, size = 1)

real <- wolves_final[1:20,]

# wolf plot
main.plot <- wolves_final %>%
  ggplot(aes(x = Year, y = Pack.Size))+
  geom_line(aes(color = Pack), show.legend = F, size = 1)+
  scale_color_Yiridis_d(option="cividis")+
  theme(strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("Year")+
  ylab("Wolf Density")+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14))

# elk plot
inset.plot <- real %>%
  ggplot(aes(x = Year, y = Elk)) +
  theme_bw()+
  geom_line(aes(col = "#BDBDBD"), show.legend = F, size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Elk Density")+
  scale_colour_grey()+
  theme(strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14))

inset_main <- ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.65, y = 0.6, width = .33, height = .33)

##############################
##############################
# Code for running Stan model
##############################
##############################

# // is the prefix needed in Stan code for commenting

general = '
 data {
 // set up data inputs
   int<lower=1> time_year;
   int<lower=1> n_groups;
   int<lower=0> X[time_year, n_groups];
   int<lower=0> Y[time_year];
   int<lower=0> w_index[time_year, n_groups]; 
 }
 parameters {
    // elk parameters
   real omega_Y;
   real lambda_Y;
   real gamma;
   // monitor random effects
   real omega_X[n_groups];
   real lambda_X[n_groups];
   real psi[n_groups];
   real delta[n_groups];
   // mean hyperparameters
   real mu_omega_X;
   real mu_lambda_X;
   real mu_psi;
   real mu_delta;
    // sd hyperparameters
   real<lower = 0> sigma_omega_X;
   real<lower = 0> sigma_lambda_X;
   real<lower = 0> sigma_psi;
   real<lower = 0> sigma_delta;
 }
 model {
  for (i in 1:n_groups){
  // random effects
    omega_X[i] ~ normal(mu_omega_X, sigma_omega_X);
    lambda_X[i] ~ normal(mu_lambda_X, sigma_lambda_X);
    psi[i] ~ normal(mu_psi, sigma_psi);
    delta[i] ~ normal(mu_delta, sigma_delta);
  }
  for(t in 2:time_year){
  // elk mean process
    Y[t] ~  poisson_log(omega_Y + lambda_Y*log(Y[t-1] + 1) + gamma*log(sum(X[t-1,]) + 1));
    for(i in 1:n_groups){
    //(x = a ? b : c returns b if a is non-zero and c otherwise)
      // wolf pack mean processes
      X[t,i] ~ poisson_log(w_index[t,i] ? omega_X[i] + lambda_X[i]*log(X[t-1,i] + 1) + psi[i]*log(sum(X[t-1,]) - X[t-1,i] + 1) + delta[i]*log(Y[t-1] + 1) : 0.000000000001);
    // Stan does not allow degenerate distributions - cannot take 0 as Poisson mean 
    // (so use sufficiently small number,  0.000000000001, in replace of 0)
    
    }
  }
  
  // priors
  omega_Y ~ normal(0, 100);
  lambda_Y ~ normal(0, 100);
  gamma ~ normal(0, 100);
  mu_omega_X ~ normal(0, 100);
  mu_lambda_X ~ normal(0, 100);
  mu_psi ~ normal(0, 100);
  mu_delta ~ normal(0, 100);
  sigma_omega_X ~ normal(0, 100);
  sigma_lambda_X ~ normal(0, 100);
  sigma_psi ~ normal(0, 100);
  sigma_delta ~ normal(0, 100);
 }
'

# convert data-set to time series matrix
# w is the indicator C in paper
Y <- wolves_final$Elk[1:20]
X <- w <- matrix(0, ncol = 42, nrow = 20)
j <- 1
for (i in 1:42){
  X[, i] <- wolves_final$Pack.Size[j:(j+19)]
  w[,i] <- wolves_final$index[j:(j+19)] 
  j <- (j+20)
}

# run stan model (default is NUTS MCMC)
# can change chains etc as needed
stan_run = stan(model_code = general,
                data = list(
                  time_year = length(Y),
                  n_groups = ncol(X),
                  X = X, 
                  Y = Y,
                  w_index = w),
                chains = 4,
                iter = 20000,
                thin = 20,
                warmup = 2000,
                # algorithm = "HMC",
                pars = c("omega_Y",
                         "lambda_Y",
                         "gamma",
                         "mu_omega_X",
                         "mu_lambda_X",
                         "mu_psi",
                         "mu_delta",
                         "sigma_omega_X",
                         "sigma_lambda_X",
                         "sigma_psi",
                         "sigma_delta",
                         "omega_X",
                         "lambda_X",
                         "psi",
                         "delta")
)

# save data
# save(stan_run, file = "stan_wolf_elk_feb_2023.RData")

##########################
##########################
# Convergence diagnostics
##########################
##########################

# Gelman-Rubin statistic and effective sample size given in Stan output

# Trace plots
traceplot(stan_run,c("omega_Y",
                   "lambda_Y",
                   "gamma",
                   "mu_omega_X",
                   "mu_lambda_X",
                   "mu_psi",
                   "mu_delta",
                   "sigma_omega_X",
                   "sigma_lambda_X",
                   "sigma_psi",
                   "sigma_delta"), include = TRUE)

# Autocorrelation plots
stan_ac(stan_run, c("omega_Y",
                  "lambda_Y",
                  "gamma",
                  "mu_omega_X",
                  "mu_lambda_X",
                  "mu_psi",
                  "mu_delta",
                  "sigma_omega_X",
                  "sigma_lambda_X",
                  "sigma_psi",
                  "sigma_delta"), lags = 10)

####################
####################
# Empirical results
####################
####################

load("~/stan_output_wolf_elk.RData")

# retrieve mean estimates
output <- as.data.frame(get_posterior_mean(stan_run))
names_row <- rownames(output)
output <- output %>%
  rename(mean = "mean-all chains")

# retrieve estimated parameter values
index_omega_Y <- tail(which(grepl("omega_Y", names_row) == TRUE), n = 1)
omega_Y <- output$mean[index_omega_Y]
index_lambda_Y <- tail(which(grepl("lambda_Y", names_row) == TRUE), n = 1)
lambda_Y <- output$mean[index_lambda_Y]
index_gamma <- tail(which(grepl("gamma", names_row) == TRUE), n = 1)
gamma <- output$mean[index_gamma]
index_mu_psi <- tail(which(grepl("mu_psi", names_row) == TRUE), n = 1)
mu_psi <- output$mean[index_mu_psi]
index_mu_delta <- tail(which(grepl("mu_delta", names_row) == TRUE), n = 1)
mu_delta <- output$mean[index_mu_delta]
index_mu_omega_X <- tail(which(grepl("mu_omega_X", names_row) == TRUE), n = 1)
mu_omega_X <- output$mean[index_mu_omega_X]
index_sigma_omega_X <- tail(which(grepl("sigma_omega_X", names_row) == TRUE), n = 1)
sigma_omega_X <- output$mean[index_sigma_omega_X]
index_mu_lambda_X <- tail(which(grepl("mu_lambda_X", names_row) == TRUE), n = 1)
mu_lambda_X <- output$mean[index_mu_lambda_X]
index_sigma_lambda_X <- tail(which(grepl("sigma_lambda_X", names_row) == TRUE), n = 1)
sigma_lambda_X <- output$mean[index_sigma_lambda_X]
index_sigma_psi <- tail(which(grepl("sigma_psi", names_row) == TRUE), n = 1)
sigma_psi <- output$mean[index_sigma_psi]
index_sigma_delta <- tail(which(grepl("sigma_delta", names_row) == TRUE), n = 1)
sigma_delta <- output$mean[index_sigma_delta]

# retrieve random effect realisations
omega_X_true <- output$mean[4:(42+3)]
lambda_X_true <- output$mean[44:(44+43)]
psi_true <- output$mean[((44+43)+1):(((44+43)+1)+41)]
delta_true <- output$mean[(((44+43)+1)+42):((((44+43)+1)+41)+42)]

# construct parameter table with true values and estimates
Parameters <- data.frame(Parameter = c("omega_Y", "lambda_Y", "gamma", 
                                       "mu_omega_X", "sigma_omega_X",
                                       "mu_lambda_X", "sigma_lambda_X",
                                       "mu_psi", "sigma_psi",
                                       "mu_delta", "sigma_delta"),
                         Estimate = c(omega_Y, lambda_Y, gamma, 
                                      mu_omega_X,sigma_omega_X,
                                      mu_lambda_X, sigma_lambda_X,
                                      mu_psi, sigma_psi,
                                      mu_delta, sigma_delta))

# write.table(x = Parameters, file = paste("Parameters_Wolf_Elk", ".txt", sep = ""), sep = "\t")

################################
################################
# Prior-posterior overlap plots
################################
################################

# generate sufficiently large draws from normal distribution
# as set in paper as N(0,100^2)
PRn <- rnorm(30000, 0, 100)

# generate PPO plot with PRn and
# each parameter posterior densities
# to check for weak identifiability

MCMCtrace(stan_run,
          params = c("omega_Y", "lambda_Y", "gamma", 
                     "mu_omega_X", "sigma_omega_X",
                     "mu_lambda_X", "sigma_lambda_X",
                     "mu_psi", "sigma_psi",
                     "mu_delta", "sigma_delta"),
          type = 'density',
          ISB = FALSE,
          exact = TRUE,
          priors = PRn,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = FALSE,
          lwd_den = 2,
          lty_Xr = 1,
          lwd_Xr = 2,
          col_Xr = 'black',
          col_txt = 'black',
          col_den = 'grey',
          sz_txt = 1.2,
          sz_ax = 2,
          sz_ax_txt = 1.2,
          sz_tick_txt = 1.3,
          sz_main_txt = 1.4,
          xlab_den = '')

###############
###############
# Violin plots
###############
###############

#code below is for the random effect \psi
# code is the same for other three random effects

psi_draws <- tidybayes::tidy_draws(stan_run) %>%
  select(.chain, .iteration, .draw, starts_with("psi"))

both_draws_psi <- bind_rows(psi_draws[, #index_groups
] %>%
  gather(key = key, value = value, starts_with("psi")) %>%
  mutate(model = ""))

both_draws_psi <- both_draws_psi %>%
  rename("Estimate" = value,
         "Group" = key)

both_draws_psi %>%
  ggplot(mapping = aes(x = Group, y = Estimate)) +
  geom_Yiolin(draw_quantiles = c(0.025, 0.5, 0.975), col = "black", fill = "grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  #stat_summary(fun.y=mean, geom="point", shape=16, size=2, col = "#C55A11") +
  scale_x_discrete(labels = 1:42) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 9)) +
  ylab(expression(psi))

##########################
##########################
# Simulate from the model
##########################
##########################

# the following was used for each scenario S1 to S5
# 100 simulations for each scenario
# followed by model fitting using stan code above

# example (from S2):
n_groups <- 42
time.year <- 20

# fixed effects
omega_Y <- 0.9
lambda_Y <- 0.9
gamma <- -0.05
mu_psi <- -0.5
mu_delta <- 0.2
mu_omega_X <- 0.9
mu_lambda_X <- 0.06
sigma_lambda_X <- 0.2
sigma_omega_X <- 0.2
sigma_delta <- 0.04
sigma_psi <- 0.04

# initial population values
V_initial <- 10000
P_initial <- sample(1:20, n_groups, replace = T)
P <- mu_p <- a <- b <- c <- matrix(NA, ncol = n_groups, nrow = time.year)
V <- mu_Y <- psi <- delta <- omega_X <- lambda_X <- vector()
X[1,] <- P_initial
mu_p[1,] <- 0
Y[1] <- V_initial

# simulate random effects for each group
for (i in 1:n_groups){
  omega_X[i] <- rnorm(1,mu_omega_X, sigma_omega_X);
  lambda_X[i] <- rnorm(1,mu_lambda_X, sigma_lambda_X);
  psi[i] <- rnorm(1, mu_psi, sigma_psi);
  delta[i] <- rnorm(1,mu_delta, sigma_delta);
}

# loop over known observation period 
for(t in 2:time.year){
  # auxiliary abundance
  mu_Y[t] <- exp(omega_Y + lambda_Y*log(Y[t-1] + 1) + gamma*log(sum(X[t-1,]) + 1))
  Y[t] <- rpois(1, mu_Y[t])
  # group abundances
  for(i in 1:n_groups){
    a[t,i] <- lambda_X[i]*log(X[t-1,i] + 1)
    b[t,i] <- psi[i]*log(sum(X[t-1,]) - X[t-1,i] + 1)
    c[t,i] <- delta[i]*log(Y[t-1] + 1)
    mu_p[t,i] <- exp(omega_X[i] + a[t,i] + b[t,i] + c[t,i])
    X[t,i] <- rpois(1, mu_p[t,i])
  }
}

# SIMULATION DATA FROM PAPER CAN BE PROVIDED UPON REQUEST
# file sizes are large and so were not included in Data files

##############################
##############################
# Calculate RMSE, bias and RB
##############################
##############################

# parameter estimates
Est_matrix <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  Est_matrix[i,] <- param_table$Estimate
}

# true values
Par_matrix <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  Par_matrix[i,] <- param_table$Value
}

# rmse
RMSE_Yal <- vector()
for (i in 1:11){
  RMSE_Yal[i] <- sqrt(sum((Est_matrix[,i]-Par_matrix[,i])^2)/length(Par_matrix[,i]))
}

# median relative bias over 100 simulations
bias_Yal <- vector()
for (i in 1:11){
  bias_Yal[i] <- median((Par_matrix[,i]-Est_matrix[,i])/Par_matrix[,i])
}

# bias
raw_bias <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  raw_bias[i,] <- param_table$Est - param_table$Value
}

# relative bias
rel_bias <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  rel_bias[i,] <- (param_table$Est - param_table$Value)/param_table$Value
}

bias_matrix <- data.frame(Parameter = c("omega_Y", "lambda_Y", "gamma", 
                                        "mu_omega_X", "sigma_omega_X",
                                        "mu_lambda_X", "sigma_lambda_X",
                                        "mu_psi", "sigma_psi",
                                        "mu_delta", "sigma_delta"),
                          bias = bias_Yal,
                          RMSE = RMSE_Yal)
write.table(x = bias_matrix, file = "bias_large.txt", sep = "\t")

####################
####################
# RB and bias plots
####################
####################

# data matrix for bias plots
rawbox_data <- matrix(NA, nrow = 100, ncol = 11)

# the following loops over each simulation to generate plots needed
# to assess convergence and retrieve parameter values
for (i in 1:100){
  
  # load simulation data i from Stan model
  
  output_run <- as.data.frame(get_posterior_mean(stan_run))
  names_row_run <- rownames(output_run)
  output_run <- output_run %>%
    rename(mean = "mean-all chains")
  
  # retrieves relevant parameter from Stan output
  index_omega_Y_run <- tail(which(grepl("omega_Y", names_row_run) == TRUE), n = 1)
  omega_Y_run <- output_run$mean[index_omega_Y_run]
  index_lambda_Y_run <- tail(which(grepl("lambda_Y", names_row_run) == TRUE), n = 1)
  lambda_Y_run <- output_run$mean[index_lambda_Y_run]
  index_gamma_run <- tail(which(grepl("gamma", names_row_run) == TRUE), n = 1)
  gamma_run <- output_run$mean[index_gamma_run]
  index_mu_psi_run <- tail(which(grepl("mu_psi", names_row_run) == TRUE), n = 1)
  mu_psi_run <- output_run$mean[index_mu_psi_run]
  index_mu_delta_run <- tail(which(grepl("mu_delta", names_row_run) == TRUE), n = 1)
  mu_delta_run <- output_run$mean[index_mu_delta_run]
  index_mu_omega_X_run <- tail(which(grepl("mu_omega_X", names_row_run) == TRUE), n = 1)
  mu_omega_X_run <- output_run$mean[index_mu_omega_X_run]
  index_sigma_omega_X_run <- tail(which(grepl("sigma_omega_X", names_row_run) == TRUE), n = 1)
  sigma_omega_X_run <- output_run$mean[index_sigma_omega_X_run]
  index_mu_lambda_X_run <- tail(which(grepl("mu_lambda_X", names_row_run) == TRUE), n = 1)
  mu_lambda_X_run <- output_run$mean[index_mu_lambda_X_run]
  index_sigma_lambda_X_run <- tail(which(grepl("sigma_lambda_X", names_row_run) == TRUE), n = 1)
  sigma_lambda_X_run <- output_run$mean[index_sigma_lambda_X_run]
  index_sigma_psi_run <- tail(which(grepl("sigma_psi", names_row_run) == TRUE), n = 1)
  sigma_psi_run <- output_run$mean[index_sigma_psi_run]
  index_sigma_delta_run <- tail(which(grepl("sigma_delta", names_row_run) == TRUE), n = 1)
  sigma_delta_run <- output_run$mean[index_sigma_delta_run]
  
  boxplot_data[i,] <- c(omega_Y_run, lambda_Y_run, gamma_run, 
                        mu_omega_X_run,sigma_omega_X_run,
                        mu_lambda_X_run, sigma_lambda_X_run,
                        mu_psi_run, sigma_psi_run,
                        mu_delta_run, sigma_delta_run)
}

raw_bias <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  raw_bias[i,] <- param_table$Est - param_table$Value
}

rel_bias <- matrix(NA, ncol = 11, nrow = 100)
for (i in 1:100){
  param_table <- read.table(paste("Parameters_", i, ".txt", sep = ""))
  rel_bias[i,] <- (param_table$Est - param_table$Value)/param_table$Value
}

colnames(rel_bias) <- c("omega_Y", "lambda_Y", "gamma", 
                        "mu_omega", "sigma_omega_X",
                        "mu_lambda_X", "sigma_lambda_X",
                        "mu_psi", "sigma_psi",
                        "mu_delta", "sigma_delta")

rel_data_long <- as.data.frame(rel_bias) %>%
  gather(key = "parameter", value = "estimate", 1:11)

x_axis_labels <- expression(gamma, lambda[Y], mu[delta], mu[lambda], mu[omega], mu[psi],
                            omega[Y], sigma[delta], sigma[lambda], sigma[omega], sigma[psi])

boxplot_mod_rel <- ggplot(rel_data_long, aes(x = parameter, y = estimate)) +
  geom_boxplot(outlier.colour = "darkgreen",
               outlier.color = NULL,
               outlier.fill = NULL,
               outlier.shape = 16,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = NULL,
               col = "black",
               notch = F) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_y_continuous(limits = c(-1, 2.1), 
                     expand = c(0, 0.05)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.5, color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.spacing.y = unit(0.2, "lines"),
        legend.key = element_blank(),
        strip.background = element_blank()) +
  ggtitle("") +
  xlab("")+
  ylab(expression(frac(hat(theta) - theta,theta)))+
  scale_x_discrete(labels = x_axis_labels)+
  xlab("Parameter")

####################
####################
# Bivariate plots
####################
####################

# to check for confouning

# In paper example is given in Appendix for 
# simulation no. 20 for the parameters below

mcmc_pairs(as.array(stan_run), pars = c("omega_Y",
                                      "lambda_Y",
                                      "mu_psi",
                                      "mu_delta"),
           off_diag_args = list(size = 0.75))

