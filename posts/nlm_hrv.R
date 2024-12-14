library(data.table)
library(ggplot2)
library(brms)

# Simulate time points (e.g., every minute for 100 minutes)
time_points <- seq(0, 300, length.out = 500)

hrv_sim <- function(t, ## Vector of time measurements
                    baseline, ## Baseline HRV
                    drop, ## Initial drop after exercise
                    rate_drop, ## Rate of initial at exercise
                    recov_magnitude, ## Magnitude of recovery when exercise is stopped
                    recov_rate, ## Rate at which the HRV recovery occurred
                    t_0, ## Time at which the exercise beggins
                    t_1) { ## Time at which the exercise stops
  hrv <- 
    baseline + 
    (drop / (1 + exp(-rate_drop * (t - t_0)))) + 
    ((-drop*recov_magnitude) / (1 + exp(-recov_rate * (t - t_1))))
  
  out <- data.frame(hrv, time = t)
  return(out)
}

n_time <- 100
n_sample <- 10

set.seed(1234)
params <- list(
  baseline = rnorm(n_sample, 100, 5), 
  drop = rnorm(n_sample, -30, 5), 
  rate_drop = rnorm(n_sample, 0.25, 0.05),
  recov_magnitude = rnorm(n_sample, 1.5, 0.1), 
  recov_rate = rnorm(n_sample, 0.1, 0.025), 
  t_0 = rnorm(n_sample, 100, 5), 
  t_1 = rnorm(n_sample, 170, 5)
)

data <- vector("list", length = n_sample)
time_points <- seq(0, 300, length.out = n_time)

for (i in 1:n_sample) {
  data[[i]] <- hrv_sim(
    t = time_points,
    baseline = params$baseline[i],
    drop = params$drop[i],
    rate_drop = params$rate_drop[i],
    recov_magnitude = params$recov_magnitude[i],
    recov_rate = params$recov_rate[i],
    t_0 = params$t_0[i],
    t_1 = params$t_1[i]
  )
  
  data[[i]] <- data[[i]] + rnorm(n_time, 0, 2)
}

data <- rbindlist(data, idcol = "id")
data[, id := factor(id)]

## Plot data 
ggplot(data, aes(time, hrv, col = id)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "gam")

# Nonlinear formula: we want to estimate N (population size) and p (probability of sighting)
formula <- bf(hrv ~ baseline + 
                (drop / (1 + exp(-rateDrop * (time - t0)))) + 
                ((-drop*recov) / (1 + exp(-rateRecov * (time - t1)))),
              baseline + drop + rateDrop + recov + rateRecov + t0 + t1 ~ 1,
              nl = TRUE)  

priors <- c(
  prior(student_t(3, 100, 10), nlpar = "baseline", lb = 0),
  prior(student_t(3, -30, 5), nlpar = "drop", ub = 0),
  prior(normal(1, 0.5), nlpar = "recov", lb = 0),
  prior(normal(0.1, 0.1), nlpar = "rateDrop", lb = 0),
  prior(normal(0.1, 0.1), nlpar = "rateRecov", lb = 0),
  prior(student_t(3, 75, 30), nlpar = "t0", lb = 0, ub = 300),
  prior(student_t(3, 200, 30), nlpar = "t1", lb = 0, ub = 300)
)

# Fit the model
fit <- brm(
  formula = formula,
  data = data,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  seed = 1234,
  iter = 10000
)

# Check model summary
summary(fit)
conditional_effects(fit)
pp_check(fit, ndraws = 100)

# Plot posterior distributions
out <- plot(fit, nvariables = 4, plot = FALSE)

pairs(fit, diag_fun = "dens", off_diag_fun = "hex", regex_pars = "^b_")

