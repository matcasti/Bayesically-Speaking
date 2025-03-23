# Full NUTS sampler
nuts <- function(log_posterior, grad_log_posterior, 
                 q0, Sigma, mu, 
                 epsilon, lambda_ta = 0.8, gamma = 0.5,
                 n_samples = 1000, n_warmup = 500) {
  # Initialization
  nuts_state <- nuts_init(q0, Sigma, mu, epsilon)
  samples <- matrix(0, nrow = n_samples + n_warmup, ncol = length(q0))
  
  # Warmup phase for step size adaptation
  for (i in 1:n_warmup) {
    nuts_state$p <- rnorm(length(q0), 0, 1) # Resample momentum
    trajectory <- build_trajectory(nuts_state$q, nuts_state$p, grad_log_posterior, nuts_state$epsilon, Sigma, mu)
    
    #Metropolis Hastings
    current_logp <- log_posterior(nuts_state$q, Sigma, mu)
    proposed_logp <- log_posterior(trajectory$q, Sigma, mu)
    
    acc_rate <- acceptance_prob(current_logp, proposed_logp)
    
    if(runif(1) < acc_rate){
      nuts_state$q <- trajectory$q
    }
    
    ## Adaptative epsilon
    nuts_state$epsilon <- adapt_epsilon(epsilon, acc_rate, lambda_ta, gamma)
    
    samples[i, ] <- nuts_state$q
  }
  
  # Sampling phase
  for (i in 1:n_samples) {
    nuts_state$p <- rnorm(length(q0), 0, 1) # Resample momentum
    trajectory <- build_trajectory(nuts_state$q, nuts_state$p, grad_log_posterior, nuts_state$epsilon, Sigma, mu)
    
    #Metropolis Hastings
    current_logp <- log_posterior(nuts_state$q, Sigma, mu)
    proposed_logp <- log_posterior(trajectory$q, Sigma, mu)
    
    if(runif(1) < acceptance_prob(current_logp, proposed_logp)){
      nuts_state$q <- trajectory$q
    }
    
    samples[i + n_warmup, ] <- nuts_state$q
  }
  
  samples
}

adapt_epsilon <- function(epsilon, lambda_ca, lambda_ta = 0.8, gamma = 0.5) {
  delta <- lambda_ca - lambda_ta
  epsilon <- epsilon * exp(gamma * delta)
  epsilon
}

# Leapfrog integration
leapfrog <- function(q, p, grad_func, epsilon, Sigma, mu) {
  grad_q <- grad_func(q, Sigma, mu) # Gradient calculation
  p <- p + (epsilon / 2) * grad_q # Half-step for momentum
  q <- q + epsilon * p # Full-step for position
  grad_q <- grad_func(q, Sigma, mu)
  p <- p + (epsilon / 2) * grad_q # Half-step for momentum
  list(q = q, p = p)
}

# U-turn check (dot product of momentum vectors)
is_uturn <- function(q_left, q_right, p_left, p_right) {
  (sum((q_right - q_left) * p_left) < 0) || 
    (sum((q_right - q_left) * p_right) < 0)
}

# Build a trajectory (doubling until U-turn)
build_trajectory <- function(q, p, grad_func, epsilon, Sigma, mu) {
  q_left <- q
  p_left <- p
  q_right <- q
  p_right <- p
  trajectory <- list(q = q, p = p)
  
  j <- 0
  while (!is_uturn(q_left, q_right, p_left, p_right) && j < 10) {
    if (runif(1) < 0.5) {
      # Expand left
      leapfrog_result <- leapfrog(q_left, p_left, grad_func, -epsilon, Sigma, mu)
      q_left <- leapfrog_result$q
      p_left <- leapfrog_result$p
    } else {
      # Expand right
      leapfrog_result <- leapfrog(q_right, p_right, grad_func, epsilon, Sigma, mu)
      q_right <- leapfrog_result$q
      p_right <- leapfrog_result$p
    }
    trajectory$q <- cbind(trajectory$q, leapfrog_result$q)
    trajectory$p <- cbind(trajectory$p, leapfrog_result$p)
    j <- j + 1
  }
  index <- sample(1:ncol(trajectory$q), 1)
  list(q = trajectory$q[, index], p = trajectory$p[, index])
}

# Initialization function of NUTS state
nuts_init <- function(q0, Sigma, mu, epsilon) {
  list(
    q = q0, # Initial position
    mu = mu,
    p = rnorm(length(q0), 0, 1), # Initial random momentum
    epsilon = epsilon # Initial step size (this will be adapted later)
  )
}

# Metropolis-Hastings acceptance probability
acceptance_prob <- function(current_logp, proposed_logp) {
  min(1, exp(proposed_logp - current_logp))
}

# Log-posterior function
log_posterior <- function(q, Sigma, mu) {
  diff <- q - mu
  -0.5 * t(diff) %*% solve(Sigma) %*% diff - 0.5 * log(det(Sigma))
}

## Gradient of Log-posterior
grad_log_posterior <- function(q, Sigma, mu) {
  -solve(Sigma, q - mu)
}

## Testing ====

q0 <- c(0, 0) ## Initial position
mu <- c(-1, 1) ## Target mu1 and mu2

Sigma <- matrix(c(1, 0, 0, 1), nrow = 2) ## Var-Cov Matrix

set.seed(12345) ## Seed for reproducibility
samples <- nuts(
  log_posterior = log_posterior, 
  grad_log_posterior = grad_log_posterior, 
  q0 = q0, Sigma = Sigma, mu = mu, epsilon = 0.1,
  n_samples = 1000, n_warmup = 1000
)

# Basic plot ----

## Traceplots
plot(samples[, 1], type = "l", col = rgb(0,0.5,0.5,1), ylim = c(-3,3), 
     axes = FALSE, xlab = "Iterations", ylab = expression(theta))
lines(samples[, 2], type = "l", col = rgb(0.5,0,0.5,1))
axis(1); axis(2); abline(h = c(-1,1), v = 1000, lty = 3)

## 2D posterior of mu1 and mu2
plot(samples, pch = 16, col = rgb(0,0.5,0.5,1/3), axes = FALSE, 
     xlab = expression(mu[1]), ylab = expression(mu[2]), 
     main = expression(P(theta~"|"~mu[1]*","*mu[2])))
axis(1); axis(2); abline(h = 1, v = -1, lty = 3)

