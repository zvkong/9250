X <- matrix(c(1,0,0,0,
              1,1,0,0,
              0,0,0,0,
              0,0,0,0), nrow = 4, ncol = 4)

## calculate N(X) and D(X)
calculate_ND <- function(X) {
  N <- sum(X)
  D <- sum(abs(diff(X, 1, 1))) + sum(abs(diff(t(X), 1, 1)))
  return(list(N = N, D = D))
}

# Metropolis-Hastings algorithm with acceptance rate, N(X), and D(X) Q1
metropolis_hastings_q1 <- function(p, lambda, iterations, initial_X) {
  # Initialize
  X <- initial_X
  samples <- list()
  samples[[1]] <- X
  accept_count <- 0
  
  # Initialize vectors to track N(X) and D(X)
  N_values <- numeric(iterations)
  D_values <- numeric(iterations)
  
  # Calculate initial N(X) and D(X)
  ND <- calculate_ND(X)
  N_values[1] <- ND$N
  D_values[1] <- ND$D
  
  for (i in 2:iterations) {
    # Propose a new candidate matrix X' from Q1 proposal distribution
    X_candidate <- matrix(rbinom(16, 1, p), nrow = 4)
    
    # Calculate acceptance probability alpha
    candidate_ND <- calculate_ND(X_candidate)
    current_pi <- p^ND$N * (1-p)^(16-ND$N) * exp(-lambda * ND$D)
    candidate_pi <- p^candidate_ND$N * (1-p)^(16-candidate_ND$N) * exp(-lambda * candidate_ND$D)
    alpha <- min(1, candidate_pi / current_pi)
    
    # Accept or reject the candidate matrix
    if (runif(1) < alpha) {
      X <- X_candidate
      accept_count <- accept_count + 1
      ND <- candidate_ND
    }
    
    # Store N(X) and D(X) for the current sample
    N_values[i] <- ND$N
    D_values[i] <- ND$D
    
    # Store the sample
    samples[[i]] <- X
  }
  
  acceptance_rate <- accept_count / (iterations - 1)
  return(list(samples = samples, acceptance_rate = acceptance_rate, N_values = N_values, D_values = D_values))
}

# Metropolis-Hastings algorithm with acceptance rate, N(X), and D(X) Q2
metropolis_hastings_q2 <- function(p, lambda, iterations, initial_X) {
  # Initialize
  X <- initial_X
  samples <- list()
  samples[[1]] <- X
  accept_count <- 0
  
  # Initialize vectors to track N(X) and D(X)
  N_values <- numeric(iterations)
  D_values <- numeric(iterations)
  
  # Calculate initial N(X) and D(X)
  ND <- calculate_ND(X)
  N_values[1] <- ND$N
  D_values[1] <- ND$D
  
  for (i in 2:iterations) {
    # Propose a new candidate matrix X' from Q2 proposal distribution
    X_candidate <- X
    row <- sample(1:4, 1)
    col <- sample(1:4, 1)
    X_candidate[row, col] <- rbinom(1, 1, p)
    
    # Calculate acceptance probability alpha
    candidate_ND <- calculate_ND(X_candidate)
    current_pi <- p^ND$N * (1-p)^(16-ND$N) * exp(-lambda * ND$D)
    candidate_pi <- p^candidate_ND$N * (1-p)^(16-candidate_ND$N) * exp(-lambda * candidate_ND$D)
    alpha <- min(1, candidate_pi / current_pi)
    
    # Accept or reject the candidate matrix
    if (runif(1) < alpha) {
      X <- X_candidate
      accept_count <- accept_count + 1
      ND <- candidate_ND
    }
    
    # Store N(X) and D(X) for the current sample
    N_values[i] <- ND$N
    D_values[i] <- ND$D
    
    # Store the sample
    samples[[i]] <- X
  }
  
  acceptance_rate <- accept_count / (iterations - 1)
  return(list(samples = samples, acceptance_rate = acceptance_rate, N_values = N_values, D_values = D_values))
}

# Simulation parameters
p <- 0.8
lambda <- 0.5 # Try different values of lambda, like 0.5 and 1
iterations <- 100000
initial_X <- matrix(rbinom(16, 1, p), nrow = 4, ncol = 4)

# Run the simulation
set.seed(123)
results_q1 <- metropolis_hastings_q1(p, lambda, iterations, initial_X)
results_q2 <- metropolis_hastings_q2(p, lambda, iterations, initial_X)

# Estimate the probability of having all ones on the diagonal
diagonal_ones <- sapply(results$samples, function(X) { all(diag(X) == 1) })
probability_diagonal_ones <- mean(diagonal_ones)

# Print results
cat("Probability of all ones on the diagonal:", probability_diagonal_ones, "\n")
cat("Acceptance rate:", results$acceptance_rate, "\n")

results$acceptance_rate
plot(results_q1$N_values, type = 'l')
plot(results_q1$D_values, type = 'l')

plot(results_q2$N_values, type = 'l')
plot(results_q2$D_values, type = 'l')
