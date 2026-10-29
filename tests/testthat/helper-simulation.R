# Simulators used only by the tests. Neither produces realistic data.

# Closure along the first (feature) dimension of a D x N matrix.
miniclo <- function(X) {
  sweep(X, 2, colSums(X), FUN = "/")
}

# Linear-model simulation with the truth in CLR coordinates. No regression
# noise is added. Returns Y (D x N counts), X (P x N), W, and Lambda (P x D).
aldex.lm.sim.clr <- function(D = 10, N = 11, P = 2, depth = 10000) {
  X <- matrix(runif(P * N, -2, 2), P, N)
  Lambda <- matrix(rnorm(P * D), P, D)
  Lambda <- sweep(Lambda, 1, rowMeans(Lambda), FUN = `-`)
  W <- t(Lambda) %*% X
  pi <- miniclo(2^W)
  Y <- matrix(NA, D, N)
  for (n in 1:N) {
    Y[, n] <- rmultinom(1, 1000, pi[, n])
  }
  list(Y = Y, X = X, W = W, Lambda = Lambda)
}

# Mixed-effects simulation: subject random intercepts, an optional correlated
# random slope (corr), an optional second random intercept (location > 1),
# and optional AR(1) residual correlation (rho_ar1 > 0).
aldex.mem.sim <- function(D, days, subjects, depth = 10000, location = 1,
                          random_slope = FALSE, corr = 0, rho_ar1 = 0,
                          sd_resid = 0.1) {
  N <- days * subjects
  if (subjects < 1) {
    stop("subjects must be int > 0")
  }
  subject_ids <- rep(1:subjects, each = days)
  time <- rep(1:days, times = subjects)
  location_ids <- sample(1:location, N, replace = TRUE)
  treatment <- sample(c(0, 1), N, replace = TRUE)
  X <- cbind(1, treatment)
  meta <- data.frame(cbind(treatment, subject_ids, time, location_ids))

  W <- c()
  fixed_effs <- c()
  for (d in 1:D) {
    fixed_intercept <- runif(1, 2, 8)
    fixed_treatment <- runif(1, -2, 4)
    fixed_effs <- rbind(fixed_effs, c(fixed_intercept, fixed_treatment))
    sd_i <- 1.5
    sd_s <- ifelse(random_slope, 0.5, 0)
    sd_l <- ifelse(location > 1, 0.75, 0)
    Sigma_b <- rbind(c(sd_i^2, corr * sd_i * sd_s),
                     c(corr * sd_i * sd_s, sd_s^2))
    b <- MASS::mvrnorm(subjects, mu = c(0, 0), Sigma = Sigma_b)
    r_i <- b[, 1]
    r_s <- b[, 2]
    r_l <- rnorm(location, 0, sd_l^2)
    fixed <- fixed_intercept + fixed_treatment * treatment
    subject_ind <- as.integer(subject_ids)
    location_ind <- as.integer(location_ids)
    rand_eff <- r_i[subject_ind] + r_s[subject_ind] * treatment +
      r_l[location_ind]
    eps <- numeric(length(time))
    if (rho_ar1 > 0) {
      for (i in 1:subjects) {
        idx <- which(subject_ind == i)
        eps[idx] <- arima.sim(n = days, list(ar = rho_ar1), sd = sd_resid)
      }
    } else {
      eps <- rnorm(N, 0, sd_resid)
    }
    y <- fixed + rand_eff + eps
    W <- rbind(W, 2^y)
  }

  Y <- matrix(NA, D, N)
  for (n in 1:N) {
    Y[, n] <- rmultinom(1, depth, W[, n] / sum(W[, n]))
  }
  row.names(Y) <- paste0("taxa_", 1:D)
  colnames(Y) <- paste0("sample_", 1:N)
  row.names(meta) <- paste0("sample_", 1:N)
  row.names(fixed_effs) <- paste0("taxa_", 1:D)

  list(W = W, X = X, meta = meta, fixed_effs = fixed_effs, Y = Y)
}
