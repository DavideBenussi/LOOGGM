rm(list = ls())
gc()
set.seed(123)

library(MASS)
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(parallel)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
LRT <- function(data, grp) {
  X  <- as.matrix(data)
  n1 <- sum(grp == 1); n2 <- sum(grp == 2)
  p  <- ncol(X)
  
  if (p == 1) {
    s1 <- var(X[grp == 1]) * (n1 - 1) / n1
    s2 <- var(X[grp == 2]) * (n2 - 1) / n2
    sT <- var(X)           * (n1 + n2 - 1) / (n1 + n2)
    return(n1 * log(sT / s1) + n2 * log(sT / s2))
  } else {
    S1 <- var(X[grp == 1, , drop = FALSE]) * (n1 - 1) / n1
    S2 <- var(X[grp == 2, , drop = FALSE]) * (n2 - 1) / n2
    ST <- var(X) * (n1 + n2 - 1) / (n1 + n2)
    return(n1 * log(det(ST) / det(S1)) + n2 * log(det(ST) / det(S2)))
  }
}

mu_n <- function(p, n1, n2) {
  n  <- n1 + n2
  rx <- function(p, x) sqrt(-log(1 - p / x))
  (1/4) * (-(4*p + p/n1 + p/n2)
           + n  * rx(p, n)^2     * (2*p - 2*n  + 3)
           - (n1 * rx(p, n1-1)^2 * (2*p - 2*n1 + 3)
              + n2 * rx(p, n2-1)^2 * (2*p - 2*n2 + 3)))
}

bartlett_adj <- function(W, p, n1, n2) {
  df0  <- p * (p + 3) / 2
  corr <- df0 / (-2 * mu_n(p, n1, n2))
  W * corr
}

df_delta <- function(l, q) l * (2*q - l + 3) / 2

is_pd <- function(A) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) > 1e-8
}

# ------------------------------------------------------------
# Dense covariance generator
# Sigma = A A^T + rho * I
# ------------------------------------------------------------
make_A <- function(q = 8, A_seed = 123) {
  set.seed(A_seed)
  matrix(rnorm(q * q), nrow = q, ncol = q)
}

make_S1 <- function(q = 8,
                    rho = 0.4,
                    A = NULL,
                    standardize = TRUE) {
  
  if (is.null(A)) {
    A <- matrix(rnorm(q * q), nrow = q, ncol = q)
  }
  
  S <- A %*% t(A) + rho * diag(q)
  
  if (isTRUE(standardize)) {
    Dinv <- diag(1 / sqrt(diag(S)), q, q)
    S <- Dinv %*% S %*% Dinv
  }
  
  if (!is_pd(S)) stop("Generated S1 is not positive definite.")
  S
}

# ------------------------------------------------------------
# Build parameters
# mode = "node": mean/variance node perturbation
# mode = "edge": pure precision-edge change
# ------------------------------------------------------------
build_params <- function(q, k, dm, xi, S1, mu0 = 0,
                         perturbed_set = NULL,
                         mode = c("node", "edge"),
                         edge_change = NULL,
                         delta_edge = 0) {
  
  mode <- match.arg(mode)
  
  mu1 <- rep(mu0, q)
  mu2 <- mu1
  
  if (mode == "node") {
    if (k > 0 && is.null(perturbed_set)) perturbed_set <- 1:k
    if (!is.null(perturbed_set) && length(perturbed_set) != k) {
      stop("length(perturbed_set) must equal k")
    }
    
    if (k > 0 && dm != 0) {
      mu2[perturbed_set] <- mu2[perturbed_set] + dm
    }
    
    r <- rep(1, q)
    if (k > 0 && xi != 1) {
      r[perturbed_set] <- sqrt(xi)
    }
    
    D  <- diag(r, q, q)
    S2 <- D %*% S1 %*% D
    affected <- if (k > 0) perturbed_set else integer(0)
  }
  
  if (mode == "edge") {
    if (is.null(edge_change) || length(edge_change) != 2) {
      stop("For mode = 'edge', provide edge_change = c(i,j).")
    }
    
    i <- edge_change[1]
    j <- edge_change[2]
    
    if (i == j) stop("edge_change must contain two distinct nodes.")
    if (any(edge_change < 1) || any(edge_change > q)) {
      stop("edge_change indices must be between 1 and q.")
    }
    
    Theta1 <- solve(S1)
    Theta2 <- Theta1
    
    Theta2[i, j] <- Theta2[i, j] + delta_edge
    Theta2[j, i] <- Theta2[i, j]
    
    if (!is_pd(Theta2)) {
      stop("Theta2 is not positive definite. Try a smaller delta_edge.")
    }
    
    S2 <- solve(Theta2)
    affected <- sort(edge_change)
  }
  
  list(
    mu1 = mu1,
    mu2 = mu2,
    S1 = S1,
    S2 = S2,
    S = affected
  )
}

# ------------------------------------------------------------
# Compute statistics
# ------------------------------------------------------------
compute_stats <- function(Y, grp, max_subset_size = 1) {
  q      <- ncol(Y)
  n1     <- sum(grp == 1)
  n2     <- sum(grp == 2)
  full_W <- LRT(Y, grp)
  full_T <- bartlett_adj(full_W, q, n1, n2)
  
  out <- list(
    full = list(
      W = full_W,
      T = full_T,
      pW = pchisq(full_W, df = q*(q+3)/2, lower.tail = FALSE),
      pT = pchisq(full_T, df = q*(q+3)/2, lower.tail = FALSE)
    ),
    deltas = list(singleton = NULL, pair = NULL, triplet = NULL)
  )
  
  idx <- 1:q
  
  if (max_subset_size >= 1) {
    dW <- sapply(idx, function(j) {
      full_W - LRT(Y[, setdiff(idx, j), drop = FALSE], grp)
    })
    
    dT <- sapply(idx, function(j) {
      W_red <- LRT(Y[, setdiff(idx, j), drop = FALSE], grp)
      T_red <- bartlett_adj(W_red, q - 1, n1, n2)
      full_T - T_red
    })
    
    df1 <- df_delta(1, q)
    
    out$deltas$singleton <- data.frame(
      node = idx,
      dW = dW,
      dT = dT,
      pW = pchisq(dW, df = df1, lower.tail = FALSE),
      pT = pchisq(dT, df = df1, lower.tail = FALSE),
      df = df1
    )
  }
  
  if (max_subset_size >= 2) {
    pairs <- combn(idx, 2, simplify = FALSE)
    
    dW <- vapply(pairs, function(pp) {
      full_W - LRT(Y[, setdiff(idx, pp), drop = FALSE], grp)
    }, 0.0)
    
    dT <- vapply(pairs, function(pp) {
      W_red <- LRT(Y[, setdiff(idx, pp), drop = FALSE], grp)
      T_red <- bartlett_adj(W_red, q - 2, n1, n2)
      full_T - T_red
    }, 0.0)
    
    df2 <- df_delta(2, q)
    
    out$deltas$pair <- data.frame(
      set = sapply(pairs, paste, collapse = ","),
      dW = dW,
      dT = dT,
      pW = pchisq(dW, df = df2, lower.tail = FALSE),
      pT = pchisq(dT, df = df2, lower.tail = FALSE),
      df = df2
    )
  }
  
  if (max_subset_size >= 3) {
    trips <- combn(idx, 3, simplify = FALSE)
    
    dW <- vapply(trips, function(tp) {
      full_W - LRT(Y[, setdiff(idx, tp), drop = FALSE], grp)
    }, 0.0)
    
    dT <- vapply(trips, function(tp) {
      W_red <- LRT(Y[, setdiff(idx, tp), drop = FALSE], grp)
      T_red <- bartlett_adj(W_red, q - 3, n1, n2)
      full_T - T_red
    }, 0.0)
    
    df3 <- df_delta(3, q)
    
    out$deltas$triplet <- data.frame(
      set = sapply(trips, paste, collapse = ","),
      dW = dW,
      dT = dT,
      pW = pchisq(dW, df = df3, lower.tail = FALSE),
      pT = pchisq(dT, df = df3, lower.tail = FALSE),
      df = df3
    )
  }
  
  out
}

# ------------------------------------------------------------
# One replicate
# ------------------------------------------------------------
one_replicate <- function(n, q = 8, k = 2, dm = 1, xi = 0.5,
                          rho = 0.4, A = NULL,
                          standardize_Sigma = TRUE,
                          max_subset_size = 1,
                          perturbed_set = NULL, alpha = 0.05,
                          fixed_nodes = TRUE,
                          mode = c("node", "edge"),
                          edge_change = NULL,
                          delta_edge = 0) {
  
  mode <- match.arg(mode)
  S1 <- make_S1(q = q, rho = rho, A = A,
                standardize = standardize_Sigma)
  
  if (mode == "node") {
    if (k == 0) {
      S_fix <- integer(0)
    } else if (!is.null(perturbed_set)) {
      S_fix <- perturbed_set
    } else if (isTRUE(fixed_nodes)) {
      S_fix <- 1:k
    } else {
      S_fix <- sample(1:q, k)
    }
  } else {
    S_fix <- integer(0)
  }
  
  prm <- build_params(
    q = q, k = k, dm = dm, xi = xi, S1 = S1,
    mu0 = 0, perturbed_set = S_fix,
    mode = mode,
    edge_change = edge_change,
    delta_edge = delta_edge
  )
  
  Y1  <- MASS::mvrnorm(n, prm$mu1, prm$S1)
  Y2  <- MASS::mvrnorm(n, prm$mu2, prm$S2)
  Y   <- rbind(Y1, Y2)
  grp <- rep(1:2, each = n)
  
  stats <- compute_stats(Y, grp, max_subset_size = max_subset_size)
  
  if (mode == "node") {
    is_H0 <- as.integer(dm == 0 && xi == 1)
  } else {
    is_H0 <- as.integer(delta_edge == 0)
  }
  
  sing <- stats$deltas$singleton
  null_nodes <- setdiff(1:q, prm$S)
  
  fwerW <- if (length(null_nodes) == 0) 0L else as.integer(any(sing$pW[null_nodes] < alpha))
  fwerT <- if (length(null_nodes) == 0) 0L else as.integer(any(sing$pT[null_nodes] < alpha))
  
  selW <- which(sing$pW < alpha)
  selT <- which(sing$pT < alpha)
  
  list(
    params = list(
      q = q, n = n, k = k, dm = dm, xi = xi,
      S = prm$S, is_H0 = is_H0,
      mode = mode,
      edge_change = edge_change,
      delta_edge = delta_edge
    ),
    full = stats$full,
    singleton = sing,
    pair = stats$deltas$pair,
    triplet = stats$deltas$triplet,
    sel = list(nodesW = selW, nodesT = selT),
    fwer = list(W = fwerW, T = fwerT)
  )
}

# ------------------------------------------------------------
# Many replicates
# ------------------------------------------------------------
run_sim <- function(Sim = 1000,
                    n_vec = c(10, 50, 100, 250),
                    k_vec = c(0, 1, 2),
                    dm_vec = c(0, 0.5, 1.0, 1.5),
                    xi_vec = c(0.5, 1.0, 1.5),
                    delta_edge_vec = c(0, 0.1, 0.2),
                    q = 8, rho = 0.4,
                    A = NULL,
                    A_seed = 123,
                    standardize_Sigma = TRUE,
                    max_subset_size = 1,
                    alpha = 0.05,
                    fixed_nodes = TRUE,
                    perturbed_set = NULL,
                    mode = c("node", "edge"),
                    edge_change = NULL,
                    seed = 2,
                    n_cores = max(1, parallel::detectCores() - 1)) {
  
  mode <- match.arg(mode)
  
  # Generate the dense covariance factor once.
  if (is.null(A)) {
    A_fixed <- make_A(q = q, A_seed = A_seed)
  } else {
    A_fixed <- A
  }
  
  S1_fixed <- make_S1(q = q, rho = rho, A = A_fixed,
                      standardize = standardize_Sigma)
  Theta1_fixed <- solve(S1_fixed)
  
  message("Baseline covariance S1:")
  print(round(S1_fixed, 3))
  message("Baseline precision Theta1:")
  print(round(Theta1_fixed, 3))
  
  if (mode == "node") {
    grid <- expand.grid(
      n = n_vec, k = k_vec, dm = dm_vec, xi = xi_vec,
      KEEP.OUT.ATTRS = FALSE
    )
  } else {
    grid <- expand.grid(
      n = n_vec, delta_edge = delta_edge_vec,
      KEEP.OUT.ATTRS = FALSE
    )
    grid$k <- 0
    grid$dm <- 0
    grid$xi <- 1
  }
  
  res <- vector("list", nrow(grid))
  
  for (g in seq_len(nrow(grid))) {
    n  <- grid$n[g]
    k  <- grid$k[g]
    dm <- grid$dm[g]
    xi <- grid$xi[g]
    
    if (mode == "edge") {
      delta_edge <- grid$delta_edge[g]
      tag <- sprintf(
        "edge_%d-%d_n=%03d_delta=%.3f",
        edge_change[1], edge_change[2], n, delta_edge
      )
    } else {
      delta_edge <- 0
      tag <- sprintf("node_n=%03d_k=%d_dm=%.3f_xi=%.2f", n, k, dm, xi)
    }
    
    cat("Running:", tag, "\n")
    
    repl <- parallel::mclapply(
      seq_len(Sim),
      function(s) {
        set.seed(seed + 100000 * g + s)
        
        one_replicate(
          n = n, q = q, k = k, dm = dm, xi = xi, rho = rho,
          A = A_fixed,
          standardize_Sigma = standardize_Sigma,
          max_subset_size = max_subset_size,
          alpha = alpha,
          fixed_nodes = fixed_nodes,
          perturbed_set = perturbed_set,
          mode = mode,
          edge_change = edge_change,
          delta_edge = delta_edge
        )
      },
      mc.cores = n_cores
    )
    
    is_H0 <- vapply(repl, function(r) r$params$is_H0, 0L) == 1
    
    rej_full_W <- vapply(repl, function(r) as.integer(r$full$pW < alpha), 0L)
    rej_full_T <- vapply(repl, function(r) as.integer(r$full$pT < alpha), 0L)
    
    fwerW <- vapply(repl, function(r) r$fwer$W, 0L)
    fwerT <- vapply(repl, function(r) r$fwer$T, 0L)
    
    power_any_W <- vapply(repl, function(r) {
      S <- r$params$S
      if (length(S) == 0) return(NA_integer_)
      as.integer(any(r$sel$nodesW %in% S))
    }, 0L)
    
    power_any_T <- vapply(repl, function(r) {
      S <- r$params$S
      if (length(S) == 0) return(NA_integer_)
      as.integer(any(r$sel$nodesT %in% S))
    }, 0L)
    
    power_all_W <- vapply(repl, function(r) {
      S <- r$params$S
      if (length(S) == 0) return(NA_integer_)
      as.integer(all(S %in% r$sel$nodesW))
    }, 0L)
    
    power_all_T <- vapply(repl, function(r) {
      S <- r$params$S
      if (length(S) == 0) return(NA_integer_)
      as.integer(all(S %in% r$sel$nodesT))
    }, 0L)
    
    summary <- list(
      tag = tag,
      grid = grid[g, ],
      N = Sim,
      mode = mode,
      edge_change = edge_change,
      full = list(
        type1_W = if (any(is_H0)) mean(rej_full_W[is_H0]) else NA_real_,
        type1_T = if (any(is_H0)) mean(rej_full_T[is_H0]) else NA_real_,
        power_W = if (any(!is_H0)) mean(rej_full_W[!is_H0]) else NA_real_,
        power_T = if (any(!is_H0)) mean(rej_full_T[!is_H0]) else NA_real_
      ),
      node_level = list(
        fwer_H0_W = if (any(is_H0)) mean(fwerW[is_H0]) else NA_real_,
        fwer_H0_T = if (any(is_H0)) mean(fwerT[is_H0]) else NA_real_,
        power_any_W = if (any(!is_H0)) mean(power_any_W[!is_H0], na.rm = TRUE) else NA_real_,
        power_any_T = if (any(!is_H0)) mean(power_any_T[!is_H0], na.rm = TRUE) else NA_real_,
        power_all_W = if (any(!is_H0)) mean(power_all_W[!is_H0], na.rm = TRUE) else NA_real_,
        power_all_T = if (any(!is_H0)) mean(power_all_T[!is_H0], na.rm = TRUE) else NA_real_
      )
    )
    
    res[[g]] <- list(summary = summary, replicates = repl)
  }
  
  attr(res, "A") <- A_fixed
  attr(res, "S1") <- S1_fixed
  attr(res, "Theta1") <- Theta1_fixed
  res
}

# ------------------------------------------------------------
# Save helpers
# ------------------------------------------------------------
short_name_S <- function(cfg, prefix = "sim") {
  if (cfg$mode == "edge") {
    paste0(
      prefix, "_denseSigma_edge_",
      paste(cfg$edge_change, collapse = "-"),
      ".RData"
    )
  } else {
    s_tag <- if (!is.null(cfg$perturbed_set) && length(cfg$perturbed_set) > 0) {
      paste0("S", paste(cfg$perturbed_set, collapse = "-"))
    } else {
      "H0"
    }
    paste0(prefix, "_denseSigma_", s_tag, ".RData")
  }
}

save_sim_S <- function(out, cfg, dir = "results", prefix = "sim") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fname <- short_name_S(cfg, prefix)
  meta  <- list(
    saved_at = Sys.time(),
    note = "Full configuration is stored in `cfg`",
    file_basename = fname,
    A = attr(out, "A"),
    S1 = attr(out, "S1"),
    Theta1 = attr(out, "Theta1")
  )
  save(out, cfg, meta, file = file.path(dir, fname))
  message("Saved: ", file.path(dir, fname))
  invisible(fname)
}

# ------------------------------------------------------------
# Utility: check admissible precision-edge perturbations
# ------------------------------------------------------------
check_delta_edge <- function(q = 8, rho = 0.4, A = NULL,
                             A_seed = 123,
                             standardize_Sigma = TRUE,
                             edge_change = c(1, 2),
                             delta_grid = seq(-2, 2, by = 0.05)) {
  
  if (is.null(A)) A <- make_A(q = q, A_seed = A_seed)
  
  S1 <- make_S1(q = q, rho = rho, A = A,
                standardize = standardize_Sigma)
  Theta1 <- solve(S1)
  
  out <- sapply(delta_grid, function(delta) {
    Theta2 <- Theta1
    i <- edge_change[1]
    j <- edge_change[2]
    
    Theta2[i, j] <- Theta2[i, j] + delta
    Theta2[j, i] <- Theta2[i, j]
    
    min(eigen(Theta2, symmetric = TRUE, only.values = TRUE)$values)
  })
  
  data.frame(
    delta_edge = delta_grid,
    min_eigen = out,
    positive_definite = out > 1e-8
  )
}

# ------------------------------------------------------------
# Run simulations
# ------------------------------------------------------------
setwd("~/LOOnodeID-070526")

# ------------------------------------------------------------
# Common base configuration for node-change scenarios
# ------------------------------------------------------------
base_cfg <- list(
  mode = "node",
  Sim = 1000,
  n_vec = c(10, 50, 100, 250),
  dm_vec = c(0, 0.5, 1.0, 1.5),
  xi_vec = c(0.5, 1.0, 1.5),
  q = 8,
  rho = 0.4,
  A_seed = 123,
  standardize_Sigma = TRUE,
  max_subset_size = 3,
  alpha = 0.05,
  fixed_nodes = TRUE,
  seed = 123,
  n_cores = max(1, parallel::detectCores() - 1)
)

# ------------------------------------------------------------
# 1) Node-change scenario: S = {1}
# ------------------------------------------------------------
cfgA <- modifyList(
  base_cfg,
  list(k_vec = c(0, 1), perturbed_set = c(1))
)

outA <- run_sim(
  Sim = cfgA$Sim,
  n_vec = cfgA$n_vec,
  k_vec = cfgA$k_vec,
  dm_vec = cfgA$dm_vec,
  xi_vec = cfgA$xi_vec,
  q = cfgA$q,
  rho = cfgA$rho,
  A_seed = cfgA$A_seed,
  standardize_Sigma = cfgA$standardize_Sigma,
  max_subset_size = cfgA$max_subset_size,
  alpha = cfgA$alpha,
  fixed_nodes = cfgA$fixed_nodes,
  perturbed_set = cfgA$perturbed_set,
  mode = cfgA$mode,
  seed = cfgA$seed,
  n_cores = cfgA$n_cores
)

save_sim_S(outA, cfgA, dir = "results")

# ------------------------------------------------------------
# 2) Node-change scenario: S = {1, 2}
# ------------------------------------------------------------
cfgB <- modifyList(
  base_cfg,
  list(k_vec = c(0, 2), perturbed_set = c(1, 2))
)

outB <- run_sim(
  Sim = cfgB$Sim,
  n_vec = cfgB$n_vec,
  k_vec = cfgB$k_vec,
  dm_vec = cfgB$dm_vec,
  xi_vec = cfgB$xi_vec,
  q = cfgB$q,
  rho = cfgB$rho,
  A_seed = cfgB$A_seed,
  standardize_Sigma = cfgB$standardize_Sigma,
  max_subset_size = cfgB$max_subset_size,
  alpha = cfgB$alpha,
  fixed_nodes = cfgB$fixed_nodes,
  perturbed_set = cfgB$perturbed_set,
  mode = cfgB$mode,
  seed = cfgB$seed,
  n_cores = cfgB$n_cores
)

save_sim_S(outB, cfgB, dir = "results")

# ------------------------------------------------------------
# 3) Node-change scenario: S = {1, 2, 3, 4, 5}
# ------------------------------------------------------------
cfgC <- modifyList(
  base_cfg,
  list(k_vec = c(0, 5), perturbed_set = c(1, 2, 3, 4, 5))
)

outC <- run_sim(
  Sim = cfgC$Sim,
  n_vec = cfgC$n_vec,
  k_vec = cfgC$k_vec,
  dm_vec = cfgC$dm_vec,
  xi_vec = cfgC$xi_vec,
  q = cfgC$q,
  rho = cfgC$rho,
  A_seed = cfgC$A_seed,
  standardize_Sigma = cfgC$standardize_Sigma,
  max_subset_size = cfgC$max_subset_size,
  alpha = cfgC$alpha,
  fixed_nodes = cfgC$fixed_nodes,
  perturbed_set = cfgC$perturbed_set,
  mode = cfgC$mode,
  seed = cfgC$seed,
  n_cores = cfgC$n_cores
)

save_sim_S(outC, cfgC, dir = "results")

# ------------------------------------------------------------
# 4) Edge-change scenario: edge = (1, 2)
# ------------------------------------------------------------
edge_cfg <- list(
  mode = "edge",
  Sim = 1000,
  n_vec = c(10, 50, 100, 250),
  delta_edge_vec = c(0, 0.25, 0.75, 1),
  edge_change = c(1, 2),
  q = 8,
  rho = 0.4,
  A_seed = 123,
  standardize_Sigma = TRUE,
  max_subset_size = 3,
  alpha = 0.05,
  seed = 123,
  n_cores = max(1, parallel::detectCores() - 1)
)

A_edge <- make_A(q = edge_cfg$q, A_seed = edge_cfg$A_seed)

S1_tmp <- make_S1(
  q = edge_cfg$q,
  rho = edge_cfg$rho,
  A = A_edge,
  standardize = edge_cfg$standardize_Sigma
)

Theta1_tmp <- solve(S1_tmp)

delta_delete <- -Theta1_tmp[
  edge_cfg$edge_change[1],
  edge_cfg$edge_change[2]
]

edge_cfg$delta_edge_vec <- sort(unique(c(
  edge_cfg$delta_edge_vec,
  delta_delete
)))

print(edge_cfg$delta_edge_vec)

delta_check <- check_delta_edge(
  q = edge_cfg$q,
  rho = edge_cfg$rho,
  A = A_edge,
  standardize_Sigma = edge_cfg$standardize_Sigma,
  edge_change = edge_cfg$edge_change,
  delta_grid = seq(-2, 8, by = 0.05)
)

print(subset(delta_check, positive_definite))

out_edge <- run_sim(
  Sim = edge_cfg$Sim,
  n_vec = edge_cfg$n_vec,
  delta_edge_vec = edge_cfg$delta_edge_vec,
  q = edge_cfg$q,
  rho = edge_cfg$rho,
  A = A_edge,
  A_seed = edge_cfg$A_seed,
  standardize_Sigma = edge_cfg$standardize_Sigma,
  max_subset_size = edge_cfg$max_subset_size,
  alpha = edge_cfg$alpha,
  mode = edge_cfg$mode,
  edge_change = edge_cfg$edge_change,
  seed = edge_cfg$seed,
  n_cores = edge_cfg$n_cores
)

save_sim_S(out_edge, edge_cfg, dir = "results")

