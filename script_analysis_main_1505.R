rm(list = ls())
gc()

## =========================
## Libraries
## =========================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(purrr)
library(cowplot)

## =========================
## Paths
## =========================
setwd("~/LOOnodeID-070526")

base_dir_fig    <- file.path(getwd(), "figures")
base_dir_fig_SM <- file.path(getwd(), "figures_SM")

if (!dir.exists(base_dir_fig)) dir.create(base_dir_fig, recursive = TRUE)
if (!dir.exists(base_dir_fig_SM)) dir.create(base_dir_fig_SM, recursive = TRUE)

## =========================
## Small helpers: labels, file I/O
## =========================
add_n_label <- function(plot, n) {
  plot + ggtitle(bquote(n[1] == n[2] ~ "=" ~ .(n))) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
}

cfg_tag <- function(cfg) {
  if (!is.null(cfg$perturbed_set) && length(cfg$perturbed_set) > 0) {
    paste0("S", paste(cfg$perturbed_set, collapse = "-"))
  } else {
    "H0"
  }
}

load_sim_file <- function(path) {
  e <- new.env()
  load(path, envir = e)
  if (!exists("out", envir = e) || !exists("cfg", envir = e)) {
    stop("File ", path, " does not contain `out` and `cfg`.")
  }
  list(
    out  = e$out,
    cfg  = e$cfg,
    meta = if (exists("meta", envir = e)) e$meta else NULL
  )
}

## =========================
## Compatibility helpers
## =========================
.normalize_singleton <- function(s) {
  if (!"node" %in% names(s)) {
    if ("j" %in% names(s)) {
      s$node <- s$j
    } else if ("var" %in% names(s)) {
      s$node <- s$var
    } else if ("idx" %in% names(s)) {
      s$node <- s$idx
    } else if ("set" %in% names(s)) {
      s$node <- suppressWarnings(as.integer(sub(",.*$", "", as.character(s$set))))
    } else {
      stop(
        "Singleton table missing `node` and no compatible alternative. Columns found: ",
        paste(names(s), collapse = ", ")
      )
    }
  }
  
  s$node <- as.integer(as.character(s$node))
  
  keep <- intersect(c("node", "dW", "dT", "pW", "pT", "df"), names(s))
  s <- s[, keep, drop = FALSE]
  
  s <- s[order(s$node), , drop = FALSE]
  rownames(s) <- NULL
  s
}

.ensure_node <- function(df) {
  if ("node" %in% names(df)) {
    if (is.factor(df$node) || is.character(df$node)) {
      tmp <- gsub("^Node\\s+", "", as.character(df$node))
      if (all(grepl("^[0-9]+$", tmp))) df$node <- as.integer(tmp)
    }
    return(df)
  }
  
  if ("Node" %in% names(df)) {
    names(df)[names(df) == "Node"] <- "node"
    tmp <- gsub("^Node\\s+", "", as.character(df$node))
    if (all(grepl("^[0-9]+$", tmp))) df$node <- as.integer(tmp)
    return(df)
  }
  
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df) && all(grepl("^[0-9]+$", rn))) {
    df$node <- as.integer(rn)
    rownames(df) <- NULL
    return(df)
  }
  
  stop(
    "Expected a `node` column here but could not find or reconstruct it. Available columns: ",
    paste(names(df), collapse = ", ")
  )
}

.ensure_set_label <- function(df) {
  if ("set_label" %in% names(df)) return(df)
  if ("set" %in% names(df)) {
    df$set_label <- paste0("(", as.character(df$set), ")")
    return(df)
  }
  if ("subset" %in% names(df)) {
    df$set_label <- as.character(df$subset)
    return(df)
  }
  df$set_label <- paste0("(", seq_len(nrow(df)), ")")
  df
}

## =========================
## Core accessors
## =========================
get_cell <- function(out, n, k, dm, xi) {
  for (z in out) {
    g <- z$summary$grid
    if (g$n == n && g$k == k && abs(g$dm - dm) < 1e-12 && abs(g$xi - xi) < 1e-12) {
      return(z)
    }
  }
  stop("No matching cell for n=", n, ", k=", k, ", dm=", dm, ", xi=", xi)
}

get_q_from_cell <- function(cell) {
  q <- try(cell$replicates[[1]]$params$q, silent = TRUE)
  if (inherits(q, "try-error") || is.null(q) || is.na(q)) {
    s1 <- .normalize_singleton(cell$replicates[[1]]$singleton)
    q <- max(as.integer(s1$node))
  }
  as.integer(q)
}

## =========================
## H0: singletons tidy + plots
## =========================
extract_H0_tidy <- function(cell) {
  reps <- cell$replicates
  q    <- get_q_from_cell(cell)
  
  full_df <- dplyr::bind_rows(lapply(seq_along(reps), function(i) {
    r <- reps[[i]]$full
    tibble::tibble(sim = i, W = r$W, T = r$T, pW = r$pW, pT = r$pT)
  }))
  
  sing_df <- dplyr::bind_rows(lapply(seq_along(reps), function(i) {
    s <- .normalize_singleton(reps[[i]]$singleton)
    s$sim <- i
    s
  })) %>%
    dplyr::mutate(
      node = as.integer(node),
      node_label = factor(node, levels = 1:q, labels = paste0("Node ", 1:q))
    )
  
  list(full = full_df, singletons = sing_df, q = q)
}

make_plots_H0_for_n <- function(out, n, k = 0, dm = 0, xi = 1, alpha = 0.05) {
  cell <- get_cell(out, n = n, k = k, dm = dm, xi = xi)
  tid  <- extract_H0_tidy(cell)
  q    <- tid$q
  df_full <- q * (q + 3) / 2
  df_d1   <- q + 1
  
  p_full_W <- ggplot(tid$full, aes(x = W)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    stat_function(fun = dchisq, args = list(df = df_full)) +
    labs(x = "W", y = "Density") + theme_classic()
  
  p_full_T <- ggplot(tid$full, aes(x = T)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    stat_function(fun = dchisq, args = list(df = df_full)) +
    labs(x = "T", y = "Density") + theme_classic()
  
  qq_full_W <- ggplot(tid$full, aes(sample = W)) +
    stat_qq(distribution = qchisq, dparams = list(df = df_full)) +
    stat_qq_line(distribution = qchisq, dparams = list(df = df_full)) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()
  
  qq_full_T <- ggplot(tid$full, aes(sample = T)) +
    stat_qq(distribution = qchisq, dparams = list(df = df_full)) +
    stat_qq_line(distribution = qchisq, dparams = list(df = df_full)) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()
  
  p_pval_full <- tid$full %>%
    tidyr::pivot_longer(c(pW, pT), names_to = "stat", values_to = "pval") %>%
    dplyr::mutate(stat = dplyr::recode(stat, pW = "non-Bartlett", pT = "Bartlett")) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "grey90", color = "grey60") +
    stat_function(fun = dunif, args = list(min = 0, max = 1)) +
    facet_wrap(~ stat, ncol = 2) +
    labs(x = "p-value", y = "Density") + theme_classic()
  
  p_delta_W <- ggplot(tid$singletons, aes(x = dW)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    stat_function(fun = dchisq, args = list(df = df_d1), color = "black") +
    facet_wrap(~ node_label, ncol = 4, nrow = 2) +
    labs(x = expression(Delta[j] ~ "(non-Bartlett)"), y = "Density") + theme_classic()
  
  p_delta_T <- ggplot(tid$singletons, aes(x = dT)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    stat_function(fun = dchisq, args = list(df = df_d1), color = "black") +
    facet_wrap(~ node_label, ncol = 4, nrow = 2) +
    labs(x = expression(Delta[j] ~ "(Bartlett)"), y = "Density") + theme_classic()
  
  p_pval_nodes_B <- tid$singletons %>%
    dplyr::transmute(node_label, pT) %>%
    ggplot(aes(x = pT)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "grey90", color = "grey60") +
    stat_function(fun = dunif, args = list(min = 0, max = 1)) +
    facet_wrap(~ node_label, ncol = 4, nrow = 2) +
    labs(x = "p-value (Bartlett)", y = "Density") + theme_classic()
  
  qq_nodes_W <- ggplot(tid$singletons, aes(sample = dW)) +
    stat_qq(distribution = qchisq, dparams = list(df = df_d1)) +
    stat_qq_line(distribution = qchisq, dparams = list(df = df_d1)) +
    facet_wrap(~ node_label, ncol = 4) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()
  
  qq_nodes_T <- ggplot(tid$singletons, aes(sample = dT)) +
    stat_qq(distribution = qchisq, dparams = list(df = df_d1)) +
    stat_qq_line(distribution = qchisq, dparams = list(df = df_d1)) +
    facet_wrap(~ node_label, ncol = 4) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()
  
  type1_by_node <- tid$singletons %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(
      fpr_W = mean(pW < alpha, na.rm = TRUE),
      fpr_T = mean(pT < alpha, na.rm = TRUE),
      .groups = "drop"
    )
  
  delta_long <- tid$singletons %>%
    dplyr::select(dW, dT) %>%
    tidyr::pivot_longer(c(dW, dT), names_to = "which", values_to = "delta") %>%
    dplyr::mutate(stat = dplyr::recode(which, dW = "non-Bartlett", dT = "Bartlett"))
  
  xmax_delta <- max(delta_long$delta, na.rm = TRUE)
  xgrid <- seq(0, xmax_delta * 1.05, length.out = 400)
  overlay_chi <- tidyr::expand_grid(stat = c("non-Bartlett", "Bartlett"), x = xgrid) %>%
    dplyr::mutate(dens = dchisq(x, df = df_d1))
  
  p_delta_both <- ggplot(delta_long, aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    geom_line(data = overlay_chi, aes(x = x, y = dens), inherit.aes = FALSE,
              color = "black", linewidth = 0.8) +
    facet_wrap(~ stat, ncol = 2) +
    labs(x = "Global LRT on G", y = "Density") +
    theme_classic(base_size = 11)
  
  list(
    plots = list(
      full_W = p_full_W, full_T = p_full_T,
      qq_full_W = qq_full_W, qq_full_T = qq_full_T,
      full_pvals = p_pval_full,
      delta_both = p_delta_both,
      delta_W = p_delta_W, delta_T = p_delta_T,
      node_pvals_bartlett = p_pval_nodes_B,
      qq_nodes_W = qq_nodes_W, qq_nodes_T = qq_nodes_T
    ),
    type1 = type1_by_node,
    tidy  = tid,
    df    = list(df_full = df_full, df_d1 = df_d1)
  )
}

## =========================
## H0: pairs/triplets tidy + plots
## =========================
extract_H0_sets <- function(cell, l = 2) {
  stopifnot(l %in% c(2, 3))
  reps <- cell$replicates
  q    <- get_q_from_cell(cell)
  
  df_list <- lapply(seq_along(reps), function(i) {
    z <- if (l == 2) reps[[i]]$pair else reps[[i]]$triplet
    if (is.null(z)) stop("Requested l=", l, " but replicate lacks that component.")
    z$sim <- i
    z
  })
  raw <- dplyr::bind_rows(df_list)
  
  all_sets <- if (l == 2) combn(q, 2, simplify = FALSE) else combn(q, 3, simplify = FALSE)
  all_labels <- sapply(all_sets, function(v) paste0("(", paste(v, collapse = ","), ")"))
  
  raw %>%
    dplyr::mutate(
      set_label = paste0("(", set, ")"),
      set_label = factor(set_label, levels = all_labels)
    )
}

make_plots_H0_for_n_sets <- function(out, n, l = 2, alpha = 0.05) {
  cell <- get_cell(out, n = n, k = 0, dm = 0, xi = 1)
  tid  <- extract_H0_sets(cell, l = l)
  q    <- get_q_from_cell(cell)
  df_chi <- if (l == 2) (2 * q + 1) else (3 * q)
  
  p_delta_T_sets <- ggplot(tid, aes(x = dT)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    stat_function(fun = dchisq, args = list(df = df_chi), color = "black") +
    facet_wrap(~ set_label, ncol = 8, drop = FALSE) +
    labs(x = expression(Delta[M] ~ "(Bartlett)"), y = "Density") + theme_classic()
  
  pvals_sets_B <- ggplot(tid, aes(x = pT)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "grey90", color = "grey60") +
    stat_function(fun = dunif, args = list(min = 0, max = 1)) +
    facet_wrap(~ set_label, ncol = 8, drop = FALSE) +
    labs(x = "p-value (Bartlett)", y = "Density") + theme_classic()
  
  type1_by_set <- tid %>%
    dplyr::group_by(set_label) %>%
    dplyr::summarise(
      fpr_W = mean(pW < alpha, na.rm = TRUE),
      fpr_T = mean(pT < alpha, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(plots = list(deltaT = p_delta_T_sets, pvals = pvals_sets_B), type1 = type1_by_set)
}

## =========================
## Power helpers
## =========================
cell_node_freq <- function(cell, alpha = 0.05, use_bartlett = TRUE) {
  reps <- cell$replicates
  pv_list <- lapply(reps, function(r) {
    s <- .normalize_singleton(r$singleton)
    if (use_bartlett) s$pT else s$pW
  })
  P <- do.call(rbind, pv_list)
  freq <- colMeans(P < alpha, na.rm = TRUE)
  tibble::tibble(node = seq_along(freq), freq = as.numeric(freq))
}

get_targets <- function(out, n, k, xi, dm_vals) {
  for (dm in dm_vals) {
    z <- try(get_cell(out, n, k, dm, xi), silent = TRUE)
    if (!inherits(z, "try-error")) {
      SS <- try(z$replicates[[1]]$params$S, silent = TRUE)
      if (!inherits(SS, "try-error") && length(SS) == k) return(SS)
    }
  }
  1:k
}

power_curve_for_n <- function(out, n, k, xi, dm_vals, alpha = 0.05, use_bartlett = TRUE) {
  targets <- if (k > 0) get_targets(out, n, k, xi, dm_vals) else integer(0)
  pcs <- lapply(dm_vals, function(dm) {
    cell <- try(get_cell(out, n, k, dm, xi), silent = TRUE)
    if (inherits(cell, "try-error") && dm == 0 && abs(xi - 1) < 1e-12) {
      cell <- get_cell(out, n, 0, dm, xi)
    }
    fr <- cell_node_freq(cell, alpha = alpha, use_bartlett = use_bartlett)
    fr$d <- dm
    fr
  })
  dplyr::bind_rows(pcs) %>%
    dplyr::mutate(
      n = n,
      xi = xi,
      group = ifelse(node %in% targets, "Target Nodes", "Other Nodes")
    )
}

## =========================
## H1 singletons helpers
## =========================
extract_H1_singletons <- function(cell, bartlett = TRUE) {
  reps <- cell$replicates
  q    <- get_q_from_cell(cell)
  df <- dplyr::bind_rows(lapply(seq_along(reps), function(i) {
    s <- .normalize_singleton(reps[[i]]$singleton)
    tibble::tibble(
      sim = i,
      node = s$node,
      delta = if (bartlett) s$dT else s$dW
    )
  })) %>%
    dplyr::mutate(
      node = as.integer(node),
      node_label = factor(node, levels = 1:q, labels = paste0("Node ", 1:q))
    )
  list(df = df, q = q)
}

estimate_lambda_nodes <- function(df_delta, df_chi) {
  group_var <- if ("node_label" %in% names(df_delta)) "node_label" else "node"
  df_delta %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(lambda = pmax(mean(delta, na.rm = TRUE) - df_chi, 0), .groups = "drop") %>%
    dplyr::rename(node_label = .data[[group_var]])
}

dens_grid_for_nodes <- function(lambdas, xmax, df_chi) {
  x <- seq(0, xmax, length.out = 400)
  
  if (is.null(lambdas) || nrow(lambdas) == 0 || !all(c("node_label", "lambda") %in% names(lambdas))) {
    return(tibble::tibble(x = x, density = dchisq(x, df = df_chi), node_label = factor(NA_character_)))
  }
  
  dplyr::bind_rows(lapply(seq_len(nrow(lambdas)), function(i) {
    tibble::tibble(
      x = x,
      density = stats::dchisq(x, df = df_chi, ncp = lambdas$lambda[i]),
      node_label = lambdas$node_label[i]
    )
  }))
}

plot_H1_singletons_for_cell <- function(cell, bartlett = TRUE) {
  q       <- get_q_from_cell(cell)
  df_chi  <- q + 1
  ex      <- extract_H1_singletons(cell, bartlett = bartlett)
  ddf     <- ex$df
  S_try   <- try(cell$replicates[[1]]$params$S, silent = TRUE)
  targets <- if (!inherits(S_try, "try-error") && length(S_try)) paste0("Node ", S_try) else character(0)
  
  lambdas <- estimate_lambda_nodes(ddf, df_chi)
  xmax    <- max(ddf$delta, na.rm = TRUE)
  xmax    <- xmax + 0.05 * xmax
  dens    <- dens_grid_for_nodes(lambdas, xmax, df_chi)
  
  ddf$group <- ifelse(ddf$node_label %in% targets, "Target Nodes", "Other Nodes")
  dens$group <- ifelse(!is.na(dens$node_label) & dens$node_label %in% targets, "Target Nodes", "Other Nodes")
  
  p_targets <- ggplot(dplyr::filter(ddf, group == "Target Nodes"), aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "white", color = "grey60") +
    geom_line(data = dplyr::filter(dens, group == "Target Nodes"),
              aes(x = x, y = density), color = "black", linewidth = 0.8) +
    facet_wrap(~ node_label, ncol = min(max(1, length(targets)), 4), scales = "free_y") +
    labs(x = expression(Delta[j] ~ "(Bartlett)"), y = "Density") +
    theme_classic()
  
  p_others <- ggplot(dplyr::filter(ddf, group == "Other Nodes"), aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    geom_line(data = dplyr::filter(dens, group == "Other Nodes"),
              aes(x = x, y = density), color = "black", linewidth = 0.8) +
    facet_wrap(~ node_label, ncol = 4, scales = "free_y") +
    labs(x = expression(Delta[j] ~ "(Bartlett)"), y = "Density") +
    theme_classic()
  
  list(plot = p_targets / p_others, lambdas = lambdas)
}

## =========================
## FWER helpers
## =========================
.compute_fwer_cell <- function(cell, l = 1, alpha = 0.05,
                               method = c("raw", "bonferroni", "holm"),
                               bartlett = TRUE) {
  method <- match.arg(method)
  reps   <- cell$replicates
  pick_tbl <- function(r) {
    if (l == 1) r$singleton else if (l == 2) r$pair else r$triplet
  }
  any_flag <- vapply(reps, function(r) {
    tbl <- pick_tbl(r)
    if (is.null(tbl)) return(FALSE)
    p <- if (bartlett) tbl$pT else tbl$pW
    p <- p[is.finite(p)]
    if (!length(p)) return(FALSE)
    padj <- switch(
      method,
      raw = p,
      bonferroni = p.adjust(p, "bonferroni"),
      holm = p.adjust(p, "holm")
    )
    any(padj < alpha)
  }, logical(1))
  mean(any_flag)
}

.build_fwer_table <- function(out, n_vals, alpha = 0.05) {
  methods <- c("raw", "bonferroni", "holm")
  
  avail_levels_for_n <- function(n) {
    cell <- get_cell(out, n = n, k = 0, dm = 0, xi = 1)
    rr   <- cell$replicates[[1]]
    levs <- c(1L)
    if (!is.null(rr$pair))    levs <- c(levs, 2L)
    if (!is.null(rr$triplet)) levs <- c(levs, 3L)
    levs
  }
  
  rows <- lapply(n_vals, function(n) {
    cell <- get_cell(out, n = n, k = 0, dm = 0, xi = 1)
    lvls <- avail_levels_for_n(n)
    do.call(rbind, lapply(lvls, function(l) {
      do.call(rbind, lapply(methods, function(mth) {
        do.call(rbind, lapply(c(TRUE, FALSE), function(bart) {
          data.frame(
            n = n, l = l, method = mth, bartlett = bart,
            FWER = .compute_fwer_cell(cell, l = l, alpha = alpha, method = mth, bartlett = bart),
            stringsAsFactors = FALSE
          )
        }))
      }))
    }))
  })
  
  out_df <- dplyr::bind_rows(rows)
  out_df$level <- factor(out_df$l, levels = c(1, 2, 3), labels = c("Singletons", "Pairs", "Triplets"))
  out_df
}

## =========================
## Per-config analyzer
## =========================
analyze_config <- function(out, cfg, alpha = 0.05) {
  q        <- cfg$q
  n_vals   <- cfg$n_vec
  dm_vals  <- cfg$dm_vec
  xi_vals  <- cfg$xi_vec
  k_vals   <- setdiff(cfg$k_vec, 0)
  tag      <- cfg_tag(cfg)
  
  results_list <- setNames(vector("list", length(n_vals)), paste0("n=", n_vals))
  for (n in n_vals) {
    results_list[[paste0("n=", n)]] <- make_plots_H0_for_n(out, n = n, k = 0, dm = 0, xi = 1, alpha = alpha)
  }
  
  full_hist_W_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$full_W, n)), ncol = 1)
  full_hist_T_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$full_T, n)), ncol = 1)
  full_pvals_hist_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$full_pvals, n)), ncol = 1)
  qq_full_W_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$qq_full_W, n)), ncol = 1)
  qq_full_T_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$qq_full_T, n)), ncol = 1)
  combined_pval_nodes_bartlett <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$node_pvals_bartlett, n)), ncol = 1)
  combined_deltaT_nodes <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$delta_T, n)), ncol = 1)
  qq_nodes_W_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$qq_nodes_W, n)), ncol = 1)
  qq_nodes_T_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$qq_nodes_T, n)), ncol = 1)
  full_delta_hist_all <- wrap_plots(lapply(n_vals, function(n) add_n_label(results_list[[paste0("n=", n)]]$plots$delta_both, n)), ncol = 1)
  
  type1_wide_W <- bind_rows(lapply(n_vals, function(n) {
    df <- results_list[[paste0("n=", n)]]$type1
    df <- .ensure_node(df)
    df <- dplyr::select(df, node, fpr_W)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_W = round(fpr_W, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_W, names_prefix = "n=") %>%
    dplyr::arrange(node)
  
  type1_wide_T <- bind_rows(lapply(n_vals, function(n) {
    df <- results_list[[paste0("n=", n)]]$type1
    df <- .ensure_node(df)
    df <- dplyr::select(df, node, fpr_T)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_T = round(fpr_T, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_T, names_prefix = "n=") %>%
    dplyr::arrange(node)
  
  df_d1 <- q + 1
  make_qq_overlay_for_n <- function(n, node_cols = 4) {
    sing <- results_list[[paste0("n=", n)]]$tidy$singletons
    sing <- .ensure_node(sing)
    sing <- dplyr::select(sing, node, dW, dT)
    
    long <- sing %>%
      tidyr::pivot_longer(c(dW, dT), names_to = "which", values_to = "val") %>%
      dplyr::mutate(
        stat = ifelse(which == "dT", "Bartlett", "non-Bartlett"),
        node_label = factor(node, levels = sort(unique(node)), labels = paste0("Node ", sort(unique(node))))
      ) %>%
      dplyr::filter(is.finite(val))
    
    qq_pts <- long %>%
      dplyr::group_by(node_label, stat) %>%
      dplyr::arrange(val, .by_group = TRUE) %>%
      dplyr::mutate(
        r = dplyr::row_number(),
        m = dplyr::n(),
        p = (r - 0.5) / m,
        th = qchisq(p, df = df_d1),
        smp = val
      ) %>%
      dplyr::ungroup()
    
    ggplot() +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
      geom_path(
        data = dplyr::filter(qq_pts, stat == "non-Bartlett"),
        aes(x = th, y = smp, linetype = stat, color = stat),
        linewidth = 0.9, show.legend = TRUE
      ) +
      geom_path(
        data = dplyr::filter(qq_pts, stat == "Bartlett"),
        aes(x = th, y = smp, linetype = stat, color = stat),
        linewidth = 0.9, show.legend = TRUE
      ) +
      facet_wrap(~ node_label, ncol = node_cols) +
      scale_linetype_manual(values = c("Bartlett" = "solid", "non-Bartlett" = "solid")) +
      scale_color_manual(values = c("Bartlett" = "black", "non-Bartlett" = "grey60")) +
      guides(linetype = guide_legend(title = NULL), color = guide_legend(title = NULL)) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      theme_classic(base_size = 13) +
      theme(legend.position = "top", strip.text = element_text(face = "plain"))
  }
  
  n_small <- n_vals[n_vals %in% c(10, 50)]
  n_large <- n_vals[n_vals %in% c(100, 250)]
  qq_nodes_overlay_small <- wrap_plots(lapply(n_small, function(n) add_n_label(make_qq_overlay_for_n(n), n)), ncol = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "top")
  qq_nodes_overlay_large <- wrap_plots(lapply(n_large, function(n) add_n_label(make_qq_overlay_for_n(n), n)), ncol = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "top")
  
  H0 <- list(
    per_n = results_list,
    combined = list(
      full_hist_W_all = full_hist_W_all,
      full_hist_T_all = full_hist_T_all,
      full_pvals_hist_all = full_pvals_hist_all,
      full_delta_hist_all = full_delta_hist_all,
      qq_full_W_all = qq_full_W_all,
      qq_full_T_all = qq_full_T_all,
      node_pvals_bartlett_all = combined_pval_nodes_bartlett,
      deltaT_nodes_all = combined_deltaT_nodes,
      qq_nodes_W_all = qq_nodes_W_all,
      qq_nodes_T_all = qq_nodes_T_all,
      qq_overlay_small = qq_nodes_overlay_small,
      qq_overlay_large = qq_nodes_overlay_large
    ),
    type1_singletons_wide = list(W = type1_wide_W, T = type1_wide_T)
  )
  
  pair_plots_by_n <- setNames(lapply(n_vals, function(n) {
    res <- make_plots_H0_for_n_sets(out, n = n, l = 2, alpha = alpha)
    list(deltaT = add_n_label(res$plots$deltaT, n), pvals = add_n_label(res$plots$pvals, n))
  }), paste0("n=", n_vals))
  
  triplet_plots_by_n <- setNames(lapply(n_vals, function(n) {
    res <- make_plots_H0_for_n_sets(out, n = n, l = 3, alpha = alpha)
    list(deltaT = add_n_label(res$plots$deltaT, n), pvals = add_n_label(res$plots$pvals, n))
  }), paste0("n=", n_vals))
  
  type1_pairs_wide_W <- bind_rows(lapply(n_vals, function(n) {
    tmp <- make_plots_H0_for_n_sets(out, n = n, l = 2, alpha = alpha)$type1
    tmp <- .ensure_set_label(tmp)
    df <- dplyr::select(tmp, set_label, fpr_W)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_W = round(fpr_W, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_W, names_prefix = "n=") %>%
    dplyr::arrange(set_label)
  
  type1_pairs_wide_T <- bind_rows(lapply(n_vals, function(n) {
    tmp <- make_plots_H0_for_n_sets(out, n = n, l = 2, alpha = alpha)$type1
    tmp <- .ensure_set_label(tmp)
    df <- dplyr::select(tmp, set_label, fpr_T)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_T = round(fpr_T, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_T, names_prefix = "n=") %>%
    dplyr::arrange(set_label)
  
  type1_triplets_wide_W <- bind_rows(lapply(n_vals, function(n) {
    tmp <- make_plots_H0_for_n_sets(out, n = n, l = 3, alpha = alpha)$type1
    tmp <- .ensure_set_label(tmp)
    df <- dplyr::select(tmp, set_label, fpr_W)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_W = round(fpr_W, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_W, names_prefix = "n=") %>%
    dplyr::arrange(set_label)
  
  type1_triplets_wide_T <- bind_rows(lapply(n_vals, function(n) {
    tmp <- make_plots_H0_for_n_sets(out, n = n, l = 3, alpha = alpha)$type1
    tmp <- .ensure_set_label(tmp)
    df <- dplyr::select(tmp, set_label, fpr_T)
    df$n <- n
    df
  })) %>%
    dplyr::mutate(fpr_T = round(fpr_T, 3)) %>%
    tidyr::pivot_wider(names_from = n, values_from = fpr_T, names_prefix = "n=") %>%
    dplyr::arrange(set_label)
  
  pairs <- list(per_n = pair_plots_by_n, type1_wide = list(W = type1_pairs_wide_W, T = type1_pairs_wide_T))
  triplets <- list(per_n = triplet_plots_by_n, type1_wide = list(W = type1_triplets_wide_W, T = type1_triplets_wide_T))
  
  power <- list(by_k = list())
  for (k_power in k_vals) {
    power_plots_by_xi <- lapply(xi_vals, function(xi) {
      curves_xi <- bind_rows(lapply(n_vals, function(n) {
        power_curve_for_n(out, n, k_power, xi, dm_vals, alpha, use_bartlett = TRUE)
      }))
      plots_row <- lapply(n_vals, function(n) {
        df_n <- filter(curves_xi, n == !!n)
        ggplot() +
          geom_line(data = filter(df_n, group == "Other Nodes"), aes(x = d, y = freq, group = node), color = "grey70", linewidth = 0.8) +
          geom_point(data = filter(df_n, group == "Other Nodes"), aes(x = d, y = freq, group = node), color = "grey70", size = 1.5) +
          geom_line(data = filter(df_n, group == "Target Nodes"), aes(x = d, y = freq, group = node), color = "black", linewidth = 0.8) +
          geom_point(data = filter(df_n, group == "Target Nodes"), aes(x = d, y = freq, group = node), color = "black", size = 1.5) +
          geom_hline(yintercept = alpha, linetype = 2, linewidth = 0.5) +
          scale_x_continuous(breaks = dm_vals) +
          scale_y_continuous(limits = c(0, 1)) +
          xlab(expression(delta[mu])) + ylab("Power") +
          theme_classic(base_size = 12) +
          theme(legend.position = "none", panel.grid.minor = element_blank(),
                axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
      })
      wrap_plots(lapply(seq_along(n_vals), function(i) add_n_label(plots_row[[i]], n_vals[i])), nrow = 1) &
        plot_annotation(title = bquote(xi == .(xi)),
                        theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                                      plot.margin = margin(b = 8)))
    })
    names(power_plots_by_xi) <- paste0("xi=", formatC(xi_vals, format = "f", digits = 1))
    power$by_k[[paste0("k=", k_power)]] <- list(by_xi = power_plots_by_xi)
  }
  
  fwer_tbl <- .build_fwer_table(out, n_vals, alpha = alpha)
  
  mk_fwer_plot <- function(df_level, alpha,
                           shape_unadj = 9, shape_holm = 16,
                           col_bart = "black", col_nonbart = "grey60") {
    df_plot <- dplyr::filter(df_level, method %in% c("raw", "holm")) %>%
      dplyr::mutate(
        Bartlett = ifelse(bartlett, "Bartlett", "non-Bartlett"),
        Method = ifelse(method == "raw", "Unadjusted", "Holm")
      )
    
    ggplot(df_plot,
           aes(x = n, y = FWER,
               group = interaction(Method, Bartlett),
               linetype = Bartlett,
               color = Bartlett,
               shape = Method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.8) +
      geom_hline(yintercept = alpha, linetype = "dotted", color = "black") +
      scale_color_manual(values = c("Bartlett" = col_bart, "non-Bartlett" = col_nonbart)) +
      scale_linetype_manual(values = c("Bartlett" = "solid", "non-Bartlett" = "solid")) +
      scale_shape_manual(values = c("Unadjusted" = shape_unadj, "Holm" = shape_holm)) +
      guides(
        color = "none",
        linetype = guide_legend(override.aes = list(shape = NA, color = c(col_bart, col_nonbart), linetype = c("solid", "solid")), order = 1),
        shape = guide_legend(override.aes = list(linetype = "blank", color = "black"), order = 2)
      ) +
      scale_x_continuous(breaks = sort(unique(df_plot$n))) +
      labs(x = bquote(n[1] == n[2]), y = "FWER", linetype = NULL, shape = NULL) +
      theme_classic(base_size = 12) +
      theme(legend.position = "top")
  }
  
  fwer_plots <- list()
  if (any(fwer_tbl$l == 1)) fwer_plots$singletons <- mk_fwer_plot(dplyr::filter(fwer_tbl, l == 1), alpha)
  if (any(fwer_tbl$l == 2)) fwer_plots$pairs      <- mk_fwer_plot(dplyr::filter(fwer_tbl, l == 2), alpha)
  if (any(fwer_tbl$l == 3)) fwer_plots$triplets   <- mk_fwer_plot(dplyr::filter(fwer_tbl, l == 3), alpha)
  
  FWER <- list(table = fwer_tbl, plots = fwer_plots)
  
  H1_singletons <- list(by_k = list())
  for (k_alt in k_vals) {
    by_xi <- lapply(xi_vals, function(xi) {
      per_n_plots   <- setNames(vector("list", length(n_vals)), paste0("n=", n_vals))
      per_n_lambdas <- setNames(vector("list", length(n_vals)), paste0("n=", n_vals))
      for (n in n_vals) {
        dm_list  <- list()
        lam_list <- list()
        for (dm in setdiff(dm_vals, 0)) {
          cell_try <- try(get_cell(out, n = n, k = k_alt, dm = dm, xi = xi), silent = TRUE)
          if (inherits(cell_try, "try-error")) next
          h1 <- plot_H1_singletons_for_cell(cell_try, bartlett = TRUE)
          p <- add_n_label(h1$plot, n) +
            plot_annotation(title = bquote(delta[mu] == .(dm)),
                            theme = theme(plot.title = element_text(size = 12, hjust = 0.5)))
          dm_list[[paste0("dm=", formatC(dm, format = "f", digits = 2))]] <- p
          lam_list[[paste0("dm=", formatC(dm, format = "f", digits = 2))]] <- h1$lambdas
        }
        per_n_plots[[paste0("n=", n)]] <- dm_list
        per_n_lambdas[[paste0("n=", n)]] <- lam_list
      }
      list(per_n = per_n_plots, lambdas = per_n_lambdas)
    })
    names(by_xi) <- paste0("xi=", formatC(xi_vals, format = "f", digits = 2))
    H1_singletons$by_k[[paste0("k=", k_alt)]] <- by_xi
  }
  
  list(
    label = tag,
    cfg = cfg,
    H0 = H0,
    pairs = pairs,
    triplets = triplets,
    power = power,
    FWER = FWER,
    H1_singletons = H1_singletons
  )
}

## =========================
## Analyze saved simulations
## =========================
analyze_all <- function(paths, alpha = 0.05) {
  res <- list()
  for (path in paths) {
    sim <- load_sim_file(path)
    tag <- cfg_tag(sim$cfg)
    res[[tag]] <- analyze_config(sim$out, sim$cfg, alpha = alpha)
  }
  res
}

save_deltaT_pvals <- function(res, type, n_vals, prefix) {
  for (n in n_vals) {
    deltaT <- res[[type]]$per_n[[paste0("n=", n)]]$deltaT
    pvals  <- res[[type]]$per_n[[paste0("n=", n)]]$pvals
    plot_n <- wrap_plots(deltaT, pvals, ncol = 1, heights = c(1, 1))
    ggsave(
      filename = sprintf("%s_n%d.pdf", prefix, n),
      plot = plot_n,
      path = base_dir_fig_SM,
      width = 10,
      height = 10,
      units = "in",
      dpi = 300,
      device = cairo_pdf
    )
  }
}

xi_label_panel <- function(xi_val) {
  lab <- paste0("xi = ", formatC(xi_val, format = "f", digits = 1))
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = lab, fontface = "bold", size = 4.2) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

.make_xi_tag <- function(x) paste0("xi=", formatC(x, format = "f", digits = 1))

make_power_block <- function(res_scenario, k_power, label_width = 0.15) {
  xi_vals <- res_scenario$cfg$xi_vec
  xi_tags <- vapply(xi_vals, .make_xi_tag, character(1))
  bx <- res_scenario$power$by_k[[paste0("k=", k_power)]]$by_xi
  plots_xi <- lapply(xi_tags, function(tag) {
    p <- bx[[tag]]
    if (is.null(p)) stop("Missing power plot for tag: ", tag, ". Check digits used when naming by_xi.")
    p
  })
  rows <- Map(function(xi, plt) {
    wrap_plots(xi_label_panel(xi), plt, nrow = 1, widths = c(label_width, 1))
  }, xi_vals, plots_xi)
  wrap_plots(rows, ncol = 1)
}



## =========================
## Load dense-Sigma simulations
## =========================
files <- c(
  "results/sim_denseSigma_S1.RData",
  "results/sim_denseSigma_S1-2.RData",
  "results/sim_denseSigma_S1-2-3-4-5.RData"
)

files <- files[file.exists(files)]
if (!length(files)) stop("No dense-Sigma result files found in `results/`.")

RES <- analyze_all(files)

## =========================
## Scenario 1: S1
## =========================
if ("S1" %in% names(RES)) {
  RES$S1$H0$combined$full_hist_W_all
  RES$S1$H0$combined$full_hist_T_all
  RES$S1$H0$combined$full_pvals_hist_all
  RES$S1$H0$combined$full_delta_hist_all
  
  plot_pvals_delta_H0 <- RES$S1$H0$combined$full_delta_hist_all |
    RES$S1$H0$combined$full_pvals_hist_all
  print(plot_pvals_delta_H0)
  
  ggsave(
    filename = "plot_pvals_delta_H0.pdf",
    plot = plot_pvals_delta_H0,
    path = base_dir_fig_SM,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  RES$S1$H0$combined$qq_full_W_all
  RES$S1$H0$combined$qq_full_T_all
  RES$S1$H0$combined$node_pvals_bartlett_all
  RES$S1$H0$combined$deltaT_nodes_all
  
  plot_singletons_H0 <- RES$S1$H0$combined$deltaT_nodes_all |
    RES$S1$H0$combined$node_pvals_bartlett_all
  print(plot_singletons_H0)
  
  ggsave(
    filename = "plot_singletons_H0.pdf",
    plot = plot_singletons_H0,
    path = base_dir_fig_SM,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  RES$S1$H0$combined$qq_nodes_W_all
  RES$S1$H0$combined$qq_nodes_T_all
  RES$S1$H0$combined$qq_overlay_small
  RES$S1$H0$combined$qq_overlay_large
  
  qq_small_n_singletons <- RES$S1$H0$combined$qq_overlay_small
  ggsave(
    filename = "qq_small_n_singletons.pdf",
    plot = qq_small_n_singletons,
    path = base_dir_fig_SM,
    width = 8,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  print(RES$S1$H0$type1_singletons_wide$W)
  print(RES$S1$H0$type1_singletons_wide$T)
  
  n_vals <- RES$S1$cfg$n_vec
  for (n in n_vals) {
    cat("\n== S1: n =", n, "==\n")
    P <- RES$S1$H0$per_n[[paste0("n=", n)]]$plots
    print(P$full_W); print(P$full_T); print(P$full_pvals)
    print(P$qq_full_W); print(P$qq_full_T)
    print(P$delta_W); print(P$delta_T)
    print(P$node_pvals_bartlett)
    print(P$qq_nodes_W); print(P$qq_nodes_T)
  }
  
  for (n in n_vals) {
    cat("\n== S1 pairs: n =", n, "==\n")
    print(RES$S1$pairs$per_n[[paste0("n=", n)]]$deltaT)
    print(RES$S1$pairs$per_n[[paste0("n=", n)]]$pvals)
  }
  
  for (n in n_vals) {
    cat("\n== S1 triplets: n =", n, "==\n")
    print(RES$S1$triplets$per_n[[paste0("n=", n)]]$deltaT)
    print(RES$S1$triplets$per_n[[paste0("n=", n)]]$pvals)
  }
  
  save_deltaT_pvals(RES$S1, "pairs", RES$S1$cfg$n_vec, "plot_pairs")
  save_deltaT_pvals(RES$S1, "triplets", RES$S1$cfg$n_vec, "plot_triplets")
  
  print(RES$S1$pairs$type1_wide$W)
  print(RES$S1$pairs$type1_wide$T)
  print(RES$S1$triplets$type1_wide$W)
  print(RES$S1$triplets$type1_wide$T)
  
  xi_tags <- paste0("xi=", formatC(RES$S1$cfg$xi_vec, format = "f", digits = 1))
  for (x in xi_tags) print(RES$S1$power$by_k[["k=1"]]$by_xi[[x]])
  
  plot_block_power_k1 <- make_power_block(RES$S1, k_power = 1)
  print(plot_block_power_k1)
  
  ggsave(
    filename = "plot_block_power_k1.pdf",
    plot = plot_block_power_k1,
    path = base_dir_fig_SM,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  print(RES$S1$FWER$table)
  print(RES$S1$FWER$plots$singletons)
  print(RES$S1$FWER$plots$pairs)
  print(RES$S1$FWER$plots$triplets)
  
  print(RES$S1$H1_singletons$by_k[["k=1"]][["xi=1.00"]]$per_n[["n=250"]][["dm=1.50"]])
  print(RES$S1$H1_singletons$by_k[["k=1"]][["xi=1.00"]]$lambdas[["n=250"]][["dm=1.50"]])
}

## =========================
## Scenario 2: S1-2
## =========================
if ("S1-2" %in% names(RES)) {
  RES[["S1-2"]]$H0$combined$full_hist_W_all
  RES[["S1-2"]]$H0$combined$full_hist_T_all
  RES[["S1-2"]]$H0$combined$full_pvals_hist_all
  RES[["S1-2"]]$H0$combined$qq_full_W_all
  RES[["S1-2"]]$H0$combined$qq_full_T_all
  RES[["S1-2"]]$H0$combined$node_pvals_bartlett_all
  RES[["S1-2"]]$H0$combined$deltaT_nodes_all
  RES[["S1-2"]]$H0$combined$qq_nodes_W_all
  RES[["S1-2"]]$H0$combined$qq_nodes_T_all
  RES[["S1-2"]]$H0$combined$qq_overlay_small
  RES[["S1-2"]]$H0$combined$qq_overlay_large
  
  print(RES[["S1-2"]]$H0$type1_singletons_wide$W)
  print(RES[["S1-2"]]$H0$type1_singletons_wide$T)
  
  n_vals <- RES[["S1-2"]]$cfg$n_vec
  for (n in n_vals) {
    cat("\n== S1-2: n =", n, "==\n")
    P <- RES[["S1-2"]]$H0$per_n[[paste0("n=", n)]]$plots
    print(P$full_W); print(P$full_T); print(P$full_pvals)
    print(P$qq_full_W); print(P$qq_full_T)
    print(P$delta_W); print(P$delta_T)
    print(P$node_pvals_bartlett)
    print(P$qq_nodes_W); print(P$qq_nodes_T)
  }
  
  for (n in n_vals) {
    cat("\n== S1-2 pairs: n =", n, "==\n")
    print(RES[["S1-2"]]$pairs$per_n[[paste0("n=", n)]]$deltaT)
    print(RES[["S1-2"]]$pairs$per_n[[paste0("n=", n)]]$pvals)
  }
  
  for (n in n_vals) {
    cat("\n== S1-2 triplets: n =", n, "==\n")
    print(RES[["S1-2"]]$triplets$per_n[[paste0("n=", n)]]$deltaT)
    print(RES[["S1-2"]]$triplets$per_n[[paste0("n=", n)]]$pvals)
  }
  
  print(RES[["S1-2"]]$pairs$type1_wide$W)
  print(RES[["S1-2"]]$pairs$type1_wide$T)
  print(RES[["S1-2"]]$triplets$type1_wide$W)
  print(RES[["S1-2"]]$triplets$type1_wide$T)
  
  xi_tags <- paste0("xi=", formatC(RES[["S1-2"]]$cfg$xi_vec, format = "f", digits = 1))
  for (x in xi_tags) print(RES[["S1-2"]]$power$by_k[["k=2"]]$by_xi[[x]])
  
  plot_block_power_k2 <- make_power_block(RES[["S1-2"]], k_power = 2)
  print(plot_block_power_k2)
  
  ggsave(
    filename = "plot_block_power_k2.pdf",
    plot = plot_block_power_k2,
    path = base_dir_fig_SM,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  print(RES[["S1-2"]]$FWER$table)
  print(RES[["S1-2"]]$FWER$plots$singletons)
  print(RES[["S1-2"]]$FWER$plots$pairs)
  print(RES[["S1-2"]]$FWER$plots$triplets)
  
  print(RES[["S1-2"]]$H1_singletons$by_k[["k=2"]][["xi=1.00"]]$per_n[["n=250"]][["dm=1.50"]])
  print(RES[["S1-2"]]$H1_singletons$by_k[["k=2"]][["xi=1.00"]]$lambdas[["n=250"]][["dm=1.50"]])
}

## =========================
## Scenario 3: S1-2-3-4-5
## =========================
if ("S1-2-3-4-5" %in% names(RES)) {
  RES[["S1-2-3-4-5"]]$H0$combined$full_hist_W_all
  RES[["S1-2-3-4-5"]]$H0$combined$full_hist_T_all
  RES[["S1-2-3-4-5"]]$H0$combined$full_pvals_hist_all
  RES[["S1-2-3-4-5"]]$H0$combined$qq_full_W_all
  RES[["S1-2-3-4-5"]]$H0$combined$qq_full_T_all
  RES[["S1-2-3-4-5"]]$H0$combined$node_pvals_bartlett_all
  RES[["S1-2-3-4-5"]]$H0$combined$deltaT_nodes_all
  RES[["S1-2-3-4-5"]]$H0$combined$qq_nodes_W_all
  RES[["S1-2-3-4-5"]]$H0$combined$qq_nodes_T_all
  RES[["S1-2-3-4-5"]]$H0$combined$qq_overlay_small
  RES[["S1-2-3-4-5"]]$H0$combined$qq_overlay_large
  
  print(RES[["S1-2-3-4-5"]]$H0$type1_singletons_wide$W)
  print(RES[["S1-2-3-4-5"]]$H0$type1_singletons_wide$T)
  
  n_vals <- RES[["S1-2-3-4-5"]]$cfg$n_vec
  for (n in n_vals) {
    cat("\n== S1-2-3-4-5: n =", n, "==\n")
    P <- RES[["S1-2-3-4-5"]]$H0$per_n[[paste0("n=", n)]]$plots
    print(P$full_W); print(P$full_T); print(P$full_pvals)
    print(P$qq_full_W); print(P$qq_full_T)
    print(P$delta_W); print(P$delta_T)
    print(P$node_pvals_bartlett)
    print(P$qq_nodes_W); print(P$qq_nodes_T)
  }
  
  for (n in n_vals) {
    cat("\n== S1-2-3-4-5 pairs: n =", n, "==\n")
    print(RES[["S1-2-3-4-5"]]$pairs$per_n[[paste0("n=", n)]]$deltaT)
    print(RES[["S1-2-3-4-5"]]$pairs$per_n[[paste0("n=", n)]]$pvals)
  }
  
  for (n in n_vals) {
    cat("\n== S1-2-3-4-5 triplets: n =", n, "==\n")
    print(RES[["S1-2-3-4-5"]]$triplets$per_n[[paste0("n=", n)]]$deltaT)
    print(RES[["S1-2-3-4-5"]]$triplets$per_n[[paste0("n=", n)]]$pvals)
  }
  
  print(RES[["S1-2-3-4-5"]]$pairs$type1_wide$W)
  print(RES[["S1-2-3-4-5"]]$pairs$type1_wide$T)
  print(RES[["S1-2-3-4-5"]]$triplets$type1_wide$W)
  print(RES[["S1-2-3-4-5"]]$triplets$type1_wide$T)
  
  xi_tags <- paste0("xi=", formatC(RES[["S1-2-3-4-5"]]$cfg$xi_vec, format = "f", digits = 1))
  for (x in xi_tags) print(RES[["S1-2-3-4-5"]]$power$by_k[["k=5"]]$by_xi[[x]])
  
  plot_block_power_k5 <- make_power_block(RES[["S1-2-3-4-5"]], k_power = 5)
  print(plot_block_power_k5)
  
  ggsave(
    filename = "plot_block_power_k5.pdf",
    plot = plot_block_power_k5,
    path = base_dir_fig_SM,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
  
  print(RES[["S1-2-3-4-5"]]$FWER$table)
  print(RES[["S1-2-3-4-5"]]$FWER$plots$singletons)
  print(RES[["S1-2-3-4-5"]]$FWER$plots$pairs)
  print(RES[["S1-2-3-4-5"]]$FWER$plots$triplets)
  
  print(RES[["S1-2-3-4-5"]]$H1_singletons$by_k[["k=5"]][["xi=1.00"]]$per_n[["n=100"]][["dm=1.50"]])
  print(RES[["S1-2-3-4-5"]]$H1_singletons$by_k[["k=5"]][["xi=1.00"]]$lambdas[["n=100"]][["dm=1.50"]])
}

## =========================
## Combined results: S1 table
## =========================
if ("S1" %in% names(RES)) {
  S1 <- RES$S1$H0$type1_singletons_wide
  W1 <- S1$W
  T1 <- S1$T
  
  nodes_chr <- gsub("^Node\\s+", "", as.character(W1$node))
  ord <- order(as.numeric(nodes_chr))
  fmt3 <- function(x) sprintf("%.3f", as.numeric(x))
  
  tab <- tibble::tibble(
    Node = nodes_chr[ord],
    `W (10)`   = fmt3(W1$`n=10`[ord]),
    `T (10)`   = fmt3(T1$`n=10`[ord]),
    `W (50)`   = fmt3(W1$`n=50`[ord]),
    `T (50)`   = fmt3(T1$`n=50`[ord]),
    `W (100)`  = fmt3(W1$`n=100`[ord]),
    `T (100)`  = fmt3(T1$`n=100`[ord]),
    `W (250)`  = fmt3(W1$`n=250`[ord]),
    `T (250)`  = fmt3(T1$`n=250`[ord])
  )
  
  print(
    kableExtra::kbl(
      tab,
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      caption = "Empirical marginal type I error under $H_0$ for singleton tests ($l=1$) in scenario (S1). Entries are rejection frequencies at $\\alpha=0.05$ for the non-Bartlett statistic $\\Delta_j$ based on $W_n$ and its Bartlett-adjusted version based on $T_n$, varying $n_1=n_2\\in\\{10,50,100,250\\}$",
      col.names = c("Node", "$W_n$", "$T_n$", "$W_n$", "$T_n$", "$W_n$", "$T_n$", "$W_n$", "$T_n$"),
      linesep = ""
    ) |>
      kableExtra::add_header_above(c(" " = 1, "$n_1=n_2=10$" = 2, "$n_1=n_2=50$" = 2, "$n_1=n_2=100$" = 2, "$n_1=n_2=250$" = 2)) |>
      kableExtra::kable_styling(latex_options = c("hold_position"), font_size = 7, position = "center")
  )
  
  p1FWER <- RES$S1$FWER$plots$singletons + theme(legend.position = "none")
  p2FWER <- RES$S1$FWER$plots$pairs + theme(legend.position = "none")
  p3FWER <- RES$S1$FWER$plots$triplets + theme(legend.position = "none")
  
  legend <- cowplot::get_legend(RES$S1$FWER$plots$singletons + theme(legend.position = "top"))
  
  hdr <- function(txt) {
    ggplot() + theme_void() +
      ggtitle(txt) +
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
            plot.margin = margin(0, 0, 0, 0))
  }
  
  headers <- patchwork::wrap_plots(hdr("Singletons"), hdr("Pairs"), hdr("Triplets"), ncol = 3)
  plots <- patchwork::wrap_plots(p1FWER, p2FWER, p3FWER, ncol = 3)
  body <- patchwork::wrap_plots(headers, plots, ncol = 1, heights = c(0.05, 1))
  plot_combinedFWER <- cowplot::plot_grid(legend, body, ncol = 1, rel_heights = c(0.05, 1))
  print(plot_combinedFWER)
  
  ggsave(
    filename = "plot_combinedFWER.pdf",
    plot = plot_combinedFWER,
    path = base_dir_fig,
    width = 8,
    height = 4,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
}

## =========================
## Aggregate selected singleton/pair/triplet plots
## =========================
.collect_for_n <- function(out, n, j, pair_ab, trip_abc) {
  cell <- get_cell(out, n = n, k = 0, dm = 0, xi = 1)
  q    <- get_q_from_cell(cell)
  
  sing_tbl <- extract_H0_tidy(cell)$singletons
  sing_tbl <- .ensure_node(sing_tbl)
  
  sing <- sing_tbl %>%
    dplyr::filter(node == j) %>%
    dplyr::transmute(subset = paste0("Singleton ", j), pval = pT, delta = dT)
  sing$n <- n
  
  pr <- extract_H0_sets(cell, l = 2)
  lab2 <- paste0("(", paste(pair_ab, collapse = ","), ")")
  pairs <- pr %>%
    dplyr::filter(set_label == lab2) %>%
    dplyr::transmute(subset = paste0("Pair ", lab2), pval = pT, delta = dT)
  pairs$n <- n
  
  tr <- extract_H0_sets(cell, l = 3)
  lab3 <- paste0("(", paste(trip_abc, collapse = ","), ")")
  trips <- tr %>%
    dplyr::filter(set_label == lab3) %>%
    dplyr::transmute(subset = paste0("Triplet ", lab3), pval = pT, delta = dT)
  trips$n <- n
  
  list(
    df_p = dplyr::bind_rows(sing, pairs, trips),
    df_d = dplyr::bind_rows(sing, pairs, trips) %>% dplyr::select(n, subset, delta),
    q = q
  )
}

make_aggregate_plots <- function(files, scenario = "S1",
                                 singleton = 1, pair = c(1, 2), triplet = c(1, 2, 3),
                                 ns_keep = NULL,
                                 bins_p = 20, bins_d = 30) {
  scenario_file_map <- function(paths) {
    nm <- vapply(paths, function(p) {
      e <- load_sim_file(p)
      cfg_tag(e$cfg)
    }, character(1))
    stats::setNames(paths, nm)
  }
  
  fmap <- scenario_file_map(files)
  if (!scenario %in% names(fmap)) stop("Scenario not found in files: ", scenario)
  sim <- load_sim_file(fmap[[scenario]])
  out <- sim$out
  cfg <- sim$cfg
  n_all <- cfg$n_vec
  if (!is.null(ns_keep)) n_all <- intersect(n_all, ns_keep)
  n_all <- sort(unique(n_all))
  
  coll <- purrr::map(n_all, ~ .collect_for_n(out, n = .x, j = singleton, pair_ab = pair, trip_abc = triplet))
  df_p <- dplyr::bind_rows(purrr::map(coll, "df_p"))
  df_d <- dplyr::bind_rows(purrr::map(coll, "df_d"))
  q <- purrr::pluck(coll, 1, "q")
  
  df_map <- tibble::tibble(
    subset = c(
      paste0("Singleton ", singleton),
      paste0("Pair (", paste(pair, collapse = ","), ")"),
      paste0("Triplet (", paste(triplet, collapse = ","), ")")
    ),
    df_chi = c(q + 1, 2 * q + 1, 3 * q)
  )
  
  overlay_chi <- df_d %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(xmax = max(delta, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(df_map, by = "subset") %>%
    dplyr::mutate(xmax = ifelse(is.finite(xmax) & xmax > 0, xmax, 10 + df_chi)) %>%
    dplyr::group_by(subset, df_chi) %>%
    dplyr::do({
      xr <- seq(0, .$xmax[1] * 1.05, length.out = 400)
      tibble::tibble(x = xr, dens = dchisq(xr, df = .$df_chi[1]))
    }) %>%
    dplyr::ungroup()
  
  df_p$subset <- factor(df_p$subset, levels = df_map$subset)
  df_d$subset <- factor(df_d$subset, levels = df_map$subset)
  overlay_chi$subset <- factor(overlay_chi$subset, levels = df_map$subset)
  
  pvals_agg <- ggplot2::ggplot(df_p, ggplot2::aes(x = pval)) +
    ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), bins = bins_p, fill = "grey90", color = "grey60") +
    ggplot2::stat_function(fun = dunif, args = list(min = 0, max = 1), linewidth = 0.8, color = "black") +
    ggplot2::facet_grid(rows = vars(subset), cols = vars(n), scales = "free_y",
                        labeller = ggplot2::label_bquote(cols = bold(n[1] == n[2] ~ "=" ~ .(n)))) +
    ggplot2::labs(x = "p-value (Bartlett)", y = "Density") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  
  delta_agg <- ggplot2::ggplot(df_d, ggplot2::aes(x = delta)) +
    ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)), bins = bins_d, fill = "grey90", color = "grey60") +
    ggplot2::geom_line(data = overlay_chi, ggplot2::aes(x = x, y = dens), inherit.aes = FALSE, color = "black", linewidth = 0.8) +
    ggplot2::facet_grid(rows = vars(subset), cols = vars(n), scales = "free_y",
                        labeller = ggplot2::label_bquote(cols = bold(n[1] == n[2] ~ "=" ~ .(n)))) +
    ggplot2::labs(x = expression(Delta[M] ~ "(Bartlett)"), y = "Density") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  
  list(pvals = pvals_agg, delta = delta_agg)
}

if ("S1" %in% names(RES)) {
  agg_S1 <- make_aggregate_plots(
    files,
    scenario = "S1",
    singleton = 1,
    pair = c(1, 2),
    triplet = c(1, 2, 3),
    ns_keep = c(10, 50, 100, 250)
  )
  
  plot_pvals_1_12_123 <- agg_S1$pvals
  plot_deltaM_1_12_123 <- agg_S1$delta
  print(plot_deltaM_1_12_123)
  
  ggsave(
    filename = "plot_deltaM_1_12_123.pdf",
    plot = plot_deltaM_1_12_123,
    path = base_dir_fig,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
}

## =========================
## Selected power comparison across S1/S2/S3
## =========================
n_focus <- c(50, 100)
xi_focus <- c(1.00, 0.50)
alpha <- 0.05

scenarios <- list(
  S1 = list(path = "results/sim_denseSigma_S1.RData", k = 1),
  S2 = list(path = "results/sim_denseSigma_S1-2.RData", k = 2),
  S3 = list(path = "results/sim_denseSigma_S1-2-3-4-5.RData", k = 5)
)
scenarios <- scenarios[vapply(scenarios, function(x) file.exists(x$path), logical(1))]

axis_theme <- function(show_y = TRUE) {
  theme(
    legend.position = "none",
    axis.title.y = if (show_y) element_text(size = 9) else element_blank(),
    axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(2, 2, 2, 2),
    text = element_text(size = 9)
  )
}

make_power_panel <- function(out, n, k_power, xi, dm_vals, alpha = 0.05, show_y = TRUE) {
  df_nxi <- power_curve_for_n(out, n = n, k = k_power, xi = xi, dm_vals = dm_vals, alpha = alpha, use_bartlett = TRUE)
  ggplot() +
    geom_line(data = dplyr::filter(df_nxi, group == "Other Nodes"), aes(x = d, y = freq, group = node), linewidth = 0.8, color = "grey70") +
    geom_point(data = dplyr::filter(df_nxi, group == "Other Nodes"), aes(x = d, y = freq, group = node), size = 1.5, color = "grey70") +
    geom_line(data = dplyr::filter(df_nxi, group == "Target Nodes"), aes(x = d, y = freq, group = node), linewidth = 0.8, color = "black") +
    geom_point(data = dplyr::filter(df_nxi, group == "Target Nodes"), aes(x = d, y = freq, group = node), size = 1.5, color = "black") +
    geom_hline(yintercept = alpha, linetype = 2, linewidth = 0.5) +
    scale_x_continuous(breaks = dm_vals, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(delta[mu]), y = if (show_y) "Power" else NULL) +
    ggtitle(bquote(bold(n[1] == n[2] ~ "=" ~ .(n)))) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    axis_theme(show_y = show_y)
}

build_block_for_scenario <- function(path, k_power, n_focus, xi_focus, alpha = 0.05, show_y = TRUE) {
  sim <- load_sim_file(path)
  out <- sim$out
  cfg <- sim$cfg
  dm_vals <- cfg$dm_vec
  n_avail <- intersect(n_focus, cfg$n_vec)
  if (!length(n_avail)) stop("Requested n not available in: ", path)
  
  rows <- lapply(seq_along(xi_focus), function(r) {
    xi_val <- xi_focus[r]
    cols <- lapply(seq_along(n_avail), function(j) {
      show_y_col <- show_y && j == 1
      make_power_panel(out, n = n_avail[j], k_power = k_power, xi = xi_val,
                       dm_vals = dm_vals, alpha = alpha, show_y = show_y_col)
    })
    wrap_plots(cols, nrow = 1)
  })
  wrap_plots(rows, ncol = 1)
}

if (all(c("S1", "S2", "S3") %in% names(scenarios))) {
  hdr_S1 <- ggplot() + theme_void() + ggtitle("(S1)") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  hdr_S2 <- ggplot() + theme_void() + ggtitle("(S2)") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  hdr_S3 <- ggplot() + theme_void() + ggtitle("(S3)") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  block_S1 <- build_block_for_scenario(scenarios$S1$path, scenarios$S1$k, n_focus, xi_focus, alpha, show_y = TRUE)
  block_S2 <- build_block_for_scenario(scenarios$S2$path, scenarios$S2$k, n_focus, xi_focus, alpha, show_y = FALSE)
  block_S3 <- build_block_for_scenario(scenarios$S3$path, scenarios$S3$k, n_focus, xi_focus, alpha, show_y = FALSE)
  
  xi_row1 <- xi_label_panel(xi_focus[1])
  xi_row2 <- xi_label_panel(xi_focus[2])
  
  design <- "
            ABCD
            EFGH
            IFGH
            "
  
  plot_power_sel <- wrap_plots(
    plot_spacer(), hdr_S1, hdr_S2, hdr_S3,
    xi_row1, block_S1, block_S2, block_S3,
    xi_row2,
    design = design
  ) + plot_layout(widths = c(0.30, 1, 1, 1), heights = c(0.02, 1, 1))
  
  print(plot_power_sel)
  
  ggsave(
    filename = "plot_power_sel.pdf",
    plot = plot_power_sel,
    path = base_dir_fig,
    width = 14,
    height = 8,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
}

## =========================
## H1 example plot for S1-2
## =========================
if (file.exists("results/sim_denseSigma_S1-2.RData")) {
  sim <- load_sim_file("results/sim_denseSigma_S1-2.RData")
  out <- sim$out
  cell <- get_cell(out, n = 100, k = 2, dm = 1.50, xi = 0.50)
  q <- get_q_from_cell(cell)
  ex <- extract_H1_singletons(cell, bartlett = TRUE)
  ddf <- ex$df
  df_chi <- q + 1
  S_try <- try(cell$replicates[[1]]$params$S, silent = TRUE)
  lamb <- estimate_lambda_nodes(ddf, df_chi)
  
  if (!inherits(S_try, "try-error") && length(S_try) == 2) {
    targets_labels <- paste0("Node ", S_try)
  } else {
    targets_labels <- as.character(head(lamb$node_label[order(lamb$lambda, decreasing = TRUE)], 2))
  }
  
  xmax <- max(ddf$delta, na.rm = TRUE)
  xmax <- xmax + 0.05 * xmax
  dens <- dens_grid_for_nodes(lamb, xmax, df_chi)
  ddf$group <- ifelse(ddf$node_label %in% targets_labels, "Target", "Other")
  dens$group <- ifelse(dens$node_label %in% targets_labels, "Target", "Other")
  
  p_left <- ggplot(dplyr::filter(ddf, group == "Target"), aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "white", color = "grey60") +
    geom_line(data = dplyr::filter(dens, group == "Target"), aes(x = x, y = density), color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    facet_wrap(~ node_label, ncol = 1, scales = "free_y") +
    labs(x = expression(Delta[j] ~ "(Bartlett)"), y = "Density") +
    theme_classic() +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  
  p_right <- ggplot(dplyr::filter(ddf, group == "Other"), aes(x = delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey90", color = "grey60") +
    geom_line(data = dplyr::filter(dens, group == "Other"), aes(x = x, y = density), color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    facet_wrap(~ node_label, ncol = 2, scales = "free_y") +
    labs(x = expression(Delta[j] ~ "(Bartlett)"), y = "Density") +
    theme_classic() +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  
  plot_H1_S12_n100 <- p_left | p_right
  plot_H1_S12_n100 <- plot_H1_S12_n100 + plot_layout(widths = c(1, 2))
  print(plot_H1_S12_n100)
  
  ggsave(
    filename = "plot_H1_S12_n100.pdf",
    plot = plot_H1_S12_n100,
    path = base_dir_fig,
    width = 10,
    height = 6,
    units = "in",
    dpi = 300,
    device = cairo_pdf
  )
}


## =========================
## Edge-change power curves: edge = (1, 2)
## =========================

edge_file <- "results/sim_denseSigma_edge_1-2.RData"

if (!file.exists(edge_file)) {
  stop("Edge file not found: ", edge_file)
}

edge_sim <- load_sim_file(edge_file)
out_edge <- edge_sim$out
edge_cfg <- edge_sim$cfg

edge_power_holm <- do.call(rbind, lapply(out_edge, function(cell) {
  reps <- cell$replicates
  g    <- cell$summary$grid
  
  vals <- sapply(reps, function(r) {
    p_adj <- p.adjust(r$singleton$pT, method = "holm")
    sel   <- which(p_adj < edge_cfg$alpha)
    endpoints <- r$params$S
    
    if (g$delta_edge == 0) {
      any_det <- as.integer(length(sel) > 0)
      all_det <- as.integer(length(sel) > 0)
    } else {
      any_det <- as.integer(any(sel %in% endpoints))
      all_det <- as.integer(all(endpoints %in% sel))
    }
    
    c(any = any_det, all = all_det)
  })
  
  data.frame(
    n = g$n,
    delta_edge = g$delta_edge,
    prob_any = mean(vals["any", ], na.rm = TRUE),
    prob_all = mean(vals["all", ], na.rm = TRUE)
  )
}))

print(edge_power_holm)

edge_power_long <- edge_power_holm |>
  tidyr::pivot_longer(
    cols = c(prob_any, prob_all),
    names_to = "criterion",
    values_to = "prob"
  ) |>
  dplyr::mutate(
    criterion = dplyr::recode(
      criterion,
      prob_any = "At least one of the 2 nodes",
      prob_all = "Both nodes"
    )
  )

make_edge_power_panel_holm <- function(df, n_val) {
  ggplot(
    df[df$n == n_val, ],
    aes(
      x = delta_edge,
      y = prob,
      linetype = criterion,
      shape = criterion
    )
  ) +
    geom_line(color = "black", linewidth = 0.9) +
    geom_point(color = "black", size = 1.8) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(
      breaks = sort(unique(df$delta_edge)),
      labels = function(x) sprintf("%.2f", x)
    ) +
    geom_hline(yintercept = edge_cfg$alpha, linetype = "dotted") +
    labs(
      x = expression(delta[edge]),
      y = "Pr (detect perturbed edge)",
      linetype = NULL,
      shape = NULL
    ) +
    theme_classic(base_size = 12) +
    ggtitle(bquote(n[1] == n[2] ~ "=" ~ .(n_val))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "top"
    )
}

edge_plot_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 10),
    panel.spacing = unit(0.8, "lines"),
    axis.text.x = element_text(size = 8)
  )

x_breaks_edge <- sort(unique(edge_power_long$delta_edge))

## =========================
## Plot 1: edge detection power
## =========================

plot_edge_power_holm <- ggplot(
  edge_power_long,
  aes(
    x = delta_edge,
    y = prob,
    linetype = criterion,
    shape = criterion,
    group = criterion
  )
) +
  geom_line(color = "black", linewidth = 0.9) +
  geom_point(color = "black", size = 1.8) +
  facet_wrap(
    ~ n,
    nrow = 1,
    labeller = label_bquote(n[1] == n[2] ~ "=" ~ .(n))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = x_breaks_edge,
    labels = function(x) sprintf("%.2f", x)
  ) +
  geom_hline(yintercept = edge_cfg$alpha, linetype = "dotted") +
  labs(
    x = expression(delta[edge]),
    y = "Pr (detect perturbed edge)",
    linetype = NULL,
    shape = NULL
  ) +
  edge_plot_theme

print(plot_edge_power_holm)

ggsave(
  filename = "appendix_plot_edge_change_denseSigma_1_and_2.pdf",
  plot     = plot_edge_power_holm,
  path     = base_dir_fig_SM,
  width    = 11,
  height   = 4,
  units    = "in",
  dpi      = 300,
  device   = cairo_pdf
)


## =========================
## Plot 2: node-wise rejection frequencies
## =========================

edge_node_freq <- edge_node_freq |>
  dplyr::mutate(
    node_type = dplyr::case_when(
      node == 1 ~ "Node 1",
      node == 2 ~ "Node 2",
      TRUE      ~ "Other nodes"
    ),
    node_type = factor(node_type, levels = c("Node 1", "Node 2", "Other nodes"))
  )

plot_edge_node_freq <- ggplot(
  edge_node_freq,
  aes(
    x = delta_edge,
    y = rejection_freq_holm,
    group = node,
    color = node_type,
    linewidth = node_type
  )
) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ n,
    nrow = 1,
    labeller = label_bquote(n[1] == n[2] ~ "=" ~ .(n))
  ) +
  scale_color_manual(values = c(
    "Node 1"      = "black",
    "Node 2"      = "darkblue",
    "Other nodes" = "grey70"
  )) +
  scale_linewidth_manual(values = c(
    "Node 1"      = 1.0,
    "Node 2"      = 1.0,
    "Other nodes" = 0.6
  )) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = sort(unique(edge_node_freq$delta_edge)),
    labels = function(x) sprintf("%.2f", x)
  ) +
  labs(
    x = expression(delta[edge]),
    y = "Node-wise rejection frequency",
    color = NULL
  ) +
  guides(linewidth = "none") +
  edge_plot_theme

print(plot_edge_node_freq)

ggsave(
  filename = "appendix_plot_edge_nodewise_rejection_freq_denseSigma_1-2.pdf",
  plot = plot_edge_node_freq,
  path = base_dir_fig_SM,
  width = 11,
  height = 4,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)