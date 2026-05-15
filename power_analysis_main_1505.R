rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)

setwd("~/LOOnodeID-070526")

base_dir_fig_SM <- file.path(getwd(), "figures_SM")
if (!dir.exists(base_dir_fig_SM)) dir.create(base_dir_fig_SM, recursive = TRUE)

files <- list(
  S1           = "results/sim_denseSigma_S1.RData",
  `S1-2`       = "results/sim_denseSigma_S1-2.RData",
  `S1-2-3-4-5` = "results/sim_denseSigma_S1-2-3-4-5.RData"
)

alpha    <- 0.05
n_focus  <- c(50, 100)
xi_focus <- c(1.00, 0.50)

load_sim_file <- function(path) {
  e <- new.env()
  load(path, envir = e)
  if (!exists("out", envir = e) || !exists("cfg", envir = e)) {
    stop("File ", path, " does not contain `out` and `cfg`.")
  }
  list(out = e$out, cfg = e$cfg)
}

get_cell <- function(out, n, k, dm, xi) {
  for (z in out) {
    g <- z$summary$grid
    if (g$n == n && g$k == k && abs(g$dm - dm) < 1e-12 && abs(g$xi - xi) < 1e-12) return(z)
  }
  stop("No matching cell.")
}

get_q_from_cell <- function(cell) {
  q <- try(cell$replicates[[1]]$params$q, silent = TRUE)
  if (is.numeric(q) && length(q) == 1 && is.finite(q)) return(as.integer(q))
  nrow(cell$replicates[[1]]$singleton)
}

.pick_pcol <- function(s, prefer = c("pT", "pW")) {
  for (nm in prefer) if (nm %in% names(s)) return(nm)
  stop("No p-value column found.")
}

get_targets <- function(out, n, k, xi, dm_vals) {
  for (dm in dm_vals) {
    z <- try(get_cell(out, n, k, dm, xi), silent = TRUE)
    if (!inherits(z, "try-error")) {
      S <- try(z$replicates[[1]]$params$S, silent = TRUE)
      if (!inherits(S, "try-error") && length(S) == k) return(as.integer(S))
    }
  }
  seq_len(k)
}

cell_rejection_matrix <- function(cell, alpha = 0.05,
                                  which_p = "pT",
                                  adj_method = c("none", "holm", "bonferroni")) {
  adj_method <- match.arg(adj_method)
  q <- get_q_from_cell(cell)
  R <- matrix(FALSE, nrow = length(cell$replicates), ncol = q)
  
  for (i in seq_along(cell$replicates)) {
    s <- cell$replicates[[i]]$singleton
    pcol <- .pick_pcol(s, if (which_p == "pT") c("pT", "pW") else c("pW", "pT"))
    p <- as.numeric(s[[pcol]])
    
    if ("node" %in% names(s)) {
      tmp <- rep(NA_real_, q)
      idx <- as.integer(s$node)
      tmp[idx] <- p
      p <- tmp
    }
    
    if (adj_method != "none") {
      ok <- is.finite(p)
      p[ok] <- p.adjust(p[ok], method = adj_method)
    }
    
    R[i, ] <- is.finite(p) & p < alpha
  }
  R
}

any_discovery_curve_for_n <- function(out, n, k, xi, dm_vals,
                                      alpha = 0.05,
                                      adj_method = "none") {
  targets <- get_targets(out, n, k, xi, dm_vals)
  
  bind_rows(lapply(dm_vals, function(dm) {
    cell <- get_cell(out, n, k, dm, xi)
    R <- cell_rejection_matrix(cell, alpha = alpha, which_p = "pT", adj_method = adj_method)
    tibble(
      d = dm,
      prob = mean(rowSums(R[, targets, drop = FALSE]) >= 1)
    )
  })) |>
    mutate(n = n, xi = xi)
}

make_any_panel <- function(df, dm_vals) {
  ggplot(df, aes(x = d, y = prob)) +
    geom_line(color = "black", linewidth = 0.9) +
    geom_point(color = "black", size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    scale_x_continuous(breaks = dm_vals) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(delta[mu]),
      y = expression(Pr~"(reject at least one false null)")
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))
}

add_n_label <- function(p, n) {
  p + ggtitle(bquote(n[1] == n[2] ~ "=" ~ .(n))) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
}

build_any_block_core_for_file <- function(path, k_power, xi_vals, n_vals,
                                          alpha = 0.05,
                                          adj_method = "none") {
  sim <- load_sim_file(path)
  out <- sim$out
  cfg <- sim$cfg
  
  rows <- lapply(xi_vals, function(xi) {
    panels <- lapply(n_vals, function(n) {
      df <- any_discovery_curve_for_n(out, n, k_power, xi, cfg$dm_vec, alpha, adj_method)
      add_n_label(make_any_panel(df, cfg$dm_vec), n)
    })
    wrap_plots(panels, nrow = 1)
  })
  
  wrap_plots(rows, ncol = 1)
}

make_xi_column <- function(xi_vals) {
  labs <- lapply(xi_vals, function(xi) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("bold(xi == ", formatC(xi, format = "f", digits = 1), ")"),
               parse = TRUE, size = 4.2) +
      theme_void()
  })
  wrap_plots(labs, ncol = 1)
}

add_center_title <- function(p, title_txt) {
  ggdraw() +
    draw_plot(p, 0, 0, 1, 1) +
    draw_label(title_txt, x = 0.55, y = 1, hjust = 0.5,
               vjust = 0.5, fontface = "bold", size = 13)
}

make_any_final <- function(adj_method, filename) {
  pS1 <- build_any_block_core_for_file(files$S1, 1, xi_focus, n_focus, alpha, adj_method)
  pS2 <- build_any_block_core_for_file(files$`S1-2`, 2, xi_focus, n_focus, alpha, adj_method)
  pS3 <- build_any_block_core_for_file(files$`S1-2-3-4-5`, 5, xi_focus, n_focus, alpha, adj_method)
  
  final <- wrap_plots(
    make_xi_column(xi_focus),
    add_center_title(pS1, "(S1)"),
    add_center_title(pS2, "(S2)"),
    add_center_title(pS3, "(S3)"),
    nrow = 1,
    widths = c(0.18, 1, 1, 1)
  )
  
  ggsave(filename, final, path = base_dir_fig_SM,
         width = 14, height = 6, units = "in", dpi = 300, device = cairo_pdf)
  
  final
}

plot_any_none <- make_any_final("none", "plot_H1_curves_at_least_1_rej_none.pdf")
plot_any_holm <- make_any_final("holm", "plot_H1_curves_at_least_1_rej_holm.pdf")


#################################################
#################################################

all_discovery_curve_for_n <- function(out, n, k, xi, dm_vals,
                                      alpha = 0.05,
                                      adj_method = "none") {
  targets <- get_targets(out, n, k, xi, dm_vals)
  
  bind_rows(lapply(dm_vals, function(dm) {
    cell <- get_cell(out, n, k, dm, xi)
    R <- cell_rejection_matrix(cell, alpha = alpha, which_p = "pT", adj_method = adj_method)
    tibble(
      d = dm,
      prob = mean(rowSums(R[, targets, drop = FALSE]) == length(targets))
    )
  })) |>
    mutate(n = n, xi = xi)
}

make_all_panel <- function(df, dm_vals) {
  ggplot(df, aes(x = d, y = prob)) +
    geom_line(color = "black", linewidth = 0.9) +
    geom_point(color = "black", size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    scale_x_continuous(breaks = dm_vals) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(delta[mu]),
      y = expression(Pr~"(all false nulls rejected)")
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))
}

build_all_block_core_for_file <- function(path, k_power, xi_vals, n_vals,
                                          alpha = 0.05,
                                          adj_method = "none") {
  sim <- load_sim_file(path)
  out <- sim$out
  cfg <- sim$cfg
  
  rows <- lapply(xi_vals, function(xi) {
    panels <- lapply(n_vals, function(n) {
      df <- all_discovery_curve_for_n(out, n, k_power, xi, cfg$dm_vec, alpha, adj_method)
      add_n_label(make_all_panel(df, cfg$dm_vec), n)
    })
    wrap_plots(panels, nrow = 1)
  })
  
  wrap_plots(rows, ncol = 1)
}

make_all_final <- function(adj_method, filename) {
  pS1 <- build_all_block_core_for_file(files$S1, 1, xi_focus, n_focus, alpha, adj_method)
  pS2 <- build_all_block_core_for_file(files$`S1-2`, 2, xi_focus, n_focus, alpha, adj_method)
  pS3 <- build_all_block_core_for_file(files$`S1-2-3-4-5`, 5, xi_focus, n_focus, alpha, adj_method)
  
  final <- wrap_plots(
    make_xi_column(xi_focus),
    add_center_title(pS1, "(S1)"),
    add_center_title(pS2, "(S2)"),
    add_center_title(pS3, "(S3)"),
    nrow = 1,
    widths = c(0.18, 1, 1, 1)
  )
  
  ggsave(filename, final, path = base_dir_fig_SM,
         width = 14, height = 6, units = "in", dpi = 300, device = cairo_pdf)
  
  final
}

plot_all_none <- make_all_final("none", "plot_H1_curves_all_rejected_none.pdf")
plot_all_holm <- make_all_final("holm", "plot_H1_curves_all_rejected_holm.pdf")


###############################################
###############################################

## =========================
## Appendix figures: any/all discovery, all n, xi = 0.5, 1.0, 1.5
## =========================

xi_order_appendix <- c(0.5, 1.0, 1.5)

make_appendix_fig <- function(file_path, k_power, out_filename,
                              type = c("any", "all"),
                              adj_method = "holm",
                              alpha = 0.05) {
  type <- match.arg(type)
  
  sim <- load_sim_file(file_path)
  out <- sim$out
  cfg <- sim$cfg
  
  xi_vals_use <- xi_order_appendix[xi_order_appendix %in% cfg$xi_vec]
  n_vals_use  <- cfg$n_vec
  
  panels <- lapply(xi_vals_use, function(xi) {
    row_plots <- lapply(n_vals_use, function(n) {
      if (type == "any") {
        df <- any_discovery_curve_for_n(
          out, n, k_power, xi, cfg$dm_vec,
          alpha = alpha,
          adj_method = adj_method
        )
        p <- make_any_panel(df, cfg$dm_vec)
      } else {
        df <- all_discovery_curve_for_n(
          out, n, k_power, xi, cfg$dm_vec,
          alpha = alpha,
          adj_method = adj_method
        )
        p <- make_all_panel(df, cfg$dm_vec)
      }
      add_n_label(p, n)
    })
    
    xi_lab <- ggplot() +
      annotate(
        "text", x = 0.5, y = 0.5,
        label = paste0("bold(xi == ", formatC(xi, format = "f", digits = 1), ")"),
        parse = TRUE, size = 4.2
      ) +
      theme_void()
    
    wrap_plots(xi_lab, wrap_plots(row_plots, nrow = 1),
               nrow = 1, widths = c(0.18, 1))
  })
  
  p_block <- wrap_plots(panels, ncol = 1)
  
  ggsave(
    filename = out_filename,
    plot     = p_block,
    path     = base_dir_fig_SM,
    width    = 14,
    height   = 9,
    units    = "in",
    dpi      = 300,
    device   = cairo_pdf
  )
  
  p_block
}

## S1
make_appendix_fig(files$S1, 1, "appendix_any_S1_alln_holm.pdf", type = "any")
make_appendix_fig(files$S1, 1, "appendix_all_S1_alln_holm.pdf", type = "all")

## S2
make_appendix_fig(files$`S1-2`, 2, "appendix_any_S2_alln_holm.pdf", type = "any")
make_appendix_fig(files$`S1-2`, 2, "appendix_all_S2_alln_holm.pdf", type = "all")

## S3
make_appendix_fig(files$`S1-2-3-4-5`, 5, "appendix_any_S3_alln_holm.pdf", type = "any")
make_appendix_fig(files$`S1-2-3-4-5`, 5, "appendix_all_S3_alln_holm.pdf", type = "all")
