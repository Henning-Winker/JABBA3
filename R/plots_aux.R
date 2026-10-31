# Plots for catch-only assessments, prior predictive runs and auxiliary data

#' Prior predictive against posterior
#'
#' Compares a prior predictive run (\code{jabba3(..., prior.only = TRUE)}: catch and priors only, the data likelihoods
#' switched off) with a fit to data: the trajectories (median and 95\% interval) and the marginal densities of selected
#' parameters and of terminal-year status. Shows how much the data (indices, auxiliary series) move the assessment away
#' from what the priors and the catch alone imply (the prior predictive model of Kell et al. 2025). Both runs are best fitted
#' with NUTS; fits without draws use multivariate normal draws at the mode.
#' @param prior the prior predictive run (or any reference run), or a named list of two runs (e.g.
#'   \code{list(Prior = fp, Posterior = fd)}), in which case \code{posterior} is ignored
#' @param posterior the fit to data
#' @param type trajectories to show, any of "B/Bmsy", "F/Fmsy", "B/K", "B", "F"
#' @param pars marginal densities: any of "r", "K", "MSY", "Bmsy", "Fmsy", "bkinit", "sdb" and the terminal-year "B/Bmsy",
#'   "F/Fmsy", "B/K"; \code{NULL} for trajectories only
#' @param labels names of the two runs (default from the list names, or "Prior" and "Posterior")
#' @param ndraws draws used for fits without NUTS draws
#' @return a patchwork (trajectories above, densities below)
#' @examples
#' \dontrun{
#' fp <- jabba3(d0, pr0, jb3_model(shape = "fox"), prior.only = TRUE)
#' fd <- jabba3(d0, pr0, jb3_model(shape = "fox"))
#' jbplot_priorpost_trj(fp, fd)
#' jbplot_trj(list(Prior = fp, Posterior = fd), type = c("B/Bmsy", "F/Fmsy"))   # trajectories only
#' }
#' @export
jbplot_priorpost_trj <- function(prior, posterior = NULL, type = c("B/Bmsy", "F/Fmsy"), pars = c("r", "K", "B/Bmsy"),
                                 labels = NULL, ndraws = 1000){
  fl <- if(inherits(prior, "jabba3")) list(prior, posterior) else prior
  if(length(fl) != 2 || !all(vapply(fl, inherits, TRUE, "jabba3"))) stop("need two jabba3 fits", call. = FALSE)
  labels <- labels %||% names(fl) %||% c("Prior", "Posterior"); names(fl) <- labels
  cols <- stats::setNames(c("grey55", "dodgerblue3"), labels)
  p1 <- jbplot_trj(fl, type = type, refpts = TRUE) + ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) + ggplot2::theme(legend.position = "bottom")
  if(!length(pars)) return(p1)
  dd <- do.call(rbind, lapply(labels, function(n){
    f <- fl[[n]]; dr <- .get_draws(f, ndraws); ly <- as.character(max(f$years))
    do.call(rbind, lapply(pars, function(p){
      v <- if(p %in% c("B/Bmsy", "F/Fmsy", "B/K")) dr[[p]][, ly] else dr$pars[[p]]
      if(is.null(v)) return(NULL)
      data.frame(run = n, par = if(p %in% c("B/Bmsy", "F/Fmsy", "B/K")) paste(p, ly) else p, value = v)
    }))
  }))
  dd$run <- factor(dd$run, labels); dd$par <- factor(dd$par, unique(dd$par))
  logp <- c("K", "MSY", "Bmsy")
  dd$value[dd$par %in% logp] <- log10(dd$value[dd$par %in% logp])
  levels(dd$par)[levels(dd$par) %in% logp] <- paste0("log10 ", levels(dd$par)[levels(dd$par) %in% logp])
  p2 <- ggplot2::ggplot(dd, ggplot2::aes(.data$value, fill = .data$run, colour = .data$run)) +
    ggplot2::geom_density(alpha = 0.35, adjust = 1.2) + ggplot2::facet_wrap(~par, scales = "free", nrow = 1) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) + ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = NULL, y = "Density") + theme_jabba3() + ggplot2::theme(legend.position = "none")
  if(!requireNamespace("patchwork", quietly = TRUE)) return(list(p1, p2))
  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 0.8))
}

#' Fits to the auxiliary data
#'
#' Observed auxiliary series (effort, F, F/Fmsy, B/Bmsy, B/K, Z) with the fitted values, their 95\% interval and the 95\%
#' predictive interval, one facet per series; several fits are coloured by run.
#' @inheritParams jbplot_fits
#' @param title optional plot title (e.g. to flag hypothetical data)
#' @export
jbplot_auxfits <- function(x, log = FALSE, ppd = TRUE, colour = "darkorange3", title = NULL){
  fl <- .fitlist(x)
  d <- do.call(rbind, lapply(names(fl), function(n){ a <- fl[[n]]$auxiliary; if(is.null(a)) return(NULL); a$run <- n; a }))
  if(is.null(d)) stop("no auxiliary data in the fit(s)", call. = FALSE)
  d$lab <- factor(paste0(d$series, " (", d$type, ")"), unique(paste0(d$series, " (", d$type, ")")))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year))
  if(length(unique(d$run)) == 1){
    if(ppd) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ppd_lci, ymax = .data$ppd_uci), fill = "grey85")
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = colour, alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = .data$hat), colour = colour, linewidth = 0.7)
  } else {
    cols <- stats::setNames(jb3cols(length(fl)), names(fl))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$hat, colour = .data$run), linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = cols, name = NULL)
  }
  o <- d[d$run == unique(d$run)[1], ]
  p <- p + ggplot2::geom_point(data = o, ggplot2::aes(y = .data$obs), shape = 21, fill = "white", size = 1.8) +
    ggplot2::facet_wrap(~lab, scales = "free_y") + scale_x_year() + ggplot2::labs(x = "Year", y = "Auxiliary data", title = title) +
    theme_jabba3()
  if(log) p + ggplot2::scale_y_log10() else p + ggplot2::expand_limits(y = 0)
}
