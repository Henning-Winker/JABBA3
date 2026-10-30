# Further JABBA graphics in ggplot2: posterior predictive checks, standardised residuals, index overview, biplot,
# catch error, status priors, ensembles, catch risk, summary and Selex status

.fills_kobe <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")

# draws of a fit, generating MVN draws at the mode when the fit has none
.get_draws <- function(fit, ndraws = 1000) fit$draws %||% jb3_draws(fit, ndraws)

#' Posterior predictive check
#'
#' As JABBA's \code{jbplot_PPC()}: for each draw, the discrepancy D = sum(((log I - log I.hat)/sd)^2) of the observed
#' indices against that of indices replicated from the observation model; points on either side of the 1:1 line, and the
#' Bayesian p-value P(D.rep > D.obs) (near 0.5 under a well-specified model, extreme values < 0.05 or > 0.95 flag
#' misfit). Uses the NUTS draws, or MVN draws at the mode (generated if the fit has none).
#' @param fit a \code{jabba3} fit
#' @param joint pool all indices into one check
#' @param ndraws number of MVN draws when the fit has none
#' @param npoints maximum points per panel
#' @return a ggplot; the Bayesian p-values in attribute "ppc"
#' @export
jbplot_ppc <- function(fit, joint = FALSE, ndraws = 1000, npoints = 1000){
  dr <- .get_draws(fit, ndraws); ix <- fit$index
  y <- log(ix$obs); mu <- dr$Ihat; sd <- dr$sdI; nd <- nrow(mu)   # Imean is on the log scale
  grp <- if(joint) rep("Combined", nrow(ix)) else ix$index
  set.seed(1)
  yr <- mu + matrix(stats::rnorm(length(mu)), nd)*sd
  d <- do.call(rbind, lapply(unique(grp), function(g){ w <- which(grp == g)
    Dx <- rowSums(((matrix(y[w], nd, length(w), byrow = TRUE) - mu[, w, drop = FALSE])/sd[, w, drop = FALSE])^2)
    Dy <- rowSums(((yr[, w, drop = FALSE] - mu[, w, drop = FALSE])/sd[, w, drop = FALSE])^2)
    data.frame(index = g, Dx = Dx, Dy = Dy) }))
  bp <- do.call(rbind, lapply(split(d, factor(d$index, unique(d$index))), function(x)
    data.frame(index = x$index[1], Bayesian.p = mean(x$Dy > x$Dx), nobs = sum(grp == x$index[1]))))
  lab <- stats::setNames(sprintf("%s: p = %.3f", bp$index, bp$Bayesian.p), bp$index)
  dp <- do.call(rbind, lapply(split(d, d$index), function(x) x[seq_len(min(nrow(x), npoints)), ]))
  p <- ggplot2::ggplot(dp, ggplot2::aes(.data$Dx, .data$Dy)) +
    ggplot2::geom_point(shape = 21, colour = "grey40", fill = "grey70", alpha = 0.6, size = 1.4) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.7) +
    ggplot2::facet_wrap(~index, scales = "free", labeller = ggplot2::as_labeller(lab)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(x = "Discrepancy, observed data", y = "Discrepancy, replicated data",
                  subtitle = paste("Posterior predictive check (", dr$type %||% "draws", ")")) + theme_jabba3()
  attr(p, "ppc") <- bp
  p
}

#' Standardised residuals
#'
#' As JABBA's \code{jbplot_stdresiduals()}: standardised log residuals (residual / total observation SD) of all indices
#' as segments from zero, a boxplot per year, a loess smoother through the yearly means, and SDNR with its 95\%
#' critical value (chi-square; SDNR above it suggests the observation error is underestimated).
#' @param fit a \code{jabba3} fit
#' @export
jbplot_stdresiduals <- function(fit){
  d <- fit$index[is.finite(fit$index$std.residual), ]; ni <- length(unique(d$index)); n <- nrow(d)
  d$index <- factor(d$index, unique(d$index))
  off <- stats::setNames(seq(-0.2, 0.2, length.out = max(ni, 2))[seq_len(ni)], levels(d$index))
  d$x <- d$year + if(ni > 1) off[as.character(d$index)] else 0
  m <- stats::aggregate(list(r = d$std.residual), list(year = d$year), mean)
  m$fit <- if(nrow(m) > 4) stats::predict(stats::loess(r ~ year, m)) else m$r
  sdnr <- sqrt(sum(d$std.residual^2)/(n - 1)); crit <- sqrt(stats::qchisq(0.95, n - 1)/(n - 1))
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$std.residual)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$year), fill = "grey90", colour = "grey60", outlier.shape = NA, width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$std.residual, colour = .data$index)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$x, fill = .data$index), shape = 21, size = 1.8) +
    ggplot2::geom_line(data = m, ggplot2::aes(.data$year, .data$fit), linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) +
    ggplot2::scale_fill_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Standardised residual", subtitle = sprintf("SDNR = %.2f (critical value %.2f)", sdnr, crit)) +
    theme_jabba3()
}

#' Abundance indices with observation error
#'
#' As JABBA's \code{jbplot_indices()}: each index scaled to its mean, with the 95\% interval from the input log-scale
#' standard errors (dodged by index).
#' @param data a \code{jb3data} object or a fit
#' @export
jbplot_indices <- function(data){
  if(inherits(data, "jabba3")) data <- data$data
  I <- data$I; mu <- colMeans(I, na.rm = TRUE); SE <- data$SE; if(is.null(SE)) SE <- I*NA
  d <- do.call(rbind, lapply(colnames(I), function(n) data.frame(year = data$years, index = n, value = I[, n]/mu[n],
                                                                se = SE[, n])))
  d <- d[is.finite(d$value), ]; d$se[!is.finite(d$se)] <- 0
  d$lci <- d$value*exp(-1.96*d$se); d$uci <- d$value*exp(1.96*d$se)
  ni <- length(unique(d$index)); d$index <- factor(d$index, colnames(I)[colnames(I) %in% d$index])
  pd <- ggplot2::position_dodge(width = 0.5)
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value, colour = .data$index, fill = .data$index)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), width = 0.4, position = pd, colour = "grey40") +
    ggplot2::geom_line(position = pd, alpha = 0.6) +
    ggplot2::geom_point(position = pd, shape = 21, colour = "grey30", size = 1.8) +
    ggplot2::scale_colour_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) +
    ggplot2::scale_fill_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Index (scaled to mean 1)", subtitle = paste(data$stock, "- 95% intervals from the input SE")) +
    theme_jabba3()
}

# 2-d Gaussian kernel density on a grid and the density levels enclosing given probability masses
.kde2 <- function(x, y, n = 120, lims = NULL){
  lims <- lims %||% c(range(x), range(y))
  bw <- function(v) 1.06*min(stats::sd(v), stats::IQR(v)/1.34)*length(v)^(-1/5)
  hx <- max(bw(x), 1e-6); hy <- max(bw(y), 1e-6)
  gx <- seq(lims[1], lims[2], length.out = n); gy <- seq(lims[3], lims[4], length.out = n)
  ax <- outer(gx, x, function(a, b) stats::dnorm(a, b, hx)); ay <- outer(gy, y, function(a, b) stats::dnorm(a, b, hy))
  z <- tcrossprod(ax, ay)/length(x)
  list(x = gx, y = gy, z = z)
}
.hpd_levels <- function(z, probs){
  v <- sort(as.vector(z), decreasing = TRUE); cs <- cumsum(v)/sum(v)
  vapply(probs, function(p) v[which(cs >= p)[1]], 0)
}

#' Stock status biplot
#'
#' The 'post-modern' biplot of Quinn and Collie (2005), as JABBA's \code{jbplot_biplot()}: B/Bmsy against F/Fmsy with
#' the zones of a hockey-stick harvest control rule (F target \code{ftarget} x Fmsy above Bmsy, decreasing to zero at
#' \code{bthreshold} x Bmsy; green below the target line, red above F = Fmsy line scaled the same way, yellow in between),
#' the trajectory, and 50/80/95\% highest-density contours of the terminal status with the share in each zone.
#' @param fit a \code{jabba3} fit
#' @param ftarget target F/Fmsy above Bmsy
#' @param bthreshold B/Bmsy at which the target F reaches zero
#' @param npoints draws for the terminal status when the fit has none
#' @export
jbplot_biplot <- function(fit, ftarget = 0.8, bthreshold = 0.2, npoints = 2000){
  N <- length(fit$years)
  trj <- data.frame(year = fit$years, stock = fit$ts$est[fit$ts$quantity == "B/Bmsy"], harvest = fit$ts$est[fit$ts$quantity == "F/Fmsy"])
  pt <- .kobe_cloud(fit, N, npoints)
  fmax <- max(2, trj$harvest, stats::quantile(pt$harvest, 0.98), na.rm = TRUE)*1.05
  bmax <- max(2, trj$stock, stats::quantile(pt$stock, 0.98), na.rm = TRUE)*1.05
  zone <- function(b, f) ifelse(b > 1, ifelse(f < ftarget, "green", ifelse(f > 1, "red", "yellow")),
                                ifelse(b > bthreshold + (1 - bthreshold)/ftarget*f, "green",
                                       ifelse(b < bthreshold + (1 - bthreshold)*f, "red", "yellow")))
  g <- expand.grid(f = seq(0, fmax, length.out = 300), b = seq(0, bmax, length.out = 300)); g$zone <- zone(g$b, g$f)
  pz <- round(100*prop.table(table(factor(zone(pt$stock, pt$harvest), c("green", "yellow", "red")))), 1)
  k <- .kde2(pt$harvest, pt$stock, lims = c(0, fmax, 0, bmax)); lv <- .hpd_levels(k$z, c(0.5, 0.8, 0.95))
  cfill <- c("0.95" = "cornsilk4", "0.8" = "grey", "0.5" = "cornsilk2")
  cl <- do.call(rbind, lapply(seq_along(lv), function(i){ L <- grDevices::contourLines(k$x, k$y, k$z, levels = lv[i])
    if(!length(L)) return(NULL)
    do.call(rbind, lapply(seq_along(L), function(j) data.frame(level = names(cfill)[4 - i], grp = paste(i, j), f = L[[j]]$x, b = L[[j]]$y))) }))
  cl$level <- factor(cl$level, names(cfill)); cl <- cl[order(cl$level), ]
  sel <- .selyears(trj, "harvest"); sel$lab <- factor(sel$year, sel$year)
  ggplot2::ggplot() +
    ggplot2::geom_raster(data = g, ggplot2::aes(.data$f, .data$b, fill = .data$zone), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = .fills_kobe, breaks = c("green", "yellow", "red"),
                               labels = sprintf("%s%%", pz[c("green", "yellow", "red")]), name = paste("Status", max(fit$years))) +
    lapply(names(cfill), function(l) if(any(cl$level == l)) ggplot2::geom_polygon(data = cl[cl$level == l, ],
      ggplot2::aes(.data$f, .data$b, group = .data$grp), fill = cfill[[l]], colour = NA, alpha = 0.85)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3) + ggplot2::geom_vline(xintercept = 1, linetype = 3) +
    ggplot2::geom_path(data = trj, ggplot2::aes(.data$harvest, .data$stock), na.rm = TRUE) +
    ggplot2::geom_point(data = trj, ggplot2::aes(.data$harvest, .data$stock), size = 1, na.rm = TRUE) +
    ggplot2::geom_point(data = sel, ggplot2::aes(.data$harvest, .data$stock, shape = .data$lab), fill = "white", size = 3.2) +
    ggplot2::scale_shape_manual(values = c(22, 21, 24), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, fmax), ylim = c(0, bmax), expand = FALSE) +
    ggplot2::labs(x = expression(F/F[MSY]), y = expression(B/B[MSY]),
                  subtitle = sprintf("%s - HCR zones (F target %.1f Fmsy, threshold %.1f Bmsy); contours 50/80/95%%", fit$run, ftarget, bthreshold)) +
    theme_jabba3() + ggplot2::theme(legend.position = "right", legend.title = ggplot2::element_text(size = 8))
}

#' Catch error: observed and estimated catch
#'
#' As JABBA's \code{jbplot_catcherror()}: observed catch (points) and estimated catch with its 95\% interval, and the
#' log ratio of estimated to observed catch. Without catch error the catch is assumed known and the plot says so.
#' @param fit a \code{jabba3} fit
#' @export
jbplot_catcherror <- function(fit){
  d <- fit$ts[fit$ts$quantity == "Catch", ]; d$obs <- fit$data$catch[seq_len(nrow(d))]
  d$lr <- log(d$est/d$obs)
  ce <- isTRUE(fit$model$catch.error)
  p1 <- ggplot2::ggplot(d, ggplot2::aes(.data$year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey75") +
    ggplot2::geom_line(ggplot2::aes(y = .data$est, colour = "Estimated"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = .data$obs, colour = "Observed"), shape = 21, fill = "white", size = 1.8) +
    ggplot2::scale_colour_manual(values = c(Estimated = "blue", Observed = "black")) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = NULL, y = "Catch", subtitle = if(ce) sprintf("catch error, CV %.2f", fit$model$sdc) else "catch error off: catch assumed known") +
    theme_jabba3() + ggplot2::theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1))
  if(!ce || !requireNamespace("patchwork", quietly = TRUE)) return(p1)
  p2 <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$lr)) + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$year, yend = 0), colour = "grey40") + ggplot2::geom_point(size = 1.2) +
    scale_x_year() + ggplot2::labs(x = "Year", y = "log(est/obs)") + theme_jabba3(9)
  p1 / p2 + patchwork::plot_layout(heights = c(3, 1))
}

#' Status priors: prior against posterior
#'
#' As JABBA's \code{jbplot_bprior()}, for each status prior given in \code{jb3_priors(status = )}: the lognormal prior on
#' B/K, B/Bmsy or F/Fmsy in its year against the posterior (draws, or the lognormal approximation at the mode), with the
#' posterior-to-prior mean (PPMR) and variance ratios (PPVR, CV based as in JABBA).
#' @param fit a \code{jabba3} fit
#' @return a ggplot, or NULL (with a message) when the fit has no status priors
#' @export
jbplot_bprior <- function(fit){
  st <- fit$priors$status
  if(is.null(st) || !nrow(st)){ message("no status priors in this fit (jb3_priors(status = ...))"); return(invisible(NULL)) }
  qmap <- c(BK = "B/K", BBmsy = "B/Bmsy", FFmsy = "F/Fmsy")
  d <- do.call(rbind, lapply(seq_len(nrow(st)), function(i){
    q <- qmap[[as.character(st$type[i])]]; yi <- match(st$year[i], fit$years)
    t <- fit$ts[fit$ts$quantity == q & fit$ts$year == st$year[i], ]
    post <- if(!is.null(fit$draws[[q]])) fit$draws[[q]][, yi] else exp(stats::rnorm(5000, log(t$est), t$se.log))
    mu <- log(st$median[i]); s <- st$log.sd[i]
    rng <- range(c(stats::qlnorm(c(0.001, 0.995), mu, s), stats::quantile(post, c(0.001, 0.999), na.rm = TRUE)))
    x <- seq(max(rng[1], 1e-6), rng[2], length.out = 400)
    dd <- stats::density(post[is.finite(post)], adjust = 2)
    cvp <- stats::sd(post)/mean(post); cvr <- sqrt(exp(s^2) - 1)
    lab <- sprintf("%s %s: PPMR %.2f, PPVR %.2f", q, st$year[i], mean(post)/exp(mu + s^2/2), cvp^2/cvr^2)
    rbind(data.frame(panel = lab, x = x, y = stats::dlnorm(x, mu, s), type = "Prior"),
          data.frame(panel = lab, x = dd$x, y = dd$y, type = "Posterior"))
  }))
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y, fill = .data$type, colour = .data$type)) +
    ggplot2::geom_area(position = "identity", alpha = 0.4) +
    ggplot2::scale_fill_manual(values = c(Prior = "grey40", Posterior = "dodgerblue")) +
    ggplot2::scale_colour_manual(values = c(Prior = "grey30", Posterior = "dodgerblue4")) +
    ggplot2::facet_wrap(~panel, scales = "free") + ggplot2::labs(x = NULL, y = "Density") + theme_jabba3() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), legend.position = "bottom")
}

# trajectory draws of one quantity: from the fit's draws, else lognormal per year from se.log
.trj_draws <- function(fit, q, n = 1000, seed = 1){
  if(!is.null(fit$draws[[q]])){ x <- fit$draws[[q]]; return(x[sample.int(nrow(x), n, replace = nrow(x) < n), , drop = FALSE]) }
  t <- fit$ts[fit$ts$quantity == q, ]; set.seed(seed)
  s <- ifelse(is.finite(t$se.log), t$se.log, 0)
  exp(matrix(stats::rnorm(n*nrow(t)), n)*matrix(s, n, nrow(t), byrow = TRUE) + matrix(log(t$est), n, nrow(t), byrow = TRUE))
}

#' Ensemble of fits
#'
#' Compares several fits (scenarios, sensitivity runs), as JABBA's \code{jbplot_ensemble()}: trajectories coloured by
#' run with their 95\% intervals and, on the right, the Kobe plot with each run's terminal-year uncertainty cloud and its
#' median. With \code{joint = TRUE} the runs are pooled (equal weight) into one ensemble trajectory with its interval and
#' one Kobe cloud.
#' @param x a list of \code{jabba3} fits (named), or a workflow fit
#' @param type quantities for the trajectory panels
#' @param joint pool the runs into one ensemble
#' @param kobe add the Kobe panel (needs patchwork); FALSE returns the trajectory panels only
#' @param ndraws draws per run (from the fit's draws, or lognormal at the mode)
#' @export
jbplot_ensemble <- function(x, type = c("B/Bmsy", "F/Fmsy", "B/K", "Catch"), joint = FALSE, kobe = TRUE, ndraws = 1000){
  if(inherits(x, "jb3ensemble") && joint){                       # weighted ensemble: its own (weighted) draws
    p1 <- jbplot_trj(x, type = type, ncol = 2) +
      ggplot2::labs(subtitle = paste0("Weighted ensemble of ", length(x$fits), " runs (", paste(sprintf("%s %.2f", names(x$weights), x$weights), collapse = ", "), ")"))
    if(!kobe || !requireNamespace("patchwork", quietly = TRUE)) return(p1)
    return(p1 + .kobe_panel(list(Ensemble = x), ndraws) + patchwork::plot_layout(widths = c(2, 1.2)))
  }
  fl <- .fitlist(x, steps = TRUE); nr <- length(fl)
  if(!joint) p1 <- jbplot_trj(fl, type = type, ncol = 2) else {
    d <- do.call(rbind, lapply(type, function(q){
      m <- do.call(rbind, lapply(fl, .trj_draws, q = q, n = ndraws))
      qs <- apply(m, 2, stats::quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
      data.frame(year = fl[[1]]$years, quantity = q, est = qs[2, ], lci = qs[1, ], uci = qs[3, ], se.log = NA, run = "Ensemble") }))
    ens <- fl[[1]]; ens$ts <- d; ens$run <- "Ensemble"
    ens$refpts[, "est"] <- apply(sapply(fl, function(f) f$refpts[, "est"]), 1, stats::median)
    p1 <- jbplot_trj(ens, type = type, ncol = 2) + ggplot2::labs(subtitle = paste("Ensemble of", nr, "runs (equal weights)"))
  }
  if(!kobe || !requireNamespace("patchwork", quietly = TRUE)) return(p1)
  p2 <- .kobe_panel(fl, ndraws, joint)
  p1 + p2 + patchwork::plot_layout(widths = c(2, 1.2))
}

# Kobe panel: each run's terminal-year uncertainty cloud and median (or the pooled cloud with joint = TRUE)
.kobe_panel <- function(fl, ndraws = 1000, joint = FALSE){
  cl <- do.call(rbind, lapply(names(fl), function(n){ f <- fl[[n]]; z <- .kobe_cloud(f, length(f$years), ndraws); z$run <- n; z }))
  if(joint) cl$run <- "Ensemble"
  runs <- unique(cl$run); cl$run <- factor(cl$run, runs)
  med <- stats::aggregate(cbind(stock, harvest) ~ run, cl, stats::median)
  cols <- stats::setNames(jb3cols(length(runs)), runs)
  xmax <- max(2, stats::quantile(cl$stock, 0.99))*1.05; ymax <- max(2, stats::quantile(cl$harvest, 0.99))*1.05
  qd <- data.frame(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1), ymax = c(1, 1, Inf, Inf),
                   zone = c("green", "yellow", "orange", "red"))
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                       fill = .fills_kobe[qd$zone], alpha = 0.7) +
    ggplot2::geom_point(data = cl, ggplot2::aes(.data$stock, .data$harvest, colour = .data$run), alpha = 0.15, size = 0.6) +
    ggplot2::geom_point(data = med, ggplot2::aes(.data$stock, .data$harvest, fill = .data$run), shape = 21, size = 3.5) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = expression(B/B[MSY]), y = expression(F/F[MSY]), subtitle = paste("Status", max(fl[[1]]$years))) +
    theme_jabba3() + ggplot2::theme(legend.position = "bottom")
  p2
}

#' Catch risk curves
#'
#' As JABBA's \code{jbplot_Crisk()}: the probability that B < Bmsy (or B < \code{bfrac} x Bmsy) and that F > Fmsy against
#' the projected catch (or F/Fmsy) of each scenario, one line per projection year, from \code{\link{jb3_project}()}.
#' @param proj output of \code{\link{jb3_project}()} with two or more scenarios
#' @param years projection years to show (default: up to 4 spread over the horizon)
#' @param bfrac also show P(B < bfrac x Bmsy), e.g. 0.5 for a limit; NULL to omit
#' @export
jbplot_crisk <- function(proj, years = NULL, bfrac = 0.5){
  it <- proj$iters; sc <- proj$scenarios
  py <- it[[1]]$years
  years <- years %||% unique(py[unique(round(seq(1, length(py), length.out = min(4, length(py)))))])
  d <- do.call(rbind, lapply(seq_along(it), function(s) do.call(rbind, lapply(years, function(y){ j <- match(y, it[[s]]$years)
    bb <- it[[s]]$`B/Bmsy`[, j]; ff <- it[[s]]$`F/Fmsy`[, j]
    out <- data.frame(x = sc$value[s], year = y, risk = c("P(B < Bmsy)", "P(F > Fmsy)"),
                      p = c(mean(bb < 1, na.rm = TRUE), mean(ff > 1, na.rm = TRUE)))
    if(!is.null(bfrac)) out <- rbind(out, data.frame(x = sc$value[s], year = y, risk = sprintf("P(B < %s Bmsy)", bfrac),
                                                     p = mean(bb < bfrac, na.rm = TRUE)))
    out }))))
  d$risk <- factor(d$risk, unique(d$risk)); d$year <- factor(d$year)
  ny <- length(years)
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$p, colour = .data$year)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = jb3cols(max(ny, 2))[seq_len(ny)]) +
    ggplot2::facet_wrap(~risk) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = if(sc$type[1] == "Catch") "Catch" else "F/Fmsy", y = "Probability", subtitle = paste(proj$run, "- risk by projection year")) +
    theme_jabba3() + ggplot2::theme(legend.position = "bottom")
}

#' Summary of one or several fits
#'
#' As JABBA's \code{jbplot_summary()}: B, F, B/Bmsy, F/Fmsy and B/K trajectories with 95\% intervals, and the surplus
#' production curves with Bmsy and MSY marked, coloured by run.
#' @param x a \code{jabba3} fit, a list of fits, or a workflow fit
#' @param steps for a workflow fit, show all steps
#' @param kobe add the Kobe panel with the terminal-year status (FALSE: trajectories and production curves only)
#' @param ndraws points per run in the Kobe cloud
#' @export
jbplot_summary <- function(x, steps = FALSE, kobe = TRUE, ndraws = 1000){
  fl <- .fitlist(x, steps); runs <- names(fl)
  p1 <- jbplot_trj(fl, type = c("B", "F", "B/Bmsy", "F/Fmsy", "B/K"), ncol = 3)
  cv <- do.call(rbind, lapply(runs, function(n){ z <- .spcurve(fl[[n]]); data.frame(run = n, z$curve) }))
  rp <- do.call(rbind, lapply(runs, function(n){ z <- .spcurve(fl[[n]]); data.frame(run = n, Bmsy = z$rp$Bmsy, MSY = z$rp$MSY) }))
  cv$run <- factor(cv$run, runs); rp$run <- factor(rp$run, runs)
  cols <- stats::setNames(if(length(runs) == 1) "black" else jb3cols(length(runs)), runs)
  p2 <- ggplot2::ggplot(cv, ggplot2::aes(.data$b, .data$sp, colour = .data$run)) + ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_segment(data = rp, ggplot2::aes(x = .data$Bmsy, xend = .data$Bmsy, y = 0, yend = .data$MSY), linetype = 2) +
    ggplot2::geom_point(data = rp, ggplot2::aes(.data$Bmsy, .data$MSY), size = 2.5) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "Biomass", y = "Surplus production") + theme_jabba3() +
    ggplot2::theme(legend.position = if(length(runs) == 1) "none" else "right")
  if(!requireNamespace("patchwork", quietly = TRUE)) return(list(p1, p2))
  p3 <- if(kobe) .kobe_panel(fl, ndraws) + ggplot2::theme(legend.position = "none") else patchwork::plot_spacer()
  p1 / (p2 + p3 + patchwork::plot_layout(widths = c(2, 1))) + patchwork::plot_layout(heights = c(2, 1))
}

#' Selex: spawning and exploitable biomass status
#'
#' As JABBA's \code{jbplot_ebsb_status()} for a Selex fit (\code{jb3_model(ebsb = )}): SB/SBmsy (the state, with its 95\%
#' interval) and, for each EB/SB curve, EB/EBmsy = (SB/SBmsy) x rho(SB/K) / rho(SBmsy/K).
#' @param fit a \code{jabba3} Selex fit
#' @export
jbplot_ebsb_status <- function(fit){
  d <- fit$setup$dat
  if(!isTRUE(d$eb.on)) stop("jbplot_ebsb_status() needs a Selex fit: jb3_model(ebsb = ...)", call. = FALSE)
  sb <- fit$ts[fit$ts$quantity == "B/Bmsy", ]; P <- fit$ts$est[fit$ts$quantity == "B/K"]; Pm <- fit$refpts["BmsyK", "est"]
  cn <- names(fit$model$ebsb) %||% paste0("sel", seq_len(nrow(d$ebm) - 1))
  eb <- do.call(rbind, lapply(seq_along(cn), function(j) data.frame(year = sb$year, curve = paste("EB/EBmsy:", cn[j]),
    est = sb$est*.rho_curve(P, d$ebm[j, ])/.rho_curve(Pm, d$ebm[j, ]))))
  nc <- length(cn); cols <- stats::setNames(jb3cols(max(nc, 2))[seq_len(nc)], unique(eb$curve))
  ggplot2::ggplot(sb, ggplot2::aes(.data$year, .data$est)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey80") +
    ggplot2::geom_line(ggplot2::aes(colour = "SB/SBmsy"), linewidth = 0.8) +
    ggplot2::geom_line(data = eb, ggplot2::aes(colour = .data$curve), linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = c(`SB/SBmsy` = "black", cols)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 5) + ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Status relative to MSY level", subtitle = fit$run) + theme_jabba3() +
    ggplot2::theme(legend.position = "bottom")
}
