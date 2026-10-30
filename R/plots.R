# ggplot2 graphics (theme, colours and styles from JAGG / FLRef::plotAdvice)

#' JABBA3 ggplot theme
#'
#' \code{theme_bw()} with the compact legend of \code{FLRef::plotAdvice()} (as JAGG's \code{theme_jagg()}).
#' @param base_size base font size
#' @export
theme_jabba3 <- function(base_size = 11){
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.key.size = ggplot2::unit(0.3, "cm"), legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width = ggplot2::unit(0.3, "cm"), legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 7), strip.background = ggplot2::element_rect(fill = "grey92"),
                   panel.spacing.x = ggplot2::unit(1, "lines"), panel.grid.minor = ggplot2::element_blank())
}

#' Integer year axis
#' @param n approximate number of breaks
#' @param ... passed to \code{scale_x_continuous()}
#' @export
scale_x_year <- function(n = 5, ...) ggplot2::scale_x_continuous(breaks = function(x) unique(round(pretty(x, n = n))), ...)

#' Colours for runs, scenarios or peels (r4ss rule, as JAGG and ss3diags)
#' @param n number of colours
#' @param alpha transparency
#' @export
jb3cols <- function(n, alpha = 1){
  rc <- function(n){ x <- seq(0, 1, length.out = n); r <- 1/(1 + exp(20 - 35*x)); g <- pmin(pmax(0, -0.8 + 6*x - 5*x^2), 1)
    b <- stats::dnorm(x, 0.25, 0.15)/max(stats::dnorm(x, 0.25, 0.15)); grDevices::rgb(r, g, b) }
  col <- if(n > 3) rc(n + 1)[-1] else if(n < 3) rc(n) else c("blue", "red", "green3")
  if(n == 1) col <- "dodgerblue"
  if(alpha < 1) col <- grDevices::adjustcolor(col, alpha.f = alpha)
  col
}

.ref_colour <- function(x) ifelse(x %in% c("Blim", "Flim"), "red", ifelse(x %in% c("Bpa", "Fpa"), "darkorange",
                                  ifelse(x %in% c("B0", "K", "Btrigger"), "blue", "darkgreen")))

.geom_refpts <- function(refs, years, size = 2.5){
  if(is.null(refs) || !nrow(refs)) return(NULL)
  refs$colour <- .ref_colour(refs$name)
  xy <- stats::quantile(seq(min(years), max(years)), c(0.2, 0.45, 0.75, 0.6, 0.3), names = FALSE)
  refs$x <- stats::ave(seq_len(nrow(refs)), refs$quantity, FUN = function(i) xy[(seq_along(i) - 1) %% length(xy) + 1])
  list(ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = .data$value), colour = refs$colour, linetype = 2, linewidth = 0.5),
       ggplot2::geom_label(data = refs, ggplot2::aes(x = .data$x, y = .data$value, label = .data$name), colour = refs$colour,
                           size = size, fill = "white", label.size = 0, label.padding = ggplot2::unit(0.1, "lines"), inherit.aes = FALSE))
}

# a fit, a workflow fit (steps = TRUE), or a list of fits -> named list of fits
.fitlist <- function(x, steps = FALSE){
  if(inherits(x, "jabba3")){
    if(steps && !is.null(x$steps)) return(x$steps)
    return(stats::setNames(list(x), x$run))
  }
  if(is.list(x) && all(vapply(x, inherits, TRUE, "jabba3"))){
    nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- vapply(x, function(f) f$run, "")
    return(stats::setNames(x, make.unique(nm)))
  }
  stop("x must be a jabba3 fit or a list of fits", call. = FALSE)
}

#' Trajectories with 95\% intervals
#'
#' One facet per quantity, the estimate (mode, or posterior median with NUTS) with its 95\% interval, reference points as
#' labelled lines (FLRef::plotAdvice style). Several fits (a list, or the steps of a workflow fit) are coloured by run.
#' @param x a \code{jabba3} fit, a list of fits, or a workflow fit
#' @param type quantities: any of "B", "F", "B/Bmsy", "F/Fmsy", "B/K", "Catch", "Process deviation", "Process SD",
#'   "Surplus production", "EB"
#' @param steps for a workflow fit, plot all steps
#' @param refpts draw reference points
#' @param extra named vector of extra reference points (names starting with B go to biomass panels, F to F panels)
#' @param ncol facet columns
#' @param colour ribbon colour for a single run
#' @export
jbplot_trj <- function(x, type = c("B", "F", "B/Bmsy", "F/Fmsy", "Catch", "B/K"), steps = FALSE, refpts = TRUE, extra = NULL,
                       ncol = NULL, colour = "dodgerblue"){
  fl <- .fitlist(x, steps)
  d <- do.call(rbind, lapply(names(fl), function(n){ t <- fl[[n]]$ts; t$run <- n; t[t$quantity %in% type, ] }))
  d$quantity <- factor(d$quantity, levels = type); d$run <- factor(d$run, levels = names(fl))
  nr <- length(fl)
  mx <- stats::ave(d$est, d$quantity, FUN = function(x) max(abs(x), na.rm = TRUE))
  d$uci[is.finite(d$uci) & d$uci > 10*mx] <- NA                   # guard against exploding intervals
  dr <- d[is.finite(d$lci) & is.finite(d$uci), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$est))
  if(nr == 1){
    if(nrow(dr)) p <- p + ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = colour, alpha = 0.35)
    p <- p + ggplot2::geom_line(linewidth = 0.6)
  } else {
    cols <- stats::setNames(jb3cols(nr), names(fl))
    if(nrow(dr)) p <- p + ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = .data$lci, ymax = .data$uci, fill = .data$run), alpha = 0.15)
    p <- p +
      ggplot2::geom_line(ggplot2::aes(colour = .data$run), linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols)
  }
  p <- p + ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jabba3()
  if(!any(type %in% c("Process deviation"))) p <- p + ggplot2::expand_limits(y = 0)
  if(refpts){
    f <- fl[[1]]; rp <- f$refpts
    refs <- data.frame(quantity = c("B/Bmsy", "F/Fmsy", "B", "B", "F", "Catch", "B/K", "B/K", "Process deviation"),
                       name = c("Bmsy", "Fmsy", "Bmsy", "K", "Fmsy", "MSY", "Bmsy", "B0", "zero"),
                       value = c(1, 1, rp["Bmsy", "est"], rp["K", "est"], rp["Fmsy", "est"], rp["MSY", "est"], rp["BmsyK", "est"], 1, 0))
    if(length(extra)){ nm <- names(extra)
      refs <- rbind(refs, data.frame(quantity = ifelse(substr(nm, 1, 1) == "F", "F/Fmsy", "B/Bmsy"), name = nm, value = unname(extra))) }
    refs <- refs[refs$quantity %in% type & refs$name != "zero", ]
    refs$quantity <- factor(refs$quantity, levels = type)
    p <- p + .geom_refpts(refs, range(d$year))
    if("Process deviation" %in% type) p <- p + ggplot2::geom_hline(data = data.frame(quantity = factor("Process deviation", levels = type)),
                                                                   ggplot2::aes(yintercept = 0), linetype = 2)
  }
  p
}

#' Advice plot
#'
#' Stock status (B/Bmsy, F/Fmsy) and catch with MSY reference points, as \code{FLRef::plotAdvice()} and JAGG's
#' \code{jagg_advice()}; \code{rel = FALSE} for B, F and catch on absolute scales.
#' @inheritParams jbplot_trj
#' @param rel relative (default) or absolute quantities
#' @param ... passed to \code{jbplot_trj()}
#' @export
jbplot_advice <- function(x, rel = TRUE, ...) jbplot_trj(x, type = if(rel) c("B/Bmsy", "F/Fmsy", "Catch") else c("B", "F", "Catch"), ...)

#' Process error deviations and process SD
#' @inheritParams jbplot_trj
#' @export
jbplot_procdev <- function(x, steps = FALSE){
  fl <- .fitlist(x, steps)
  ty <- if(any(vapply(fl, function(f) f$model$sdb.form == "dyn", TRUE))) c("Process deviation", "Process SD") else "Process deviation"
  jbplot_trj(x, type = ty, steps = steps, ncol = 1) + ggplot2::labs(y = "log scale")
}

#' Kobe plot
#'
#' B/Bmsy against F/Fmsy: trajectory and the joint uncertainty in the last year, from the draws (NUTS or MVN) or, without
#' draws, from the bivariate normal of the terminal log B/Bmsy and log F/Fmsy (delta method), with the percentage in each
#' quadrant.
#' @param fit a \code{jabba3} fit
#' @param npoints number of points in the cloud
#' @export
jbplot_kobe <- function(fit, npoints = 2000){
  N <- length(fit$years); y <- max(fit$years)
  bb <- fit$ts$est[fit$ts$quantity == "B/Bmsy"]; ff <- fit$ts$est[fit$ts$quantity == "F/Fmsy"]
  trj <- data.frame(year = fit$years, stock = bb, harvest = ff)
  if(!is.null(fit$draws)) pt <- data.frame(stock = fit$draws$`B/Bmsy`[, N], harvest = fit$draws$`F/Fmsy`[, N]) else {
    s <- fit$sdr; nm <- names(s$value)
    ib <- which(nm == "logBBmsy")[N]; ih <- which(nm == "logFFmsy")[N]
    S <- s$cov[c(ib, ih), c(ib, ih)]; set.seed(1)
    z <- matrix(stats::rnorm(2*npoints), ncol = 2) %*% chol(S)
    pt <- data.frame(stock = exp(s$value[ib] + z[, 1]), harvest = exp(s$value[ih] + z[, 2]))
  }
  zone <- ifelse(pt$stock >= 1 & pt$harvest <= 1, "green", ifelse(pt$stock < 1 & pt$harvest <= 1, "yellow",
                 ifelse(pt$stock >= 1 & pt$harvest > 1, "orange", "red")))
  pct <- round(100*prop.table(table(factor(zone, c("green", "yellow", "orange", "red")))), 1)
  if(nrow(pt) > npoints) pt <- pt[sample(nrow(pt), npoints), ]
  xmax <- max(2, trj$stock, stats::quantile(pt$stock, 0.99))*1.05; ymax <- max(2, trj$harvest, stats::quantile(pt$harvest, 0.99), na.rm = TRUE)*1.05
  qd <- data.frame(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1), ymax = c(1, 1, Inf, Inf),
                   zone = c("green", "yellow", "orange", "red"))
  qd$lab <- sprintf("%.1f%%", pct[qd$zone])
  fills <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")
  ends <- trj[c(1, nrow(trj)), ]
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$zone), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = fills, breaks = qd$zone, labels = qd$lab, name = paste("Uncertainty", y)) +
    ggplot2::geom_point(data = pt, ggplot2::aes(.data$stock, .data$harvest), colour = "grey30", alpha = 0.15, size = 0.6) +
    ggplot2::geom_path(data = trj, ggplot2::aes(.data$stock, .data$harvest), linewidth = 0.6) +
    ggplot2::geom_point(data = trj, ggplot2::aes(.data$stock, .data$harvest), shape = 21, fill = "white", size = 1.4) +
    ggplot2::geom_point(data = ends, ggplot2::aes(.data$stock, .data$harvest), shape = c(22, 24), fill = "white", size = 3) +
    ggplot2::geom_text(data = ends, ggplot2::aes(.data$stock, .data$harvest, label = .data$year), vjust = -1.1, size = 3) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = expression(B/B[MSY]), y = expression(F/F[MSY]), subtitle = fit$run) +
    theme_jabba3() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8), legend.position = "right")
}

#' Fits to the abundance indices
#'
#' Observations, fitted index with its 95\% interval and the 95\% predictive interval (grey), one facet per index.
#' @inheritParams jbplot_trj
#' @param log log scale
#' @param ppd draw the predictive interval
#' @export
jbplot_fits <- function(x, log = FALSE, ppd = TRUE, steps = FALSE, colour = "dodgerblue"){
  fl <- .fitlist(x, steps)
  d <- do.call(rbind, lapply(names(fl), function(n){ i <- fl[[n]]$index; i$run <- n; i }))
  d$index <- factor(d$index, levels = unique(d$index))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year))
  if(length(fl) == 1){
    if(ppd) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ppd_lci, ymax = .data$ppd_uci), fill = "grey85")
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = colour, alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = .data$hat), colour = colour, linewidth = 0.7)
  } else {
    cols <- stats::setNames(jb3cols(length(fl)), names(fl))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$hat, colour = .data$run), linewidth = 0.7) + ggplot2::scale_colour_manual(values = cols)
  }
  o <- d[d$run == names(fl)[1], ]
  p <- p + ggplot2::geom_point(data = o, ggplot2::aes(y = .data$obs), shape = 21, fill = "white", size = 1.8) +
    ggplot2::facet_wrap(~index, scales = "free_y") + scale_x_year() + ggplot2::labs(x = "Year", y = "Index") + theme_jabba3()
  if(log) p + ggplot2::scale_y_log10() else p + ggplot2::expand_limits(y = 0)
}

.resid <- function(fit, type){
  d <- fit$index; d$r <- if(type == "osa") d$osa else d$residual
  d[is.finite(d$r), ]
}

#' Index residuals by year
#'
#' As JABBA's \code{jbplot_residuals()}: residuals of all indices coloured by index, a boxplot per year, a loess smoother
#' and the RMSE; OSA residuals (N(0,1) under the model) by default, or ordinary log residuals.
#' @param fit a \code{jabba3} fit
#' @param type "osa" (default) or "ordinary"
#' @export
jbplot_residuals <- function(fit, type = c("osa", "ordinary")){
  type <- match.arg(type); d <- .resid(fit, type)
  ni <- length(unique(d$index))
  sub <- if(type == "osa") sprintf("OSA residuals: mean %.2f, sd %.2f", mean(d$r), stats::sd(d$r)) else sprintf("RMSE = %.1f%%", 100*sqrt(mean(d$r^2)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$r)) + ggplot2::geom_hline(yintercept = 0, linetype = 2)
  if(type == "osa") p <- p + ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 3, colour = "grey50")
  p + ggplot2::geom_boxplot(ggplot2::aes(group = .data$year), fill = "grey90", colour = "grey60", outlier.shape = NA, width = 0.7) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$index), alpha = 0.6) + ggplot2::geom_point(ggplot2::aes(colour = .data$index), size = 1.6) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) + scale_x_year() +
    ggplot2::labs(x = "Year", y = if(type == "osa") "OSA residual" else "log residual", subtitle = sub) + theme_jabba3()
}

#' Runs test
#'
#' Residuals as segments from zero, three-sigma limits (green: runs test passed, p >= 0.05; red: failed), outliers in red,
#' after \code{ggplotFL::plotRunstest()} and JABBA's \code{jbplot_runstest()}. OSA residuals by default.
#' @inheritParams jbplot_residuals
#' @export
jbplot_runs <- function(fit, type = c("osa", "ordinary")){
  type <- match.arg(type); d <- .resid(fit, type)
  rt <- do.call(rbind, lapply(split(d, factor(d$index, unique(d$index))), function(di){
    r <- jb3_runs(di$r[order(di$year)])
    data.frame(index = di$index[1], p.runs = r$p.runs, lcl = r$sig3lim[1], ucl = r$sig3lim[2],
               test = ifelse(r$p.runs < 0.05, "Failed", "Passed"), xmin = min(di$year) - 0.5, xmax = max(di$year) + 0.5) }))
  d <- merge(d, rt[, c("index", "lcl", "ucl")], by = "index"); d$outlier <- d$r < d$lcl | d$r > d$ucl
  lab <- stats::setNames(sprintf("%s\np = %.3f", rt$index, rt$p.runs), rt$index)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$r)) +
    ggplot2::geom_rect(data = rt, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$lcl, ymax = .data$ucl, fill = .data$test),
                       inherit.aes = FALSE, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(Passed = "#cbe368", Failed = "#ef8575"), drop = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_segment(ggplot2::aes(xend = .data$year, yend = 0)) +
    ggplot2::geom_point(size = 1.8) + ggplot2::geom_point(ggplot2::aes(colour = .data$outlier), size = 1.1) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "#d50102"), guide = "none") +
    ggplot2::facet_wrap(~index, scales = "free_x", labeller = ggplot2::as_labeller(lab)) + scale_x_year(n = 4) +
    ggplot2::labs(x = "Year", y = if(type == "osa") "OSA residual" else "log residual") + theme_jabba3() +
    ggplot2::theme(legend.position = "bottom")
}

#' SPiCT-style residual diagnostics
#'
#' Per index, as \code{spict::plotspict.diagnostic()}: OSA residuals over time (bias t-test), autocorrelation function
#' (Ljung-Box test) and normal QQ plot (Shapiro-Wilk test), p-values in the panel titles, in red when < 0.05.
#' @param fit a \code{jabba3} fit (with OSA residuals)
#' @return a patchwork of ggplots (or a list without the patchwork package)
#' @export
jbplot_spdiags <- function(fit){
  dg <- jb3_diags(fit, "osa"); ix <- fit$index
  rows <- lapply(seq_len(nrow(dg)), function(k){
    nm <- dg$index[k]; d <- ix[ix$index == nm & is.finite(ix$osa), ]; d <- d[order(d$year), ]
    col <- function(p) if(is.finite(p) && p < 0.05) "red" else "darkgreen"
    tt <- function(lab, p) list(ggplot2::labs(title = sprintf("%s: %s p = %.2f", nm, lab, p)),
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9, colour = col(p))))
    a <- stats::acf(d$osa, plot = FALSE, lag.max = min(10, nrow(d) - 2))
    ad <- data.frame(lag = a$lag[-1], acf = a$acf[-1]); ci <- 1.96/sqrt(nrow(d))
    p1 <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$osa)) + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_line(colour = "grey60") + ggplot2::geom_point(colour = "dodgerblue", size = 1.6) +
      ggplot2::labs(x = "Year", y = "OSA residual") + scale_x_year(n = 4) + theme_jabba3(9) + tt("bias", dg$bias.p[k])
    p2 <- ggplot2::ggplot(ad, ggplot2::aes(.data$lag, .data$acf)) + ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = 2, colour = "blue") +
      ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0)) +
      ggplot2::labs(x = "Lag", y = "ACF") + theme_jabba3(9) + tt("Ljung-Box", dg$LjungBox.p[k])
    p3 <- ggplot2::ggplot(d, ggplot2::aes(sample = .data$osa)) + ggplot2::stat_qq(colour = "dodgerblue", size = 1.4) +
      ggplot2::stat_qq_line(linetype = 2) + ggplot2::labs(x = "Theoretical", y = "Sample") + theme_jabba3(9) + tt("Shapiro", dg$Shapiro.p[k])
    list(p1, p2, p3)
  })
  pl <- unlist(rows, recursive = FALSE)
  if(requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pl, ncol = 3) else pl
}

#' Prior and posterior distributions
#'
#' Prior densities (lognormal; flat priors not drawn) against the posterior (draws, or the lognormal approximation at the
#' mode), with the posterior-to-prior variance ratio (PPVR, log scale) after JABBA's \code{jbplot_ppdist()}. K is shown on a
#' log axis.
#' @param fit a \code{jabba3} fit
#' @param pars parameters to show
#' @export
jbplot_priorpost <- function(fit, pars = c("r", "K", "bkinit", "sdb", "m", "kdyn")){
  pr <- fit$priors; pp <- fit$pars; dat <- fit$setup$dat
  lst <- list(r = pr$r, K = c(dat$K.pr[1]*fit$mucatch, pr$K[2:3]), bkinit = pr$bkinit, sdb = pr$sdb, m = pr$m)
  pl <- lapply(pars, function(nm){
    if(!nm %in% rownames(pp)) return(NULL)
    if((nm == "kdyn" && !dat$est.k) || (nm == "m" && !dat$est.m) || (nm == "sdb" && !dat$est.sdb)) return(NULL)
    if(nm == "kdyn"){ mu <- pr$kdyn[[1]]; s <- pr$kdyn[[2]]; on <- TRUE; logsc <- FALSE } else {
      x <- lst[[nm]]; mu <- log(x[[1]]); s <- x[[2]]; on <- isTRUE(x[[3]] == 1) && is.finite(mu); logsc <- TRUE }
    if(nm == "sdb" && dat$sdb.ig) on <- FALSE
    post <- if(!is.null(fit$draws) && nm %in% names(fit$draws$pars)) fit$draws$pars[[nm]] else NULL
    tr <- function(v) if(logsc) log(v) else v
    pm <- tr(pp[nm, "est"]); ps <- if(logsc) pp[nm, "se.log"] else (pp[nm, "uci"] - pp[nm, "lci"])/3.92
    if(!is.null(post)) { pm <- stats::median(tr(post)); ps <- stats::sd(tr(post)) }
    if(!is.finite(ps)) return(NULL)
    lo <- min(if(on) mu - 2.6*s, pm - 4*ps); hi <- max(if(on) mu + 2.6*s, pm + 4*ps)
    z <- seq(lo, hi, length.out = 400)                       # transformed scale
    if(!is.null(post)){ dd <- stats::density(tr(post)); fz <- stats::approx(dd$x, dd$y, z, yleft = 0, yright = 0)$y } else fz <- stats::dnorm(z, pm, ps)
    xs <- if(logsc) exp(z) else z; jac <- if(logsc) 1/xs else 1
    d <- rbind(if(on) data.frame(x = xs, y = stats::dnorm(z, mu, s)*jac, type = "Prior"), data.frame(x = xs, y = fz*jac, type = "Posterior"))
    g <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y, fill = .data$type, colour = .data$type)) +
      ggplot2::geom_area(position = "identity", alpha = 0.35) +
      ggplot2::scale_fill_manual(values = c(Prior = "grey60", Posterior = "dodgerblue"), drop = FALSE) +
      ggplot2::scale_colour_manual(values = c(Prior = "grey40", Posterior = "dodgerblue4"), drop = FALSE) +
      ggplot2::labs(x = NULL, y = "Density", title = sprintf("%s  (PPVR = %s)", nm, if(on) sprintf("%.2f", ps^2/s^2) else "flat prior")) +
      theme_jabba3(9) + ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 9))
    if(nm == "K") g <- g + ggplot2::scale_x_log10()
    g
  })
  pl <- pl[!vapply(pl, is.null, TRUE)]
  if(requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pl) + patchwork::plot_annotation(
    subtitle = paste0("grey: prior, blue: posterior (", if(!is.null(fit$draws)) fit$draws$type else "mode, lognormal approximation", ")")) else pl
}

#' Surplus production phase and dynamics plots
#'
#' \code{jbplot_spphase()}: the surplus production curve with the catch trajectory (yield against biomass) on the
#' Kobe colours mapped to production space (as \code{JABBA::jbplot_spphase()}): green B > Bmsy and F < Fmsy (below the
#' line C = Fmsy B), orange B > Bmsy and F > Fmsy, yellow B < Bmsy and F < Fmsy, red B < Bmsy and F > Fmsy. The hatched
#' band below Bmsy marks catches between Fmsy B and the production curve, where F < Fmsy and the stock rebuilds.
#' Blue band: MSY 95\% interval.
#'
#' \code{jbplot_spdyn()}: the realised annual surplus production SP_t = B_t+1 - B_t + C_t against B_t, with the
#' expected production curve (as \code{JABBA::jbplot_spdyn()}); departures from the curve are the process deviations
#' (and the lag) in production units.
#' @param fit a \code{jabba3} fit
#' @param rel relative axes: B/K and catch (or SP) / MSY
#' @return a ggplot
#' @export
jbplot_spphase <- function(fit, rel = FALSE){
  z <- .spcurve(fit, rel); rp <- z$rp
  xmax <- max(z$curve$b, z$trj$B, na.rm = TRUE)*1.02
  ymax <- max(z$trj$C*1.05, rp$MSYu*1.05, rp$MSY*1.1, na.rm = TRUE); big <- 10*max(xmax, ymax)
  H <- rp$MSY/rp$Bmsy                                        # Fmsy (in plot units): C = H B is F = Fmsy
  zones <- rbind(
    data.frame(zone = "red", x = c(0, rp$Bmsy, rp$Bmsy, 0), y = c(0, rp$MSY, big, big)),
    data.frame(zone = "yellow", x = c(0, rp$Bmsy, rp$Bmsy), y = c(0, 0, rp$MSY)),
    data.frame(zone = "green", x = c(rp$Bmsy, big, big, rp$Bmsy), y = c(0, 0, H*big, rp$MSY)),
    data.frame(zone = "orange", x = c(rp$Bmsy, big, big, rp$Bmsy), y = c(rp$MSY, H*big, big, big)))
  fills <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")
  lo <- z$curve[z$curve$b <= rp$Bmsy, ]; lo$h <- H*lo$b
  hatch <- lo[seq(1, nrow(lo), by = max(1, round(nrow(lo)/30))), ]
  hatch$col <- rep(c("#de2d26", "#fee391"), length.out = nrow(hatch))
  sel <- .selyears(z$trj)
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = zones, ggplot2::aes(.data$x, .data$y, group = .data$zone, fill = .data$zone), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::geom_segment(data = hatch, ggplot2::aes(x = .data$b, xend = .data$b, y = .data$h, yend = .data$sp),
                          colour = hatch$col, linetype = 3, linewidth = 0.5) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = rp$MSYl, ymax = rp$MSYu, fill = "blue", alpha = 0.15) +
    ggplot2::geom_line(data = z$curve, ggplot2::aes(.data$b, .data$sp), colour = "blue", linewidth = 1) +
    ggplot2::geom_hline(yintercept = rp$MSY, colour = "blue", linetype = 5) +
    ggplot2::annotate("segment", x = rp$Bmsy, xend = rp$Bmsy, y = 0, yend = rp$MSY, colour = "blue", linetype = 2) +
    ggplot2::geom_path(data = z$trj, ggplot2::aes(.data$B, .data$C), na.rm = TRUE) +
    ggplot2::geom_point(data = z$trj, ggplot2::aes(.data$B, .data$C), size = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(data = sel, ggplot2::aes(.data$B, .data$C, shape = .data$lab), fill = "white", size = 3.2) +
    ggplot2::scale_shape_manual(values = c(22, 21, 24), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = if(rel) "B/K" else "Biomass", y = if(rel) "Catch/MSY" else "Yield", subtitle = fit$run) +
    theme_jabba3() + ggplot2::theme(legend.position = "bottom")
}

#' @rdname jbplot_spphase
#' @export
jbplot_spdyn <- function(fit, rel = FALSE){
  z <- .spcurve(fit, rel); rp <- z$rp; sel <- .selyears(z$trj, "SP")
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = rp$MSY, linetype = 3, colour = "grey40") +
    ggplot2::geom_vline(xintercept = rp$Bmsy, linetype = 3, colour = "grey40") +
    ggplot2::geom_line(data = z$curve, ggplot2::aes(.data$b, .data$sp, linetype = "Expected"), colour = "blue", linewidth = 1) +
    ggplot2::geom_path(data = z$trj, ggplot2::aes(.data$B, .data$SP), na.rm = TRUE) +
    ggplot2::geom_point(data = z$trj, ggplot2::aes(.data$B, .data$SP), size = 1.2, na.rm = TRUE) +
    ggplot2::geom_point(data = sel, ggplot2::aes(.data$B, .data$SP, shape = .data$lab), fill = "white", size = 3.2) +
    ggplot2::scale_shape_manual(values = c(22, 21, 24), name = NULL) +
    ggplot2::scale_linetype_manual(values = 1, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.03)), limits = c(0, NA)) +
    ggplot2::labs(x = if(rel) "B/K" else "Biomass", y = if(rel) "SP/MSY" else "Surplus production",
                  subtitle = paste(fit$run, "- realised SP(t) = B(t+1) - B(t) + C(t)")) +
    theme_jabba3() + ggplot2::theme(legend.position = "bottom")
}

#' @rdname jbplot_spphase
#' @param type "phase" (\code{jbplot_spphase}) or "dyn" (\code{jbplot_spdyn})
#' @export
jbplot_production <- function(fit, type = c("phase", "dyn"), rel = FALSE){
  type <- match.arg(type)
  if(type == "phase") jbplot_spphase(fit, rel) else jbplot_spdyn(fit, rel)
}

.spcurve <- function(fit, rel = FALSE){
  rp <- fit$refpts; K <- rp["K", "est"]; MSY <- rp["MSY", "est"]; Bmsy <- rp["Bmsy", "est"]
  r <- fit$pars["r", "est"]; m <- fit$pars["m", "est"]
  g <- function(q) fit$ts$est[fit$ts$quantity == q]
  B <- g("B"); b <- seq(0, max(K, B, na.rm = TRUE), length.out = 400)
  sp <- r/(m - 1)*b*(1 - (b/K)^(m - 1))
  trj <- data.frame(year = fit$years, B = B, C = g("Catch"), SP = g("Surplus production"))
  sx <- if(rel) K else 1; sy <- if(rel) MSY else 1
  trj$B <- trj$B/sx; trj$C <- trj$C/sy; trj$SP <- trj$SP/sy
  list(curve = data.frame(b = b/sx, sp = sp/sy), trj = trj,
       rp = list(Bmsy = Bmsy/sx, MSY = MSY/sy, MSYl = rp["MSY", "lci"]/sy, MSYu = rp["MSY", "uci"]/sy))
}

.selyears <- function(trj, y = "C"){
  ok <- which(is.finite(trj[[y]])); N <- length(ok)
  s <- trj[ok[c(1, round(stats::quantile(seq_len(N), 0.7)), N)], ]
  s$lab <- factor(s$year, levels = s$year); s
}

#' Data overview
#' @param data a \code{jb3data} object or a fit
#' @export
jbplot_data <- function(data){
  if(inherits(data, "jabba3")) data <- data$data
  c1 <- data.frame(year = data$years, value = data$catch, panel = "Catch", series = "Catch")
  I <- data$I; I <- sweep(I, 2, colMeans(I, na.rm = TRUE), "/")
  i1 <- do.call(rbind, lapply(colnames(I), function(n) data.frame(year = data$years, value = I[, n], panel = "Indices (scaled to mean 1)",
                                                                  series = n, se = data$SE[, n])))
  i1 <- i1[is.finite(i1$value), ]
  ni <- length(unique(i1$series))
  ggplot2::ggplot(mapping = ggplot2::aes(.data$year, .data$value)) +
    ggplot2::geom_col(data = c1, fill = "grey60") +
    ggplot2::geom_line(data = i1, ggplot2::aes(colour = .data$series)) +
    ggplot2::geom_point(data = i1, ggplot2::aes(colour = .data$series), size = 1.4) +
    ggplot2::scale_colour_manual(values = jb3cols(max(ni, 2))[seq_len(ni)]) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) + scale_x_year() + ggplot2::labs(x = "Year", y = NULL, subtitle = data$stock) +
    ggplot2::expand_limits(y = 0) + theme_jabba3()
}

#' Retrospective plot with Mohn's rho
#' @param retro output of \code{\link{jb3_retro}()}
#' @param type quantities
#' @export
jbplot_retro <- function(retro, type = c("B/Bmsy", "F/Fmsy", "B", "F")){
  ts <- retro$ts[retro$ts$quantity %in% type, ]; rho <- retro$rho
  lab <- stats::setNames(sprintf("%s\nMohn's rho = %.3f", rho$quantity, rho$rho), rho$quantity)
  ts$quantity <- factor(ts$quantity, levels = type)
  peels <- sort(unique(ts$peel)); ty <- max(ts$year[ts$peel == 0])
  plab <- ifelse(peels == 0, "Reference", as.character(ty - peels)); ts$run <- factor(plab[match(ts$peel, peels)], levels = plab)
  cols <- stats::setNames(c("black", jb3cols(max(length(peels) - 1, 1))[seq_len(length(peels) - 1)]), plab)
  ends <- ts[stats::ave(ts$year, ts$run, ts$quantity, FUN = max) == ts$year & ts$peel > 0, ]
  r0 <- ts[ts$peel == 0, ]
  ggplot2::ggplot(ts, ggplot2::aes(.data$year, .data$est, colour = .data$run)) +
    ggplot2::geom_ribbon(data = r0, ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey80", colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(data = ends, size = 2) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", labeller = ggplot2::as_labeller(lab)) + ggplot2::expand_limits(y = 0) +
    scale_x_year() + ggplot2::labs(x = "Year", y = NULL) + theme_jabba3()
}

#' Hindcast cross-validation plot with MASE
#' @param hc output of \code{\link{jb3_hindcast}()}
#' @param nyears years before the hindcast horizon to show
#' @export
jbplot_hcxval <- function(hc, nyears = 5){
  d <- hc$pred; ty <- hc$terminal; np <- hc$n; ms <- hc$mase
  keep <- ms$index[ms$n.eval > 0]; d <- d[d$index %in% keep & d$year >= ty - np - nyears, ]; ms <- ms[ms$index %in% keep, ]
  peels <- 0:np; plab <- ifelse(peels == 0, "Reference", as.character(ty - peels))
  cols <- stats::setNames(c("black", jb3cols(max(np, 1))[seq_len(np)]), plab)
  ln <- d[is.finite(d$hat) & d$year <= ty - d$peel + ifelse(d$peel > 0, 1, 0), ]; ln$run <- factor(plab[match(ln$peel, peels)], levels = plab)
  pr <- ln[ln$peel > 0 & ln$year == ty - ln$peel + 1, ]
  ob <- d[d$peel == 0 & is.finite(d$obs), ]
  lab <- stats::setNames(sprintf("%s\nMASE = %.2f", ms$index, ms$MASE), ms$index)
  ggplot2::ggplot(ln, ggplot2::aes(.data$year, .data$hat, colour = .data$run)) +
    ggplot2::geom_line(data = ob, ggplot2::aes(y = .data$obs), linetype = 2, colour = "grey60") +
    ggplot2::geom_point(data = ob, ggplot2::aes(y = .data$obs), shape = 21, colour = "black", fill = "white", size = 2.8) +
    ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(data = pr, size = 2) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::facet_wrap(~index, scales = "free_y", labeller = ggplot2::as_labeller(lab)) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Index") + theme_jabba3() + ggplot2::theme(legend.position = "bottom")
}

#' Projection plot
#' @param proj output of \code{\link{jb3_project}()}
#' @param type quantities
#' @param nyears assessment years to show before the projection
#' @export
jbplot_projection <- function(proj, type = c("B/Bmsy", "F/Fmsy", "Catch"), nyears = 15){
  ft <- proj$fit.ts; y0 <- max(ft$year)
  ft <- ft[ft$quantity %in% type & ft$year > y0 - nyears, ]
  pr <- proj$ts[proj$ts$quantity %in% type, ]
  sc <- unique(pr$scenario); cols <- stats::setNames(jb3cols(length(sc)), sc)
  last <- do.call(rbind, lapply(sc, function(s){ x <- ft[ft$year == y0, ]; x$scenario <- s; x[, c("scenario", "year", "quantity", "est", "lci", "uci")] }))
  pr <- rbind(last, pr[, names(last)])
  ft$quantity <- factor(ft$quantity, levels = type); pr$quantity <- factor(pr$quantity, levels = type)
  ggplot2::ggplot(mapping = ggplot2::aes(.data$year, .data$est)) +
    ggplot2::geom_ribbon(data = ft, ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey80") +
    ggplot2::geom_line(data = ft, linewidth = 0.7) +
    ggplot2::geom_ribbon(data = pr, ggplot2::aes(ymin = .data$lci, ymax = .data$uci, fill = .data$scenario), alpha = 0.1) +
    ggplot2::geom_line(data = pr, ggplot2::aes(colour = .data$scenario), linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::geom_hline(data = data.frame(quantity = factor(intersect(type, c("B/Bmsy", "F/Fmsy")), levels = type)), ggplot2::aes(yintercept = 1), linetype = 2) +
    ggplot2::geom_vline(xintercept = y0, linetype = 3) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = 1) + ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jabba3()
}

#' Catch: observed and estimated
#' @param fit a \code{jabba3} fit
#' @export
jbplot_catch <- function(fit){
  d <- fit$ts[fit$ts$quantity == "Catch", ]; d$obs <- fit$data$catch[seq_len(nrow(d))]
  ggplot2::ggplot(d, ggplot2::aes(.data$year)) + ggplot2::geom_col(ggplot2::aes(y = .data$obs), fill = "grey70") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "dodgerblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$est), colour = "dodgerblue4", linewidth = 0.7) +
    scale_x_year() + ggplot2::labs(x = "Year", y = "Catch", subtitle = if(fit$model$catch.error) "observed (bars) and estimated" else "catch assumed known") +
    theme_jabba3()
}

#' Workflow comparison
#'
#' Compares the steps of a workflow fit (or any list of fits): B/Bmsy and F/Fmsy trajectories, and a table of the step
#' settings and terminal status in the subtitle.
#' @param x workflow fit or list of fits
#' @export
jbplot_steps <- function(x){
  fl <- .fitlist(x, steps = TRUE)
  st <- vapply(names(fl), function(n){ f <- fl[[n]]; b <- utils::tail(f$ts$est[f$ts$quantity == "B/Bmsy"], 1)
    sprintf("%s: B/Bmsy %.2f, sdb %.3f%s", n, b, f$pars["sdb", "est"], if(f$model$lag > 0) paste0(", lag ", f$model$lag) else "") }, "")
  jbplot_trj(fl, type = c("B/Bmsy", "F/Fmsy")) + ggplot2::labs(subtitle = paste(st, collapse = "\n"))
}

#' Save the standard plot set
#'
#' Writes the plots of a fit to PNG files (data, advice, trajectories, fits, residuals, runs test, SPiCT diagnostics,
#' Kobe, prior-posterior, production, process error; retrospective, hindcast and projections when given).
#' @param fit a \code{jabba3} fit
#' @param dir output folder
#' @param retro,hc,proj optional outputs of \code{jb3_retro()}, \code{jb3_hindcast()}, \code{jb3_project()}
#' @param width,height,dpi passed to \code{ggsave()}
#' @return invisibly, the file names
#' @export
jbplot_all <- function(fit, dir = ".", retro = NULL, hc = NULL, proj = NULL, width = 8, height = 5.5, dpi = 150){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  pl <- list(data = jbplot_data(fit), advice = jbplot_advice(fit), trajectories = jbplot_trj(fit), fits = jbplot_fits(fit),
             residuals_osa = jbplot_residuals(fit, "osa"), residuals = jbplot_residuals(fit, "ordinary"), runs = jbplot_runs(fit),
             spict_diags = jbplot_spdiags(fit), kobe = jbplot_kobe(fit), priorpost = jbplot_priorpost(fit),
             spphase = jbplot_spphase(fit), spdyn = jbplot_spdyn(fit), procdev = jbplot_procdev(fit))
  if(!is.null(fit$steps)) pl$steps <- jbplot_steps(fit)
  if(!is.null(fit$nuts$stanfit)){ pl$nuts_trace <- jbplot_trace(fit); pl$nuts_rank <- jbplot_rank(fit); pl$nuts_pairs <- jbplot_pairs(fit)
    pl$nuts_rhat <- jbplot_rhat(fit); pl$nuts_energy <- jbplot_energy(fit) }
  if(!is.null(retro)) pl$retro <- jbplot_retro(retro)
  if(!is.null(hc)) pl$hindcast <- jbplot_hcxval(hc)
  if(!is.null(proj)) pl$projection <- jbplot_projection(proj)
  fn <- file.path(dir, paste0(fit$run, "_", names(pl), ".png"))
  for(i in seq_along(pl)){
    h <- if(names(pl)[i] == "spict_diags") 2.6*length(unique(fit$index$index)) else if(names(pl)[i] %in% c("nuts_pairs", "nuts_rhat")) width else height
    ggplot2::ggsave(fn[i], pl[[i]], width = width, height = h, dpi = dpi)
  }
  invisible(fn)
}
