# Multi-panel summary (spict plot() style), panels selectable as in FLRef::ploteq()

#' JABBA3 summary plot with selectable panels
#'
#' A grid of panels (3 x 3 by default, as \code{plot()} for a SPiCT fit) built from the individual JABBA3 plots, for a
#' single fit or a named list of fits (runs coloured). Panels are chosen by name, in order, as in
#' \code{FLRef::ploteq(panels = )}:
#' \describe{
#'   \item{"B", "F", "Catch", "B/Bmsy", "F/Fmsy", "B/K", "Process deviation"}{trajectories with intervals
#'     (Bthr/Blim on B/Bmsy)}
#'   \item{"kobe"}{Kobe phase plot (single fit: with the terminal-year uncertainty; several fits: one trajectory per run)}
#'   \item{"production"}{surplus production curves with Bmsy/MSY and the catch path (single fit)}
#'   \item{"fits"}{index fits}
#'   \item{"residuals"}{OSA residuals (first run)}
#'   \item{"procdev"}{process deviations}
#'   \item{"priorpost"}{prior and posterior of r}
#' }
#' \code{plot(fit)} calls \code{jabbaplot3(fit)}.
#' @param x a \code{jabba3} fit, a named list of fits or a workflow fit
#' @param type "full" (3 x 3, as SPiCT \code{plot()}) or "simple" (2 x 2: B/Bmsy, F/Fmsy, Catch, Kobe, as SPiCT
#'   \code{plot2()}); ignored when \code{panels} is given
#' @param panels panel names or their positions in the full set
#' @param ncol number of columns
#' @param nested nested 90/60\% bands on the trajectories (single fit)
#' @param bthr,blim as in \code{\link{jbplot_trj}()}
#' @param obs single fit: add the observed indices scaled to biomass (I/q) and to B/Bmsy (I/(q Bmsy)) on the "B" and
#'   "B/Bmsy" panels, coloured by index (as SPiCT)
#' @param base_size font size
#' @return a patchwork (or a list of ggplots if patchwork is not installed)
#' @examples
#' \dontrun{
#' jabbaplot3(fit)
#' plot(fit)
#' jabbaplot3(list(Fox = fit, Lag.AR1 = rec), panels = c("B/Bmsy", "F/Fmsy", "kobe", "production"), ncol = 2)
#' }
#' @export
jabbaplot3 <- function(x, type = c("full", "simple"), panels = NULL, ncol = NULL, nested = FALSE, obs = TRUE, bthr = 0.5,
                       blim = 0.3, base_size = 9){
  type <- match.arg(type)
  all <- c("B", "F", "Catch", "B/Bmsy", "F/Fmsy", "kobe", "production", "fits", "residuals")
  if(is.null(panels)) panels <- if(type == "full") all else c("B/Bmsy", "F/Fmsy", "Catch", "kobe")
  if(is.null(ncol)) ncol <- if(length(panels) <= 4) 2 else 3
  if(is.numeric(panels)) panels <- all[panels]
  ok <- c(all, "B/K", "Process deviation", "procdev", "priorpost")
  bad <- setdiff(panels, ok); if(length(bad)) stop("unknown panels: ", paste(bad, collapse = ", "), "; use any of ", paste(ok, collapse = ", "), call. = FALSE)
  fl <- .fitlist(x); multi <- length(fl) > 1; f1 <- fl[[1]]
  th <- theme_jabba3(base_size) + ggplot2::theme(legend.position = if(multi) "bottom" else "none")
  ttl <- function(p, t){                                   # strip title as on the facetted panels
    if(is.data.frame(p$data)) p$data$.panel <- t
    for(i in seq_along(p$layers)) if(is.data.frame(p$layers[[i]]$data)) p$layers[[i]]$data$.panel <- t
    p + ggplot2::labs(title = NULL, subtitle = NULL) + ggplot2::facet_wrap(~.panel) + th
  }
  one <- function(pn){
    if(pn %in% c("B", "F", "Catch", "B/Bmsy", "F/Fmsy", "B/K", "Process deviation")){
      p <- jbplot_trj(if(multi) fl else f1, type = pn, nested = nested && !multi, bthr = bthr, blim = blim, label.size = 2) +
        ggplot2::labs(x = NULL) + th
      if(obs && !multi && pn %in% c("B", "B/Bmsy")){
        io <- .index_scaled(f1, rel = pn == "B/Bmsy"); io$quantity <- factor(pn, levels = pn)
        p <- p + ggplot2::geom_point(data = io, ggplot2::aes(.data$year, .data$value, fill = .data$index), shape = 21,
                                     size = 1.3, colour = "grey20", inherit.aes = FALSE) +
          ggplot2::scale_fill_manual(values = stats::setNames(jb3cols(length(unique(io$index))), unique(io$index)), name = NULL) +
          ggplot2::theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
                         legend.background = ggplot2::element_rect(fill = "white", colour = NA))
      }
      return(p)
    }
    switch(pn,
      kobe = ttl(if(multi) .kobe_runs(fl) + ggplot2::guides(colour = "none") else jbplot_kobe(f1) + ggplot2::guides(fill = "none"), "Kobe"),
      production = ttl(.prod_runs(fl) + ggplot2::guides(colour = if(multi) "none" else "legend"), "Surplus production"),
      fits = jbplot_fits(if(multi) fl else f1) + ggplot2::labs(x = NULL) + th,
      residuals = jbplot_residuals(f1) + ggplot2::labs(x = NULL, subtitle = NULL) + th,
      procdev = jbplot_procdev(if(multi) fl else f1) + ggplot2::labs(x = NULL) + th,
      priorpost = ttl(.pp_r(f1), "r: prior and posterior"))
  }
  pl <- stats::setNames(lapply(panels, one), panels)
  if(!requireNamespace("patchwork", quietly = TRUE)) return(pl)
  p <- patchwork::wrap_plots(pl, ncol = ncol)
  if(multi) p <- p + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")
  p
}

#' @rdname jabbaplot3
#' @param ... passed to \code{jabbaplot3()}
#' @export
plot.jabba3 <- function(x, type = c("full", "simple"), ...) jabbaplot3(x, type = type, ...)

# observed indices scaled to biomass: I / (q * selectivity curve) and / Bmsy
.index_scaled <- function(f, rel = FALSE){
  d <- f$setup$dat; rp <- f$report; nm <- f$setup$index.names
  ix <- f$index[is.finite(f$index$obs), c("index", "year", "obs")]
  i <- match(ix$index, nm); qq <- rp$q[d$sets.q[i]]
  sel <- if(isTRUE(d$eb.on)) mapply(function(p, k) .rho_curve(p, d$ebm[k, ]), rp$P[match(ix$year, f$years)], d$ebci[i]) else 1
  ix$value <- ix$obs/(qq*rp$K*sel)*f$refpts["K", "est"]        # I/(q K) is on the B/K scale (q, K internally scaled)
  if(rel) ix$value <- ix$value/f$refpts["Bmsy", "est"]
  ix$index <- factor(ix$index, nm); ix
}

.kobe_runs <- function(fl){
  d <- do.call(rbind, lapply(names(fl), function(n){ f <- fl[[n]]
    data.frame(run = n, year = f$years, stock = f$ts$est[f$ts$quantity == "B/Bmsy"], harvest = f$ts$est[f$ts$quantity == "F/Fmsy"]) }))
  d <- d[is.finite(d$harvest), ]; d$run <- factor(d$run, names(fl))
  ends <- do.call(rbind, lapply(split(d, d$run), function(z) z[nrow(z), ]))
  xm <- max(2, d$stock)*1.05; ym <- max(2, d$harvest)*1.05
  qd <- data.frame(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1), ymax = c(1, 1, Inf, Inf),
                   zone = c("green", "yellow", "orange", "red"))
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
                       fill = c("#2ca25f", "#fee391", "#fd8d3c", "#de2d26"), alpha = 0.85) +
    ggplot2::geom_path(data = d, ggplot2::aes(.data$stock, .data$harvest, colour = .data$run), linewidth = 0.5) +
    ggplot2::geom_point(data = ends, ggplot2::aes(.data$stock, .data$harvest, fill = .data$run), shape = 21, size = 2.5) +
    ggplot2::scale_colour_manual(values = stats::setNames(jb3cols(length(fl)), names(fl)), name = NULL) +
    ggplot2::scale_fill_manual(values = stats::setNames(jb3cols(length(fl)), names(fl)), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, xm), ylim = c(0, ym), expand = FALSE) +
    ggplot2::labs(x = expression(B/B[MSY]), y = expression(F/F[MSY]))
}

.prod_runs <- function(fl){
  cv <- do.call(rbind, lapply(names(fl), function(n){ z <- .spcurve(fl[[n]]); data.frame(run = n, z$curve) }))
  rp <- do.call(rbind, lapply(names(fl), function(n){ z <- .spcurve(fl[[n]]); data.frame(run = n, Bmsy = z$rp$Bmsy, MSY = z$rp$MSY) }))
  cv$run <- factor(cv$run, names(fl)); rp$run <- factor(rp$run, names(fl)); cols <- stats::setNames(jb3cols(length(fl)), names(fl))
  p <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, colour = "grey60")
  if(length(fl) == 1){ z <- .spcurve(fl[[1]])
    p <- p + ggplot2::geom_path(data = z$trj, ggplot2::aes(.data$B, .data$C), colour = "grey40", na.rm = TRUE) +
      ggplot2::geom_point(data = z$trj, ggplot2::aes(.data$B, .data$C), size = 0.6, colour = "grey30", na.rm = TRUE) }
  p + ggplot2::geom_line(data = cv[cv$sp >= 0, ], ggplot2::aes(.data$b, .data$sp, colour = .data$run), linewidth = 0.8) +
    ggplot2::geom_segment(data = rp, ggplot2::aes(x = .data$Bmsy, xend = .data$Bmsy, y = 0, yend = .data$MSY, colour = .data$run), linetype = 2) +
    ggplot2::geom_point(data = rp, ggplot2::aes(.data$Bmsy, .data$MSY, colour = .data$run), size = 2) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Biomass", y = "Production, catch")
}

.pp_r <- function(f){
  pr <- f$priors$r; est <- f$pars["r", c("est", "se.log")]
  xs <- exp(seq(log(est$est) - 4*max(pr[2], est$se.log), log(est$est) + 4*max(pr[2], est$se.log), length.out = 300))
  d <- rbind(data.frame(x = xs, y = stats::dlnorm(xs, log(pr[1]), pr[2]), what = "prior"),
             data.frame(x = xs, y = stats::dlnorm(xs, log(est$est), est$se.log), what = "posterior"))
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y, fill = .data$what)) + ggplot2::geom_area(alpha = 0.5, position = "identity") +
    ggplot2::scale_fill_manual(values = c(prior = "grey60", posterior = "dodgerblue"), name = NULL) +
    ggplot2::labs(x = "r", y = NULL)
}
