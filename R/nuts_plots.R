# NUTS (tmbstan) convergence diagnostics: summary table and ggplot2 graphics (no bayesplot dependency)

.nuts_check <- function(fit){
  if(is.null(fit$nuts$stanfit)) stop("no NUTS draws: fit with jabba3(..., engine = \"nuts\")", call. = FALSE)
  .need("rstan", "NUTS diagnostics need rstan")
  fit$nuts$stanfit
}

# sampled parameters: fixed effects by default (states logP and catch deviations are many and are covered by jbplot_rhat())
.nuts_pars <- function(sf, pars = NULL, states = FALSE){
  nm <- setdiff(names(sf), "lp__")
  if(!is.null(pars)){
    hit <- unlist(lapply(pars, function(p) nm[nm == p | startsWith(nm, paste0(p, "["))]))
    if(!length(hit)) stop("none of pars found; available: ", paste(unique(sub("\\[.*", "", nm)), collapse = ", "), call. = FALSE)
    return(unique(hit))
  }
  if(states) nm else nm[!grepl("^(logP|zP|cdev)\\[", nm)]
}

# iterations x chains x parameters (post warmup unless inc_warmup)
.nuts_array <- function(sf, pars, inc_warmup = FALSE) rstan::extract(sf, pars = pars, permuted = FALSE, inc_warmup = inc_warmup)

.nuts_long <- function(sf, pars, inc_warmup = FALSE){
  a <- .nuts_array(sf, pars, inc_warmup); di <- dim(a)
  data.frame(iter = rep(seq_len(di[1]), di[2]*di[3]), chain = factor(rep(rep(seq_len(di[2]), each = di[1]), di[3])),
             par = factor(rep(dimnames(a)[[3]], each = di[1]*di[2]), levels = dimnames(a)[[3]]), value = as.vector(a))
}

.nuts_div <- function(sf, inc_warmup = FALSE){
  sp <- rstan::get_sampler_params(sf, inc_warmup = inc_warmup)
  unlist(lapply(sp, function(x) x[, "divergent__"]))
}

#' NUTS convergence summary
#'
#' Rank-normalised split R-hat, bulk and tail effective sample sizes (Vehtari et al. 2021) for every sampled parameter,
#' and the sampler diagnostics: divergent transitions, maximum tree depth hits and E-BFMI by chain.
#' @param fit a \code{jabba3} fit with \code{engine = "nuts"}
#' @param states include the log-depletion states and catch deviations
#' @return list(pars = data.frame(par, mean, sd, q2.5, q50, q97.5, Rhat, ess_bulk, ess_tail), sampler = data.frame by chain,
#'   ok = logical verdict: R-hat < 1.01, bulk and tail ESS > 100 per chain, no divergences, E-BFMI > 0.3)
#' @export
jb3_nuts_diags <- function(fit, states = TRUE){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, states = states); a <- .nuts_array(sf, pars)
  st <- do.call(rbind, lapply(pars, function(p){ m <- a[, , p, drop = TRUE]; if(is.null(dim(m))) m <- matrix(m, ncol = 1)
    data.frame(par = p, mean = mean(m), sd = stats::sd(m), q2.5 = stats::quantile(m, 0.025, names = FALSE),
               q50 = stats::median(m), q97.5 = stats::quantile(m, 0.975, names = FALSE),
               Rhat = rstan::Rhat(m), ess_bulk = rstan::ess_bulk(m), ess_tail = rstan::ess_tail(m)) }))
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE); mtd <- fit$nuts$settings$control$max_treedepth
  smp <- data.frame(chain = seq_along(sp),
                    divergent = sapply(sp, function(x) sum(x[, "divergent__"])),
                    treedepth.max = sapply(sp, function(x) sum(x[, "treedepth__"] >= mtd)),
                    accept = sapply(sp, function(x) mean(x[, "accept_stat__"])),
                    stepsize = sapply(sp, function(x) x[1, "stepsize__"]),
                    EBFMI = sapply(sp, function(x){ e <- x[, "energy__"]; sum(diff(e)^2)/sum((e - mean(e))^2) }))
  nch <- length(sp)
  ok <- all(st$Rhat < 1.01, na.rm = TRUE) && all(st$ess_bulk > 100*nch, na.rm = TRUE) && all(st$ess_tail > 100*nch, na.rm = TRUE) &&
    sum(smp$divergent) == 0 && all(smp$EBFMI > 0.3)
  rownames(st) <- NULL
  list(pars = st, sampler = smp, ok = ok)
}

#' NUTS diagnostic plots
#'
#' Convergence graphics for a fit with \code{engine = "nuts"}, on the sampled (unconstrained, mostly log) scale.
#' \describe{
#'   \item{trace}{traces by chain; divergent transitions as rug ticks}
#'   \item{rank}{rank histograms by chain (Vehtari et al. 2021): flat and overlapping when chains mix}
#'   \item{dens}{posterior densities by chain}
#'   \item{acf}{autocorrelation by chain}
#'   \item{pairs}{pairwise scatter with divergent transitions in red (funnels, ridges, non-identifiability)}
#'   \item{rhat}{R-hat and bulk ESS for every sampled parameter, states included}
#'   \item{energy}{marginal energy against energy transitions with E-BFMI by chain}
#' }
#' \code{jbplot_nutsdiags()} combines trace, rank, R-hat/ESS and energy on one page (needs patchwork).
#' @param fit a \code{jabba3} fit with \code{engine = "nuts"}
#' @param type plot type, see Details
#' @param pars parameter names (e.g. \code{c("logK", "logr")}; \code{"logq"} matches all \code{logq[i]}). Default: all
#'   fixed effects
#' @param inc_warmup include warmup iterations (trace only)
#' @param lags maximum lag for \code{acf}
#' @param npoints maximum points per panel for \code{pairs}
#' @return a ggplot
#' @examples
#' \dontrun{
#' fn <- jabba3(d, pr, jb3_model(shape = "fox"), engine = "nuts", nuts = list(chains = 3, iter = 2000, cores = 3))
#' jb3_nuts_diags(fn)$sampler
#' jbplot_nuts(fn, "trace"); jbplot_nuts(fn, "pairs", pars = c("logK", "logr", "logsdb"))
#' jbplot_nutsdiags(fn)
#' }
#' @export
jbplot_nuts <- function(fit, type = c("trace", "rank", "dens", "acf", "pairs", "rhat", "energy"), pars = NULL,
                        inc_warmup = FALSE, lags = 20, npoints = 1500){
  type <- match.arg(type)
  switch(type, trace = jbplot_trace(fit, pars, inc_warmup), rank = jbplot_rank(fit, pars), dens = jbplot_mcmcdens(fit, pars),
         acf = jbplot_mcmcacf(fit, pars, lags), pairs = jbplot_pairs(fit, pars, npoints), rhat = jbplot_rhat(fit),
         energy = jbplot_energy(fit))
}

.chaincols <- function(n) ggplot2::scale_colour_manual(values = jb3cols(n), name = "Chain")

#' @rdname jbplot_nuts
#' @export
jbplot_trace <- function(fit, pars = NULL, inc_warmup = FALSE){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, pars); x <- .nuts_long(sf, pars, inc_warmup)
  nch <- nlevels(x$chain); div <- .nuts_div(sf, inc_warmup); niter <- length(div)/nch
  dv <- data.frame(iter = rep(seq_len(niter), nch), d = div); dv <- unique(dv$iter[dv$d == 1])
  wu <- if(inc_warmup) fit$nuts$settings$warmup/fit$nuts$settings$thin else NA
  p <- ggplot2::ggplot(x, ggplot2::aes(iter, value, colour = chain)) +
    (if(!is.na(wu)) ggplot2::annotate("rect", xmin = 0, xmax = wu, ymin = -Inf, ymax = Inf, fill = "grey90")) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.8) +
    ggplot2::facet_wrap(~par, scales = "free_y") + .chaincols(nch) +
    ggplot2::labs(x = if(inc_warmup) "Iteration (warmup shaded)" else "Iteration (post warmup)", y = NULL,
                  subtitle = sprintf("%d chains; %d divergent transitions", nch, sum(div))) + theme_jabba3(9)
  if(length(dv)) p <- p + ggplot2::geom_rug(data = data.frame(iter = dv), ggplot2::aes(x = iter), inherit.aes = FALSE,
                                            sides = "b", colour = "red", alpha = 0.5)
  p
}

#' @rdname jbplot_nuts
#' @param bins number of rank bins
#' @export
jbplot_rank <- function(fit, pars = NULL, bins = 20){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, pars); x <- .nuts_long(sf, pars)
  x$rank <- stats::ave(x$value, x$par, FUN = function(v) rank(v, ties.method = "average"))
  n <- max(x$rank); nch <- nlevels(x$chain)
  x$bin <- pmin(floor((x$rank - 1)/n*bins), bins - 1)
  h <- stats::aggregate(list(count = x$value), x[c("par", "chain", "bin")], length)
  ggplot2::ggplot(h, ggplot2::aes(bin, count, colour = chain)) +
    ggplot2::geom_hline(yintercept = n/nch/bins, linetype = 2, colour = "grey50") +
    ggplot2::geom_step(direction = "mid", linewidth = 0.4) +
    ggplot2::facet_wrap(~par) + .chaincols(nch) + ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "Rank bin (pooled chains)", y = "Draws", subtitle = "Rank histograms: flat and overlapping when chains mix") +
    theme_jabba3(9) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' @rdname jbplot_nuts
#' @export
jbplot_mcmcdens <- function(fit, pars = NULL){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, pars); x <- .nuts_long(sf, pars)
  ggplot2::ggplot(x, ggplot2::aes(value, colour = chain)) + ggplot2::geom_density(linewidth = 0.4) +
    ggplot2::facet_wrap(~par, scales = "free") + .chaincols(nlevels(x$chain)) +
    ggplot2::labs(x = NULL, y = "Density") + theme_jabba3(9) + ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

#' @rdname jbplot_nuts
#' @export
jbplot_mcmcacf <- function(fit, pars = NULL, lags = 20){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, pars); a <- .nuts_array(sf, pars)
  x <- do.call(rbind, lapply(dimnames(a)[[3]], function(p) do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch){
    v <- a[, ch, p]; ac <- if(stats::var(v) > 0) stats::acf(v, lag.max = lags, plot = FALSE)$acf[, 1, 1] else rep(NA, lags + 1)
    data.frame(par = p, chain = factor(ch), lag = 0:lags, acf = ac) }))))
  x$par <- factor(x$par, levels = dimnames(a)[[3]])
  ggplot2::ggplot(x, ggplot2::aes(lag, acf, colour = chain)) + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 0.6) +
    ggplot2::facet_wrap(~par) + .chaincols(dim(a)[2]) + ggplot2::coord_cartesian(ylim = c(-0.2, 1)) +
    ggplot2::labs(x = "Lag", y = "Autocorrelation") + theme_jabba3(9)
}

#' @rdname jbplot_nuts
#' @export
jbplot_pairs <- function(fit, pars = NULL, npoints = 1500){
  sf <- .nuts_check(fit); pars <- .nuts_pars(sf, pars)
  if(length(pars) > 8) message("pairs: ", length(pars), " parameters; consider pars = c(...) for a readable plot")
  a <- .nuts_array(sf, pars); m <- matrix(a, ncol = dim(a)[3], dimnames = list(NULL, dimnames(a)[[3]]))
  div <- .nuts_div(sf) == 1
  id <- if(nrow(m) > npoints) sort(unique(c(sample.int(nrow(m), npoints), which(div)))) else seq_len(nrow(m))
  m <- m[id, , drop = FALSE]; div <- div[id]; pn <- colnames(m)
  sc <- do.call(rbind, lapply(seq_along(pn), function(i) do.call(rbind, lapply(seq_along(pn), function(j){
    if(i <= j) return(NULL)
    data.frame(px = factor(pn[j], pn), py = factor(pn[i], pn), x = m[, j], y = m[, i], div = div) }))))
  # diagonal: histogram rescaled to the panel's y range so it shares the facet scale
  dg <- do.call(rbind, lapply(seq_along(pn), function(i){
    h <- graphics::hist(m[, i], breaks = 25, plot = FALSE); r <- range(m[, i])
    data.frame(px = factor(pn[i], pn), py = factor(pn[i], pn), xmin = utils::head(h$breaks, -1), xmax = h$breaks[-1],
               ymin = r[1], ymax = r[1] + h$counts/max(h$counts)*diff(r)) }))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = dg, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha = 0.6) +
    ggplot2::geom_point(data = sc[!sc$div, ], ggplot2::aes(x, y), size = 0.3, alpha = 0.3, colour = "grey30") +
    ggplot2::facet_grid(py ~ px, scales = "free", switch = "both") +
    ggplot2::labs(x = NULL, y = NULL, subtitle = sprintf("%d divergent transitions (red)", sum(.nuts_div(sf)))) +
    theme_jabba3(8) + ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank(),
                                     panel.spacing = ggplot2::unit(0.15, "lines"), axis.text = ggplot2::element_text(size = 5))
  if(any(sc$div)) p <- p + ggplot2::geom_point(data = sc[sc$div, ], ggplot2::aes(x, y), size = 0.8, colour = "red")
  p
}

#' @rdname jbplot_nuts
#' @export
jbplot_rhat <- function(fit){
  dg <- jb3_nuts_diags(fit, states = TRUE); x <- dg$pars; nch <- nrow(dg$sampler)
  x$group <- ifelse(grepl("^(logP|zP)\\[", x$par), "states (logP or zP)", ifelse(grepl("^cdev\\[", x$par), "catch deviations", "fixed effects"))
  x <- rbind(data.frame(par = x$par, group = x$group, metric = "R-hat", value = x$Rhat),
             data.frame(par = x$par, group = x$group, metric = "Bulk ESS / chain", value = x$ess_bulk/nch),
             data.frame(par = x$par, group = x$group, metric = "Tail ESS / chain", value = x$ess_tail/nch))
  x$metric <- factor(x$metric, c("R-hat", "Bulk ESS / chain", "Tail ESS / chain"))
  x$flag <- with(x, ifelse(metric == "R-hat", value > 1.01, value < 100))
  ref <- data.frame(metric = factor(c("R-hat", "R-hat", "Bulk ESS / chain", "Tail ESS / chain"), levels(x$metric)),
                    v = c(1.01, 1.05, 100, 100), lt = c(2, 1, 2, 2))
  x$rk <- stats::ave(x$value, x$group, x$metric, FUN = function(v) rank(v, ties.method = "first"))
  lab <- x[x$flag & x$group == "fixed effects", ]
  ggplot2::ggplot(x, ggplot2::aes(value, rk, colour = flag)) +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = v, linetype = factor(lt)), colour = "grey40", show.legend = FALSE) +
    ggplot2::geom_point(size = 0.9) +
    (if(nrow(lab)) ggplot2::geom_text(data = lab, ggplot2::aes(label = par), size = 2, hjust = -0.15, show.legend = FALSE)) +
    ggplot2::facet_grid(group ~ metric, scales = "free") +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "dodgerblue", `TRUE` = "red"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, subtitle = sprintf("max R-hat %.3f; min bulk ESS %.0f; %s", max(dg$pars$Rhat, na.rm = TRUE),
                                                         min(dg$pars$ess_bulk, na.rm = TRUE), if(dg$ok) "converged" else "check")) +
    ggplot2::labs(y = "Parameters (sorted)") +
    theme_jabba3(8) + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
}

#' @rdname jbplot_nuts
#' @export
jbplot_energy <- function(fit){
  sf <- .nuts_check(fit); sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  x <- do.call(rbind, lapply(seq_along(sp), function(ch){ e <- sp[[ch]][, "energy__"]
    bf <- sum(diff(e)^2)/sum((e - mean(e))^2)
    lab <- sprintf("Chain %d: E-BFMI %.2f", ch, bf)
    rbind(data.frame(chain = lab, what = "Marginal energy (centred)", value = e - mean(e)),
          data.frame(chain = lab, what = "Energy transition", value = c(diff(e), NA))) }))
  ggplot2::ggplot(x[is.finite(x$value), ], ggplot2::aes(value, fill = what)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, position = "identity", alpha = 0.5) +
    ggplot2::facet_wrap(~chain) + ggplot2::scale_fill_manual(values = c("red", "dodgerblue")) +
    ggplot2::labs(x = "Energy", y = NULL, subtitle = "Transitions narrower than the marginal (E-BFMI < 0.3) signal poor exploration") +
    theme_jabba3(9) + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank())
}

#' @rdname jbplot_nuts
#' @export
jbplot_nutsdiags <- function(fit, pars = NULL){
  .need("patchwork", "jbplot_nutsdiags() needs patchwork")
  g <- function(p) p + ggplot2::theme(legend.position = "none")
  (g(jbplot_trace(fit, pars)) | g(jbplot_rank(fit, pars))) / (jbplot_rhat(fit) | jbplot_energy(fit)) +
    patchwork::plot_layout(heights = c(1, 1))
}

#' Effective sample size of the derived states
#'
#' Bulk and tail ESS (rank-normalised, Vehtari et al. 2021) and R-hat of the derived trajectories (e.g. B/Bmsy and
#' F/Fmsy by year) of a NUTS fit, also per 1000 post-warmup draws. Useful to compare parameterisations
#' (\code{jb3_model(param = "centred")} against \code{"noncentred"}), because the sampled states differ between them.
#' @param fit a \code{jabba3} fit with \code{engine = "nuts"}
#' @param quantities derived quantities, any of "B/Bmsy", "F/Fmsy", "B/K", "B", "F"
#' @return data.frame(quantity, year, Rhat, ess_bulk, ess_tail, bulk.per1000, tail.per1000); attribute "summary" with
#'   minimum and median by quantity and the run time and parameterisation
#' @export
jb3_nuts_ess <- function(fit, quantities = c("B/Bmsy", "F/Fmsy")){
  sf <- .nuts_check(fit); nch <- fit$nuts$settings$chains; dr <- fit$draws
  nd <- nrow(dr[[quantities[1]]]); ni <- nd/nch
  out <- do.call(rbind, lapply(quantities, function(q){ x <- dr[[q]]
    do.call(rbind, lapply(seq_len(ncol(x)), function(j){ m <- matrix(x[, j], ni, nch)
      data.frame(quantity = q, year = fit$years[j], Rhat = rstan::Rhat(m), ess_bulk = rstan::ess_bulk(m), ess_tail = rstan::ess_tail(m)) })) }))
  out$bulk.per1000 <- 1000*out$ess_bulk/nd; out$tail.per1000 <- 1000*out$ess_tail/nd
  sm <- do.call(rbind, lapply(split(out, out$quantity), function(z) data.frame(quantity = z$quantity[1],
    min.bulk.per1000 = min(z$bulk.per1000), median.bulk.per1000 = stats::median(z$bulk.per1000),
    min.tail.per1000 = min(z$tail.per1000), max.Rhat = max(z$Rhat))))
  attr(out, "summary") <- cbind(param = fit$nuts$param %||% "centred", sm, divergent = fit$nuts$divergent, time = fit$time)
  out
}
