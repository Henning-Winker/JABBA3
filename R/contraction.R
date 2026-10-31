# Prior-posterior contraction: how much the data inform each quantity relative to the prior predictive run

.cq_draws <- function(fit, quantities, years){
  d <- fit$draws
  if(is.null(d) || !length(d)) d <- jb3_draws(fit, 1000)
  out <- list()
  for(q in quantities){
    if(q %in% names(d$pars)) out[[q]] <- d$pars[[q]]
    else if(q %in% c("B/Bmsy", "F/Fmsy", "B/K", "B", "F")){
      m <- d[[q]]; yy <- if(is.null(years)) max(fit$years) else years
      for(y in yy) out[[paste(q, y)]] <- m[, as.character(y)]
    } else if(q %in% c("creep") && !is.null(d$aux)) out[[q]] <- d$aux$creep
  }
  out
}

#' Prior-posterior contraction and shift
#'
#' Compares the posterior of a fit with its prior predictive run (\code{jabba3(..., prior.only = TRUE)}, the same
#' model, priors and catch with the index and auxiliary likelihoods switched off) for parameters, reference points
#' and status. For each quantity (log scale for positive quantities):
#' \itemize{
#'   \item contraction = 1 - var(posterior)/var(prior predictive): the share of the prior uncertainty removed by
#'     the data (0 = the data add nothing, 1 = fully determined by the data);
#'   \item shift = (mean(posterior) - mean(prior predictive))/sd(prior predictive): how far the data move the
#'     estimate, in prior SDs.
#' }
#' A quantity with low contraction is prior-driven; one with a large shift is in conflict with its prior (the
#' data pull it into the tail of the prior); high contraction with a moderate shift means it is informed by the
#' data (Schad et al. 2021; Betancourt 2018). A clearly negative contraction (posterior wider than the prior
#' predictive) means the data move the quantity into a region the priors and catch constrain less, and is
#' flagged as "posterior wider". The prior predictive run keeps the catch and the status priors,
#' so contraction measures the information in the indices and auxiliary data beyond the priors and the catch;
#' a catch-only fit therefore has a contraction of zero by construction. For r, K and the initial depletion
#' \code{contraction.pdf} compares the posterior with the stated prior density instead, which also shows how much
#' the catch history (stock survival) trims the prior.
#' @param fit a \code{jabba3} fit, or a named list of fits
#' @param prior the prior predictive fit(s) (same layout as \code{fit}); NULL runs them with NUTS
#' @param quantities parameters (columns of the draws: K, r, MSY, Bmsy, Fmsy, sdb, bkinit, m) and status
#'   quantities (B/Bmsy, F/Fmsy, B/K)
#' @param years years for the status quantities (default the last year)
#' @param contraction.min below this a quantity is classed as prior-driven
#' @param shift.max above this (absolute) a quantity is classed as in conflict with its prior
#' @param nuts NUTS settings for the prior predictive runs
#' @param digits rounding (NULL for none)
#' @return data.frame(run, quantity, prior.median, post.median, contraction, shift, contraction.pdf, class);
#'   the prior predictive fits are attached as attribute "prior"
#' @references Betancourt, M. (2018). Calibrating model-based inferences and decisions. arXiv:1803.08393.
#'   Schad, D.J., Betancourt, M. and Vasishth, S. (2021). Toward a principled Bayesian workflow in cognitive
#'   science. Psychological Methods 26, 103-126.
#' @examples
#' \dontrun{
#' ct <- jb3_contraction(fit)
#' jbplot_contraction(ct)
#' }
#' @export
jb3_contraction <- function(fit, prior = NULL, quantities = c("K", "r", "MSY", "bkinit", "B/Bmsy", "F/Fmsy"),
                            years = NULL, contraction.min = 0.25, shift.max = 2, nuts = list(), digits = 3){
  if(inherits(fit, "jabba3")) fit <- stats::setNames(list(fit), fit$run %||% "run")
  if(!is.null(prior) && inherits(prior, "jabba3")) prior <- list(prior)
  if(is.null(prior)) prior <- lapply(fit, function(f){
    message("Prior predictive run for ", f$run, " (NUTS)")
    nu <- utils::modifyList(list(chains = 3), nuts)
    jabba3(f$data, f$priors, f$model, engine = "nuts", nuts = nu, prior.only = TRUE, run = paste(f$run, "(prior)"))
  })
  stopifnot(length(prior) == length(fit))
  pdf.sd <- function(f, q){ if(!q %in% c("r", "K", "bkinit")) return(NA); p <- f$priors[[q]]; if(!is.null(p) && length(p) >= 2 && (length(p) < 3 || p[3] == 1)) p[2] else NA }
  res <- do.call(rbind, lapply(seq_along(fit), function(k){
    f <- fit[[k]]; po <- .cq_draws(f, quantities, years); pr <- .cq_draws(prior[[k]], quantities, years)
    do.call(rbind, lapply(intersect(names(po), names(pr)), function(q){
      a <- po[[q]]; b <- pr[[q]]; lg <- all(a > 0, na.rm = TRUE) && all(b > 0, na.rm = TRUE)
      if(lg){ a <- log(a); b <- log(b) }
      vb <- stats::var(b, na.rm = TRUE)
      if(!is.finite(vb) || vb < 1e-12) return(NULL)
      ps <- pdf.sd(f, q)
      data.frame(run = names(fit)[k], quantity = q,
                 prior.median = stats::median(if(lg) exp(b) else b, na.rm = TRUE), post.median = stats::median(if(lg) exp(a) else a, na.rm = TRUE),
                 contraction = 1 - stats::var(a, na.rm = TRUE)/vb,
                 shift = (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE))/sqrt(vb),
                 contraction.pdf = if(is.finite(ps) && lg) 1 - stats::var(a, na.rm = TRUE)/ps^2 else NA_real_)
    }))
  }))
  res$class <- ifelse(abs(res$shift) > shift.max, "prior-data conflict",
                ifelse(res$contraction < -contraction.min, "posterior wider",
                ifelse(res$contraction < contraction.min, "prior-driven", "data-informed")))
  if(!is.null(digits)){ k <- vapply(res, is.double, TRUE); res[k] <- lapply(res[k], signif, digits) }
  rownames(res) <- NULL
  attr(res, "prior") <- stats::setNames(prior, names(fit))
  attr(res, "limits") <- c(contraction.min = contraction.min, shift.max = shift.max)
  class(res) <- c("jb3contraction", "data.frame")
  res
}

#' Plot prior-posterior contraction against shift
#'
#' One point per quantity (and run) from \code{\link{jb3_contraction}()}: contraction on the x-axis, absolute
#' shift (in prior SDs) on the y-axis. Shaded regions mark prior-driven quantities (low contraction) and
#' prior-data conflict (large shift).
#' @param x a \code{jb3contraction} table, or a fit / list of fits passed to \code{jb3_contraction()}
#' @param ... passed to \code{jb3_contraction()} when \code{x} is a fit
#' @return a ggplot
#' @export
jbplot_contraction <- function(x, ...){
  if(!inherits(x, "jb3contraction")) x <- jb3_contraction(x, ...)
  lim <- attr(x, "limits") %||% c(contraction.min = 0.25, shift.max = 2)
  d <- as.data.frame(x); d$ashift <- abs(d$shift); xl <- max(-1, min(-0.2, d$contraction, na.rm = TRUE) - 0.05)
  d$contraction <- pmax(pmin(d$contraction, 1), xl)
  d$run <- factor(d$run, unique(d$run))
  ym <- max(3, lim[["shift.max"]] + 1, d$ashift, na.rm = TRUE)
  g <- ggplot2::ggplot(d, ggplot2::aes(contraction, ashift)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = lim[["contraction.min"]], ymin = -Inf, ymax = Inf, fill = "grey90") +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = lim[["shift.max"]], ymax = Inf, fill = "#fde0dd") +
    ggplot2::annotate("text", x = (xl + lim[["contraction.min"]])/2, y = lim[["shift.max"]]*0.93, label = "prior-driven", size = 3, colour = "grey40") +
    ggplot2::annotate("text", x = 0.5, y = ym*0.97, label = "prior-data conflict", size = 3, colour = "#c51b8a") +
    ggplot2::geom_point(ggplot2::aes(colour = run), size = 2.2)
  if(requireNamespace("ggrepel", quietly = TRUE))
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = quantity, colour = run), size = 2.8, show.legend = FALSE,
                                      max.overlaps = 30, seed = 1)
  else g <- g + ggplot2::geom_text(ggplot2::aes(label = quantity, colour = run), size = 2.8, vjust = -0.8, show.legend = FALSE)
  if(nlevels(d$run) > 1) g <- g + ggplot2::facet_wrap(~run)
  g + ggplot2::scale_x_continuous(limits = c(xl, 1), breaks = seq(-1, 1, 0.5)) + ggplot2::coord_cartesian(ylim = c(0, ym)) +
    ggplot2::labs(x = "Posterior contraction (1 - var post / var prior predictive)", y = "|Shift| (prior SDs)", colour = NULL) +
    theme_jabba3() + ggplot2::theme(legend.position = if(nlevels(d$run) > 1) "none" else "bottom")
}

#' @export
print.jb3contraction <- function(x, ...){ y <- x; class(y) <- "data.frame"; attr(y, "prior") <- NULL; print(y, ...); invisible(x) }
