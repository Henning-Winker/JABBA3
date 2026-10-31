# Index comparison tools for trend (JARA) fits: relative catchability, its distribution, and pairwise index correlations

.rel_q_draws <- function(fit, ndraws = 4000){
  # log q draws (reference index = column 1, fixed at 0): NUTS samples if available, else the lognormal approximation
  nI <- length(fit$indices)
  if(!is.null(fit$nuts$stanfit)){
    pm <- as.matrix(fit$nuts$stanfit); cn <- grep("^logq", colnames(pm))
    return(cbind(0, pm[, cn, drop = FALSE]))
  }
  fx <- summary(fit$sdr, "fixed"); lq <- fx[rownames(fx) == "logq", , drop = FALSE]
  set.seed(1)
  cbind(0, sapply(seq_len(nrow(lq)), function(k) stats::rnorm(ndraws, lq[k, 1], lq[k, 2])))
}

#' Catchability relative to the reference index
#'
#' For a relative trend fit (\code{\link{jb3_trend}()}), the catchability of each index relative to the first
#' (reference) index, whose q is 1 by definition: the scaling that aligns each index with the common trend.
#' Estimates and 95\% intervals come from the sdreport of log q (lognormal), or from the NUTS draws for
#' \code{engine = "nuts"}.
#' @param fit a relative \code{jb3trend} fit
#' @param digits rounding (NULL for none)
#' @return data.frame(index, q, lci, uci, log.se, n, years, reference)
#' @examples
#' \dontrun{
#' tr <- jb3_trend(iccat$bet$cpue, iccat$bet$se)
#' jb3_qtable(tr); jbplot_q(tr)
#' }
#' @export
jb3_qtable <- function(fit, digits = NULL){
  stopifnot(inherits(fit, "jb3trend"))
  if(fit$type != "relative") stop("jb3_qtable() needs a relative trend fit (census indices have q = 1)", call. = FALSE)
  nm <- fit$indices; l <- fit$data$long
  if(!is.null(fit$nuts$stanfit)){
    d <- .rel_q_draws(fit); q <- apply(exp(d), 2, stats::quantile, c(0.5, 0.025, 0.975))
    est <- q[1, ]; lci <- q[2, ]; uci <- q[3, ]; se <- apply(d, 2, stats::sd)
  } else {
    fx <- summary(fit$sdr, "fixed"); lq <- fx[rownames(fx) == "logq", , drop = FALSE]
    m <- c(0, lq[, 1]); se <- c(0, lq[, 2])
    est <- exp(m); lci <- exp(m - 1.96*se); uci <- exp(m + 1.96*se)
  }
  se[1] <- NA; lci[1] <- uci[1] <- 1
  yr <- vapply(nm, function(v){ y <- l$year[l$index == v]; paste0(min(y), "-", max(y)) }, "")
  out <- data.frame(index = nm, q = est, lci = lci, uci = uci, log.se = se, n = as.vector(table(factor(l$index, nm))),
                    years = yr, reference = seq_along(nm) == 1, row.names = NULL)
  if(!is.null(digits)) out[c("q", "lci", "uci", "log.se")] <- lapply(out[c("q", "lci", "uci", "log.se")], signif, digits)
  out
}

#' Distribution of catchability relative to the reference index
#'
#' Density of q relative to the first (reference) index for each other index of a relative trend fit, with the
#' reference as a vertical line at 1: the lognormal approximation from the sdreport, or the NUTS draws for
#' \code{engine = "nuts"}. Narrow distributions mean the indices overlap enough to be scaled against each other;
#' wide distributions mean the scaling relies on the continuity of the trend.
#' @param fit a relative \code{jb3trend} fit
#' @param log plot q on a log axis
#' @return a ggplot
#' @export
jbplot_q <- function(fit, log = TRUE){
  stopifnot(inherits(fit, "jb3trend"))
  if(fit$type != "relative") stop("jbplot_q() needs a relative trend fit", call. = FALSE)
  nm <- fit$indices; if(length(nm) < 2) stop("needs at least two indices", call. = FALSE)
  d <- .rel_q_draws(fit)[, -1, drop = FALSE]
  dd <- data.frame(index = factor(rep(nm[-1], each = nrow(d)), nm[-1]), q = exp(as.vector(d)))
  cols <- stats::setNames(jb3cols(length(nm))[seq_along(nm)][-1], nm[-1])
  if(length(nm) == 2) cols[] <- "dodgerblue"
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data$q, fill = .data$index, colour = .data$index)) +
    ggplot2::geom_density(alpha = 0.35, adjust = 1.2) +
    ggplot2::geom_vline(xintercept = 1, linewidth = 0.8) +
    ggplot2::annotate("label", x = 1, y = Inf, label = paste("reference:", nm[1]), vjust = 1.3, hjust = -0.05, size = 3, label.size = 0) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) + ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = paste("q relative to", nm[1]), y = "Density",
                  subtitle = if(!is.null(fit$nuts$stanfit)) "NUTS posterior" else "Lognormal approximation at the mode") +
    theme_jabba3() + ggplot2::theme(legend.position = "right")
  if(log) p <- p + ggplot2::scale_x_log10()
  p
}

#' Pairwise correlation of abundance indices
#'
#' Correlation matrix of the log indices over the years in which each pair overlaps, with p-values (Pearson test) and
#' the number of shared years: a quick check for conflicting indices before fitting a joint index. Pairs with fewer than
#' \code{min.n} shared years are left blank.
#' @param x a \code{jb3trend} fit, or data.frame(year, index1, ...)
#' @param min.n minimum number of shared years
#' @param alpha significance level marked in the plot
#' @return a ggplot; the table (index1, index2, r, p, n) is in attribute \code{"cor"}
#' @export
jbplot_indexcor <- function(x, min.n = 5, alpha = 0.05){
  I <- if(inherits(x, "jb3trend")) x$data$I else as.data.frame(x)
  nm <- names(I)[-1]; L <- log(as.matrix(I[, -1, drop = FALSE])); L[!is.finite(L)] <- NA
  tab <- do.call(rbind, lapply(seq_along(nm), function(i) do.call(rbind, lapply(seq_along(nm), function(j){
    ok <- is.finite(L[, i]) & is.finite(L[, j]); n <- sum(ok)
    r <- p <- NA_real_
    if(i != j && n >= min.n){ ct <- stats::cor.test(L[ok, i], L[ok, j]); r <- unname(ct$estimate); p <- ct$p.value }
    if(i == j){ r <- 1; p <- 0 }
    data.frame(index1 = nm[i], index2 = nm[j], r = r, p = p, n = n) }))))
  tab$index1 <- factor(tab$index1, nm); tab$index2 <- factor(tab$index2, rev(nm))
  low <- tab[match(as.character(tab$index1), nm) <= match(as.character(tab$index2), nm), ]
  low$lab <- ifelse(is.finite(low$r) & low$index1 != low$index2, paste0(sprintf("%.2f", low$r), ifelse(low$p < alpha, "*", ""),
                                                                       "\nn=", low$n), "")
  p <- ggplot2::ggplot(low, ggplot2::aes(.data$index1, .data$index2, fill = .data$r)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 2.6, lineheight = 0.85) +
    ggplot2::scale_fill_gradient2(low = "#de2d26", mid = "white", high = "#3182bd", midpoint = 0, limits = c(-1, 1),
                                  na.value = "grey92", name = "r") +
    ggplot2::labs(x = NULL, y = NULL, subtitle = sprintf("Correlation of log indices over shared years (* p < %.2f)", alpha)) +
    theme_jabba3() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "right",
                                    panel.grid = ggplot2::element_blank())
  attr(p, "cor") <- tab[as.character(tab$index1) != as.character(tab$index2), c("index1", "index2", "r", "p", "n")]
  p
}
