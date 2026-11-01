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

#' Observation and process noise of indices
#'
#' Three measures of how noisy each index is: the mean input CV (sampling or observation error, from the log-SEs),
#' the process-error CV of Francis (2011), from the residuals of a loess smoother fitted to the log index over time
#' (the year-to-year noise the index carries beyond a smooth trend), and, for a trend fit, the estimated additional
#' SD \eqn{\sigma_{p,i}} of each index against the common trend (Conn 2010) with the RMSE of its residuals.
#' A loess CV much larger than the input CV, or a large \eqn{\sigma_{p,i}}, points to variation in catchability or
#' availability that the sampling error does not capture; a survey index estimated with a design-based estimator
#' typically carries more such noise than one standardised with a (geostatistical) model (Winker et al. 2017).
#' @param x a \code{jb3trend} fit, or a data.frame of indices (year in the first column)
#' @param se data.frame of log-SEs in the layout of \code{x} (when \code{x} is a data.frame)
#' @param span loess span (Francis 2011 smoother)
#' @param digits rounding (NULL for none)
#' @return data.frame(index, n, years, obs.cv, loess.cv, and for a trend fit sd.extra, sd.total, rmse)
#' @references Francis, R.I.C.C. (2011) Data weighting in statistical fisheries stock assessment models. CJFAS 68, 1124-1138.
#' @examples
#' \dontrun{
#' tr <- jb3_trend(iccat$bum$cpue, iccat$bum$se)
#' jb3_noise(tr); jbplot_noise(tr)
#' }
#' @export
jb3_noise <- function(x, se = NULL, span = 0.6, digits = 3){
  if(inherits(x, "jb3trend")){ I <- x$data$I; S <- x$data$SE } else { I <- as.data.frame(x); S <- if(is.null(se)) NULL else as.data.frame(se) }
  yrs <- I[[1]]; nm <- names(I)[-1]
  out <- do.call(rbind, lapply(nm, function(v){
    ok <- is.finite(I[[v]]) & I[[v]] > 0; y <- yrs[ok]; li <- log(I[[v]][ok])
    lcv <- if(sum(ok) >= 6) tryCatch({ f <- stats::loess(li ~ y, span = max(span, 4/sum(ok)), degree = 1)
                                      stats::sd(stats::residuals(f)) }, error = function(e) NA_real_) else NA_real_
    ocv <- if(!is.null(S) && v %in% names(S)) mean(S[[v]][ok], na.rm = TRUE) else NA_real_
    data.frame(index = v, n = sum(ok), years = paste0(min(y), "-", max(y)), obs.cv = ocv, loess.cv = lcv)
  }))
  if(inherits(x, "jb3trend")){
    l <- x$index; sdi <- x$sdi
    out$sd.extra <- sdi
    out$sd.total <- sqrt(out$obs.cv^2 + sdi^2)
    out$rmse <- vapply(nm, function(v){ r <- l$residual[l$index == v]; sqrt(mean(r^2, na.rm = TRUE)) }, 0)
  }
  rownames(out) <- NULL
  if(!is.null(digits)){ k <- vapply(out, is.double, TRUE); out[k] <- lapply(out[k], round, digits) }
  out
}

#' Plot observation and process noise of indices
#'
#' Bars of the measures of \code{\link{jb3_noise}()} by index: input observation CV, Francis (2011) loess process
#' CV and, for a trend fit, the estimated additional SD against the common trend.
#' @inheritParams jb3_noise
#' @return a ggplot
#' @export
jbplot_noise <- function(x, se = NULL, span = 0.6){
  z <- jb3_noise(x, se, span, digits = NULL)
  ms <- c(obs.cv = "Input observation CV", loess.cv = "Loess process CV (Francis 2011)", sd.extra = "Additional SD vs common trend")
  ms <- ms[names(ms) %in% names(z)]
  d <- do.call(rbind, lapply(names(ms), function(m) data.frame(index = z$index, measure = ms[[m]], value = z[[m]])))
  d$index <- factor(d$index, z$index); d$measure <- factor(d$measure, ms)
  ggplot2::ggplot(d, ggplot2::aes(index, value, fill = measure)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.75, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("grey70", "#2b8cbe", "#d95f02")[seq_along(ms)]) +
    ggplot2::labs(x = NULL, y = "CV / log-scale SD", fill = NULL) + theme_jabba3() +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Split an index at a change in vessel, gear or method
#'
#' Splits one index into two series at a known change in catchability (a new survey vessel, a new trawl gear, a
#' change in the fishery), so that a trend or production model estimates a separate q for each period. Years
#' fished with both set-ups (calibration or paired tows) can be given to both series; they identify the q ratio
#' directly, while without overlap the ratio rests on the continuity of the trend across the change.
#' @param I data.frame of indices (year in the first column)
#' @param se data.frame of log-SEs in the layout of \code{I} (optional)
#' @param index name of the index to split
#' @param at first year of the new set-up
#' @param overlap years in which the index is measured with both set-ups: the values go to both series (for
#'   calibration years with separate values, supply them in \code{new})
#' @param new optional values of the new set-up in the \code{overlap} years (same length)
#' @param labels names of the two series
#' @return list(I, se) with the split index replaced by two columns
#' @export
jb3_splitindex <- function(I, se = NULL, index, at, overlap = NULL, new = NULL, labels = paste0(index, c(".old", ".new"))){
  I <- as.data.frame(I); yr <- I[[1]]; v <- I[[index]]
  a <- ifelse(yr < at | yr %in% overlap, v, NA); b <- ifelse(yr >= at, v, NA)
  if(!is.null(new)){ stopifnot(length(new) == length(overlap)); b[match(overlap, yr)] <- new }
  pos <- match(index, names(I))
  ins <- function(df, x1, x2){ x <- stats::setNames(data.frame(x1, x2), labels)
    cbind(df[seq_len(pos - 1)], x, df[-seq_len(pos)]) }
  out <- list(I = ins(I, a, b))
  if(!is.null(se)){ se <- as.data.frame(se); s <- se[[index]]; out$se <- ins(se, ifelse(is.na(a), NA, s), ifelse(is.na(b), NA, s)) }
  out
}
