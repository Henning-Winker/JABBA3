# Summary tables: parameters, time series of quantities, projections and the Kobe strategy matrix

.fitlist0 <- function(x){
  if(inherits(x, "jabba3")) x <- list(x)
  if(is.null(names(x))) names(x) <- vapply(x, function(f) f$run %||% "run", "")
  names(x) <- make.unique(names(x)); x
}

.fmt_ci <- function(est, lci, uci, digits){
  f <- function(v) ifelse(is.finite(v), formatC(signif(v, digits), format = "fg", digits = digits, big.mark = ","), "")
  ifelse(is.finite(lci) & is.finite(uci), paste0(f(est), " (", f(lci), "-", f(uci), ")"), f(est))
}

#' Parameter and reference point table
#'
#' Estimates with 95\% intervals (posterior median and credible interval for NUTS fits, mode and delta-method interval
#' otherwise) for the model parameters, reference points and the terminal status.
#' @param x a \code{jabba3} fit or a (named) list of fits
#' @param digits significant digits for the formatted table
#' @param format "long" (one row per run and parameter, numeric columns) or "wide" (one column per run, "est (lci-uci)")
#' @return data.frame
#' @examples
#' \dontrun{
#' jb3_partable(fit)
#' jb3_partable(list(Fox = fit, Lag.AR1 = rec), format = "wide")
#' }
#' @export
jb3_partable <- function(x, digits = 3, format = c("long", "wide")){
  format <- match.arg(format); fl <- .fitlist0(x)
  tab <- do.call(rbind, lapply(names(fl), function(n){
    f <- fl[[n]]; p <- f$pars; rp <- f$refpts
    row <- function(nm, src, lab = nm) if(nm %in% rownames(src)) data.frame(par = lab, est = src[nm, "est"], lci = src[nm, "lci"],
                                                                            uci = src[nm, "uci"], cv = src[nm, "se.log"]) else NULL
    st <- function(q, lab){ t <- f$ts[f$ts$quantity == q & is.finite(f$ts$est), ]; t <- t[nrow(t), ]
      data.frame(par = paste0(lab, " (", t$year, ")"), est = t$est, lci = t$lci, uci = t$uci, cv = t$se.log) }
    out <- rbind(row("K", p), row("r", p), row("m", p),
                 data.frame(par = "Bmsy/K", est = rp["BmsyK", "est"], lci = NA, uci = NA, cv = NA),
                 row("Bmsy", rp), row("Fmsy", rp), row("MSY", rp), row("bkinit", p, "B1/K"),
                 row("sdb", p, "sigma.proc"),
                 if("rho" %in% rownames(p)) row("rho", p) else data.frame(par = "rho", est = f$report$rho %||% 0, lci = NA, uci = NA, cv = NA),
                 row("kdyn", p),
                 do.call(rbind, lapply(grep("^q[0-9]+$", rownames(p), value = TRUE), row, src = p)),
                 do.call(rbind, lapply(grep("^sdi[0-9]+$", rownames(p), value = TRUE), row, src = p)),
                 do.call(rbind, lapply(grep("^(qA|sdA)\\.", rownames(p), value = TRUE), row, src = p)), row("creep", p),
                 st("B/Bmsy", "B/Bmsy"), st("F/Fmsy", "F/Fmsy"), st("B/K", "B/K"))
    cbind(run = n, out)
  }))
  rownames(tab) <- NULL
  if(format == "long") return(tab)
  tab$txt <- .fmt_ci(tab$est, tab$lci, tab$uci, digits)
  tab$par0 <- sub(" \\(.*", "", tab$par)                           # align terminal-year rows across runs
  w <- stats::reshape(tab[, c("run", "par0", "txt")], idvar = "par0", timevar = "run", direction = "wide")
  names(w) <- c("Parameter", names(fl)); rownames(w) <- NULL
  w[is.na(w)] <- ""; w
}

#' Time series table of estimated quantities
#'
#' @param fit a \code{jabba3} fit
#' @param quantities any of "Catch", "B", "F", "B/Bmsy", "F/Fmsy", "B/K", "Surplus production", "Process deviation",
#'   "Process SD"
#' @param ci include the 95\% interval columns
#' @param years optional subset of years
#' @return data.frame, one row per year and columns \code{<quantity>}, \code{<quantity>.lci}, \code{<quantity>.uci}
#' @export
jb3_tstable <- function(fit, quantities = c("Catch", "B", "F", "B/Bmsy", "F/Fmsy", "B/K"), ci = TRUE, years = NULL){
  yrs <- if(is.null(years)) fit$years else intersect(years, fit$years)
  out <- data.frame(year = yrs)
  for(q in quantities){
    t <- fit$ts[fit$ts$quantity == q, ]; i <- match(yrs, t$year)
    if(!nrow(t)) { warning("quantity not found: ", q, call. = FALSE); next }
    out[[q]] <- t$est[i]
    if(ci){ out[[paste0(q, ".lci")]] <- t$lci[i]; out[[paste0(q, ".uci")]] <- t$uci[i] }
  }
  out
}

#' Projection tables
#'
#' \code{jb3_projtable()}: catch options table by scenario, in the style of ICES advice sheets: F/Fmsy and F, the catch in
#' the first projection year, median B/Bmsy in selected years with the change from the terminal assessment year, and the
#' probabilities of B >= Bmsy, F <= Fmsy and of the Kobe green quadrant.
#'
#' \code{jb3_kobematrix()}: Kobe II strategy matrix (ICCAT): probability by scenario (rows) and projection year (columns) of
#' being in the green quadrant (B >= Bmsy and F <= Fmsy), or of B >= Bmsy or F <= Fmsy alone.
#' @param proj output of \code{\link{jb3_project}()}
#' @param years projection years to report (calendar years); default the first, fifth and last projection year
#' @param digits rounding
#' @param blim limit as a fraction of Bmsy for the risk P(B < Blim) (default 0.3, SPiCT benchmark default; NULL to omit)
#' @return data.frame
#' @examples
#' \dontrun{
#' pj <- jb3_project(fit, fmsy = c(0, 0.5, 0.75, 1, 1.25))
#' jb3_projtable(pj)
#' jb3_kobematrix(pj)
#' }
#' @export
jb3_projtable <- function(proj, years = NULL, digits = 3, blim = 0.3){
  py <- proj$iters[[1]]$years
  if(is.null(years)) years <- unique(py[c(1, min(5, length(py)), length(py))])
  bb0 <- stats::median(proj$status$`B/Bmsy`)
  out <- do.call(rbind, lapply(seq_along(proj$iters), function(s){
    it <- proj$iters[[s]]; sc <- proj$scenarios[s, ]; j <- match(years, it$years)
    r <- data.frame(scenario = sc$label, `F/Fmsy` = stats::median(it$`F/Fmsy`[, 1]), F = stats::median(it$F[, 1]),
                    Catch = stats::median(it$Catch[, 1]), check.names = FALSE)
    for(k in seq_along(j)){ y <- years[k]; b <- it$`B/Bmsy`[, j[k]]; h <- it$`F/Fmsy`[, j[k]]
      r[[paste0("B/Bmsy ", y)]] <- stats::median(b)
      r[[paste0("dB% ", y)]] <- 100*(stats::median(b)/bb0 - 1)
      r[[paste0("P(B>=Bmsy) ", y)]] <- mean(b >= 1)
      if(!is.null(blim)) r[[paste0("P(B<Blim) ", y)]] <- mean(b < blim)
      r[[paste0("P(F<=Fmsy) ", y)]] <- mean(h <= 1 + 1e-8)
      r[[paste0("P(green) ", y)]] <- mean(b >= 1 & h <= 1 + 1e-8) }
    r }))
  names(out)[names(out) == "Catch"] <- paste("Catch", py[1])
  num <- vapply(out, is.numeric, TRUE); out[num] <- lapply(out[num], signif, digits)
  rownames(out) <- NULL; out
}

#' Kobe II strategy matrix
#'
#' Kobe II strategy matrix (K2SM) in the ICCAT format: the probability, by scenario (rows, e.g. TAC levels) and projection
#' year (columns), that F < Fmsy, that B > Bmsy, and of both (the green Kobe quadrant). Several projections, e.g. from
#' alternative models or a JABBA and an SS3 base case, can be combined with equal or given weights, as in the ICCAT
#' combined matrices.
#'
#' Probabilities are the proportions of projection draws in each condition; the conditions are evaluated as
#' B/Bmsy >= 1 and F/Fmsy <= 1 (with a small tolerance, so a scenario at exactly Fmsy counts as F <= Fmsy), which for
#' continuous draws is the same as the strict inequalities of the ICCAT tables.
#' @param proj output of \code{\link{jb3_project}()} or \code{\link{jb3_fwd}()}, or a (named) list of them with the same
#'   scenarios and years, to be combined
#' @param what \code{"green"} (F < Fmsy and B > Bmsy), \code{"B"} (B > Bmsy), \code{"F"} (F < Fmsy), or \code{"all"} for
#'   the three matrices (a, b, c of the ICCAT tables) as a list
#' @param years projection years to show (default all)
#' @param percent if TRUE, probabilities in percent, rounded to integers (ICCAT style)
#' @param weights weights for combining a list of projections (default equal)
#' @param digits rounding of probabilities when \code{percent = FALSE}
#' @return a data.frame with the scenario value in the first column (\code{TAC} for constant catch scenarios,
#'   \code{F/Fmsy} for F scenarios, otherwise \code{scenario}) and one column per year; for \code{what = "all"} a named list
#'   of three such data.frames, of class \code{jb3k2sm}
#' @seealso \code{\link{jbplot_kobematrix}()}, \code{\link{jb3_projtable}()}
#' @examples
#' \dontrun{
#' pc <- jb3_project(fit, catch = seq(1000, 3500, 250), nyears = 10)
#' jb3_kobematrix(pc, percent = TRUE)                       # P(green)
#' k2 <- jb3_kobematrix(pc, what = "all", percent = TRUE)   # ICCAT a), b), c)
#' k2$`F<Fmsy`
#' jbplot_kobematrix(pc)
#' # combined matrix of two models, equal weight
#' jb3_kobematrix(list(base = pc1, alt = pc2), what = "all", percent = TRUE)
#' }
#' @export
jb3_kobematrix <- function(proj, what = c("green", "B", "F", "all"), years = NULL, percent = FALSE, weights = NULL,
                           digits = 2){
  what <- match.arg(what)
  pl <- if(inherits(proj, "jb3proj")) list(proj) else proj
  if(!all(vapply(pl, inherits, TRUE, "jb3proj"))) stop("proj must be a jb3_project() output or a list of them", call. = FALSE)
  p1 <- pl[[1]]; yrs <- p1$iters[[1]]$years; sc <- p1$scenarios
  for(p in pl[-1]){
    if(!identical(p$scenarios$label, sc$label) || !identical(p$iters[[1]]$years, yrs))
      stop("projections to combine must have the same scenarios and years", call. = FALSE)
  }
  if(is.null(weights)) weights <- rep(1, length(pl))
  if(length(weights) != length(pl)) stop("weights must have one value per projection", call. = FALSE)
  weights <- weights/sum(weights)
  if(is.null(years)) years <- yrs
  j <- match(years, yrs)
  if(anyNA(j)) stop("years not in the projection: ", paste(years[is.na(j)], collapse = ", "), call. = FALSE)
  # first column: the scenario value
  type <- unique(sc$type)
  key <- if(identical(type, "Catch")) "TAC" else if(identical(type, "F/Fmsy")) "F/Fmsy" else "scenario"
  lab <- if(key == "scenario") sc$label else sc$value
  one <- function(w){
    m <- Reduce(`+`, Map(function(p, wt) wt*t(vapply(p$iters, function(it){
      b <- it$`B/Bmsy`[, j, drop = FALSE]; h <- it$`F/Fmsy`[, j, drop = FALSE]
      colMeans(switch(w, green = b >= 1 & h <= 1 + 1e-8, B = b >= 1, F = h <= 1 + 1e-8))
    }, numeric(length(j)))), pl, weights))
    m <- matrix(m, nrow = nrow(sc))
    m <- if(percent) round(100*m) else round(m, digits)
    out <- data.frame(lab, m, check.names = FALSE, stringsAsFactors = FALSE)
    names(out) <- c(key, years); rownames(out) <- NULL
    out
  }
  if(what != "all") return(one(what))
  structure(list(`F<Fmsy` = one("F"), `B>Bmsy` = one("B"), `F<Fmsy & B>Bmsy` = one("green")),
            class = c("jb3k2sm", "list"), percent = percent)
}

#' @export
print.jb3k2sm <- function(x, ...){
  ttl <- c(`F<Fmsy` = "a) Probability that F < Fmsy", `B>Bmsy` = "b) Probability that B > Bmsy",
           `F<Fmsy & B>Bmsy` = "c) Probability that F < Fmsy and B > Bmsy")
  for(n in names(x)){ cat(ttl[n], if(isTRUE(attr(x, "percent"))) "(%)", "\n"); print(x[[n]], row.names = FALSE); cat("\n") }
  invisible(x)
}

.k2lab <- function(v) if(is.numeric(v)) format(v, scientific = FALSE, trim = TRUE, drop0trailing = TRUE) else as.character(v)

#' Plot a Kobe II strategy matrix
#'
#' Shaded Kobe II strategy matrices (ICCAT style): probabilities by scenario (rows) and year (columns) for F < Fmsy,
#' B > Bmsy and both, with the probabilities printed in the cells.
#' @param x output of \code{\link{jb3_project}()} (or a list of them, combined as in \code{\link{jb3_kobematrix}()}), or of
#'   \code{jb3_kobematrix(..., what = "all")}
#' @param what which matrices: any of \code{"F"}, \code{"B"}, \code{"green"} (default all three)
#' @param years projection years (default all)
#' @param weights weights for combining a list of projections
#' @param colours \code{"grey"} (ICCAT report tables), \code{"kobe"} (red to green), or a vector of colours, one per
#'   probability class
#' @param breaks probability classes in percent (default 0, 50, 60, 70, 80, 90, 100)
#' @param size text size of the cell labels
#' @param pct if TRUE (default) cell labels with a percent sign, as in the ICCAT tables
#' @param ylab row header (default "Catch (t)" for TAC scenarios)
#' @return a ggplot
#' @examples
#' \dontrun{
#' pc <- jb3_project(fit, catch = seq(1000, 3500, 250), nyears = 10)
#' jbplot_kobematrix(pc)
#' jbplot_kobematrix(pc, what = "green", colours = "kobe")
#' }
#' @export
jbplot_kobematrix <- function(x, what = c("F", "B", "green"), years = NULL, weights = NULL, colours = "grey",
                              breaks = c(0, 50, 60, 70, 80, 90, 100), size = 3, pct = TRUE, ylab = NULL){
  if(inherits(x, "jb3k2iotc")) return(.plot_k2iotc(x, colours = colours, size = size))
  what <- match.arg(what, c("F", "B", "green"), several.ok = TRUE)
  k2 <- if(inherits(x, "jb3k2sm")) x else jb3_kobematrix(x, what = "all", years = years, percent = TRUE, weights = weights)
  if(!isTRUE(attr(k2, "percent"))) k2[] <- lapply(k2, function(d){ d[-1] <- lapply(d[-1], function(v) round(100*v)); d })
  nm <- c(F = "F<Fmsy", B = "B>Bmsy", green = "F<Fmsy & B>Bmsy")[what]
  ttl <- c(`F<Fmsy` = "P(F < Fmsy)", `B>Bmsy` = "P(B > Bmsy)", `F<Fmsy & B>Bmsy` = "P(F < Fmsy and B > Bmsy)")
  key <- names(k2[[1]])[1]
  df <- do.call(rbind, lapply(nm, function(n){
    d <- k2[[n]]; yr <- names(d)[-1]
    data.frame(row = rep(.k2lab(d[[1]]), length(yr)), year = rep(yr, each = nrow(d)), p = unlist(d[-1]),
               panel = ttl[n], stringsAsFactors = FALSE)
  }))
  df$row <- factor(df$row, levels = rev(unique(.k2lab(k2[[1]][[1]]))))
  df$year <- factor(df$year, levels = names(k2[[1]])[-1])
  df$panel <- factor(df$panel, levels = ttl[nm])
  nb <- length(breaks) - 1
  df$class <- cut(df$p, breaks = breaks, include.lowest = TRUE, right = FALSE)
  levels(df$class) <- paste0(breaks[-length(breaks)], "-", breaks[-1], "%")
  df$class <- factor(df$class, levels = levels(df$class))
  pal <- if(identical(colours, "grey")) grDevices::grey(seq(1, 0.35, length.out = nb))
         else if(identical(colours, "kobe")) grDevices::colorRampPalette(c("#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641"))(nb)
         else rep_len(colours, nb)
  names(pal) <- levels(df$class)
  dark <- if(identical(colours, "grey")) names(pal)[seq_len(nb) > 0.6*nb] else character(0)
  ggplot2::ggplot(df, ggplot2::aes(x = year, y = row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = class), colour = "grey60") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(p, if(pct) "%"), colour = class %in% dark), size = size, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~panel, ncol = 1) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = ylab %||% if(key == "TAC") "Catch (t)" else key) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"), legend.position = "right")
}

#' Kobe II strategy matrix, IOTC format
#'
#' Kobe II strategy matrix in the format of the IOTC stock assessment summaries: for each catch scenario (columns,
#' labelled as a percentage of a reference catch and in tonnes) the probability (\%) of violating the MSY-based reference
#' points, P(B < Bmsy) and P(F > Fmsy), and the limit reference points, P(B < Blim) and P(F > Flim), in selected
#' projection years (by default 3 and 10 years after the last assessment year).
#' @param proj output of \code{\link{jb3_project}()} (constant catch scenarios), or a list of them to combine
#' @param years projection years (default the last assessment year + 3 and + 10, or the first and last projection year)
#' @param ref.catch reference catch for the column labels (default the catch in the last assessment year)
#' @param blim,flim limit reference points relative to Bmsy and Fmsy (IOTC defaults for tropical tunas differ by stock,
#'   e.g. Blim = 0.4 or 0.5 Bmsy, Flim = 1.3 or 1.4 Fmsy)
#' @param weights weights for combining a list of projections (the "weighted probability" of the IOTC tables)
#' @param digits decimals of the percentages
#' @return a data.frame of class \code{jb3k2iotc}: columns \code{block} ("MSY" or "Limit"), \code{refpt}, \code{year} and
#'   one column per scenario
#' @seealso \code{\link{jb3_kobematrix}()} (ICCAT format), \code{\link{jbplot_kobematrix}()}
#' @examples
#' \dontrun{
#' cref <- tail(fit$data$catch, 1)
#' pc <- jb3_project(fit, catch = cref*c(0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2), nyears = 10)
#' k2 <- jb3_kobematrix_iotc(pc, ref.catch = cref)
#' k2
#' jbplot_kobematrix(k2)                       # grey
#' jbplot_kobematrix(k2, colours = "kobe")     # colour
#' }
#' @export
jb3_kobematrix_iotc <- function(proj, years = NULL, ref.catch = NULL, blim = 0.5, flim = 1.3, weights = NULL, digits = 1){
  pl <- if(inherits(proj, "jb3proj")) list(proj) else proj
  p1 <- pl[[1]]; yrs <- p1$iters[[1]]$years; sc <- p1$scenarios
  for(p in pl[-1]) if(!identical(p$scenarios$label, sc$label) || !identical(p$iters[[1]]$years, yrs))
    stop("projections to combine must have the same scenarios and years", call. = FALSE)
  w <- if(is.null(weights)) rep(1, length(pl)) else weights; w <- w/sum(w)
  if(is.null(years)){ years <- p1$status$year + c(3, 10); if(!all(years %in% yrs)) years <- range(yrs) }
  j <- match(years, yrs); if(anyNA(j)) stop("years not in the projection: ", paste(years[is.na(j)], collapse = ", "), call. = FALSE)
  if(is.null(ref.catch)){
    ct <- p1$fit.ts[p1$fit.ts$quantity == "Catch", ]; ref.catch <- ct$est[ct$year == p1$status$year]
  }
  prob <- function(fun, q) 100*Reduce(`+`, Map(function(p, wt) wt*vapply(p$iters, function(it) fun(it[[q]][, j, drop = FALSE]),
                                                                            numeric(length(j))), pl, w))
  rows <- list(list("MSY", "B < Bmsy", "B/Bmsy", function(m) colMeans(m < 1)),
               list("MSY", "F > Fmsy", "F/Fmsy", function(m) colMeans(m > 1 + 1e-8)),
               list("Limit", "B < Blim", "B/Bmsy", function(m) colMeans(m < blim)),
               list("Limit", "F > Flim", "F/Fmsy", function(m) colMeans(m > flim)))
  out <- do.call(rbind, lapply(rows, function(r){
    m <- matrix(prob(r[[4]], r[[3]]), nrow = length(j))          # years x scenarios
    data.frame(block = r[[1]], refpt = r[[2]], year = years, round(m, digits), check.names = FALSE)
  }))
  out <- out[order(out$block != "MSY", out$year), ]
  lab <- if(all(sc$type == "Catch")){
    paste0(round(100*sc$value/ref.catch), "% (", format(round(sc$value), big.mark = ",", scientific = FALSE, trim = TRUE), "t)")
  } else sc$label
  names(out)[-(1:3)] <- lab; rownames(out) <- NULL
  structure(out, class = c("jb3k2iotc", "data.frame"), blim = blim, flim = flim, ref.catch = ref.catch)
}

#' @export
print.jb3k2iotc <- function(x, ...){
  d <- as.data.frame(unclass(x), check.names = FALSE, stringsAsFactors = FALSE)
  rp <- paste0(sub(" .*", "", d$refpt), d$year, sub("^[BF]", "", d$refpt))
  cat("Probability (%) of violating MSY-based reference points\n")
  m <- d$block == "MSY"; print(data.frame(`Reference point` = rp[m], d[m, -(1:3)], check.names = FALSE), row.names = FALSE)
  cat("\nProbability (%) of violating limit reference points (Blim = ", attr(x, "blim"), " Bmsy; Flim = ",
      attr(x, "flim"), " Fmsy)\n", sep = "")
  print(data.frame(`Reference point` = rp[!m], d[!m, -(1:3)], check.names = FALSE), row.names = FALSE)
  invisible(x)
}

.plot_k2iotc <- function(x, colours = "grey", size = 3.5, breaks = c(0, 10, 25, 50, 75, 100)){
  d <- as.data.frame(unclass(x), check.names = FALSE, stringsAsFactors = FALSE)
  sc <- names(d)[-(1:3)]
  rp <- paste0(sub(" .*", "", d$refpt), "[", d$year, "]", sub("^[BF]", "", d$refpt))
  df <- data.frame(row = rep(rp, length(sc)), scenario = rep(sc, each = nrow(d)), p = unlist(d[sc]),
                   block = rep(ifelse(d$block == "MSY", "MSY-based reference points",
                                      paste0("Limit reference points (Blim = ", attr(x, "blim"), " Bmsy, Flim = ", attr(x, "flim"), " Fmsy)")),
                               length(sc)), stringsAsFactors = FALSE)
  df$row <- factor(df$row, levels = rev(unique(rp)))
  df$scenario <- factor(sub(" (", "\n(", df$scenario, fixed = TRUE), levels = sub(" (", "\n(", sc, fixed = TRUE))
  df$block <- factor(df$block, levels = unique(df$block))
  nb <- length(breaks) - 1
  df$class <- cut(df$p, breaks = breaks, include.lowest = TRUE, right = FALSE)
  levels(df$class) <- paste0(breaks[-length(breaks)], "-", breaks[-1], "%")
  pal <- if(identical(colours, "grey")) grDevices::grey(seq(1, 0.35, length.out = nb))
         else if(identical(colours, "kobe")) grDevices::colorRampPalette(c("#1a9641", "#a6d96a", "#ffffbf", "#fdae61", "#d7191c"))(nb)
         else rep_len(colours, nb)
  names(pal) <- levels(df$class)
  dark <- if(identical(colours, "grey")) names(pal)[seq_len(nb) > 0.6*nb] else character(0)
  lbl <- function(l) parse(text = gsub(" ", "~", gsub("([<>])", "~'\\1'~", gsub("(msy|lim)", "[\\U\\1]", l, perl = TRUE))))
  ggplot2::ggplot(df, ggplot2::aes(x = scenario, y = row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = class), colour = "grey60") +
    ggplot2::geom_text(ggplot2::aes(label = format(p, nsmall = 1), colour = class %in% dark), size = size, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = "P(violation)") +
    ggplot2::facet_wrap(~block, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0), labels = function(l) parse(text = sub("^([BF])\\[(\\d+)\\](.*)$", "\\1[\\2]~'\\3'", l))) +
    ggplot2::labs(x = paste0("Catch relative to ", format(round(attr(x, "ref.catch")), big.mark = ","), " t"), y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, face = "bold"))
}
