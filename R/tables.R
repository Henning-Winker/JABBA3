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

#' @rdname jb3_projtable
#' @param what "green" (B >= Bmsy and F <= Fmsy), "B" (B >= Bmsy) or "F" (F <= Fmsy)
#' @export
jb3_kobematrix <- function(proj, what = c("green", "B", "F"), digits = 2){
  what <- match.arg(what)
  m <- t(vapply(proj$iters, function(it){
    b <- it$`B/Bmsy`; h <- it$`F/Fmsy`
    colMeans(switch(what, green = b >= 1 & h <= 1 + 1e-8, B = b >= 1, F = h <= 1 + 1e-8))
  }, numeric(length(proj$iters[[1]]$years))))
  out <- data.frame(scenario = proj$scenarios$label, round(m, digits), check.names = FALSE)
  names(out)[-1] <- proj$iters[[1]]$years
  out
}
