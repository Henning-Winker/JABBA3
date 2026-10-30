# Internal helpers

`%||%` <- function(a, b) if(is.null(a)) b else a

# smooth bounds: exact away from the bound, differentiable at it (needed for a well-conditioned Laplace Hessian)
.smax <- function(x, a, eps = 1e-8) a + 0.5*((x - a) + sqrt((x - a)^2 + eps))
.smin <- function(x, a, eps = 1e-8) a - 0.5*((a - x) + sqrt((a - x)^2 + eps))
.clamp <- function(x, lo, hi, eps = 1e-5) .smin(.smax(x, lo, eps), hi, eps)
# clamp of a positive quantity on the log scale (smoothing width ~1%)
.lclamp <- function(x, lo, hi, eps = 1e-4) exp(.smin(.smax(log(x), log(lo), eps), log(hi), eps))

# log-scale sd from a CV
cv2sd <- function(cv) sqrt(log(cv^2 + 1))

# Pella-Tomlinson shape m from Bmsy/K
.m_from_bmsyk <- function(bmsyk){
  if(abs(bmsyk - 0.5) < 1e-8) return(2)
  if(abs(bmsyk - exp(-1)) < 1e-4) return(1.001)                              # Fox
  if(bmsyk < exp(-1)) return(stats::uniroot(function(m) m^(-1/(m - 1)) - bmsyk, c(0.01, 0.9999), tol = 1e-10)$root)
  stats::uniroot(function(m) m^(-1/(m - 1)) - bmsyk, c(1.0001, 50), tol = 1e-10)$root
}

# EB/SB ratio curve rho(P) = u1 + (u2 - u1) * (1 - exp(-u3 (P - P1))) / (1 - exp(-u3 (P2 - P1))), P clamped to [P1, P2]
.rho_curve <- function(p, e){
  pc <- pmin(pmax(p, e[4]), e[5])
  e[1] + (e[2] - e[1])*(1 - exp(-e[3]*(pc - e[4])))/(1 - exp(-e[3]*(e[5] - e[4])))
}

#' Runs test with three-sigma control limits
#'
#' Wald-Wolfowitz runs test on the signs of residuals (normal approximation) and three-sigma limits from the average
#' moving range (Montgomery), as in JABBA \code{jbruns_sig3()}, ss3diags and JAGG.
#' @param x residuals
#' @param mixing "less" (default, one-sided: too few runs), "two.sided" or "greater"
#' @return list with \code{sig3lim} (lcl, ucl) and \code{p.runs}
#' @export
jb3_runs <- function(x, mixing = c("less", "two.sided", "greater")){
  mixing <- match.arg(mixing)
  x <- x[is.finite(x)]
  if(length(x) < 3) return(list(sig3lim = c(NA, NA), p.runs = NA))
  mr <- abs(diff(x)); amr <- mean(mr); mr <- mr[mr < 3.267*amr]
  stdev <- mean(mr)/1.128; lim <- c(-3*stdev, 3*stdev)
  s <- sign(x); s <- s[s != 0]; n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
  if(n1 == 0 || n2 == 0) return(list(sig3lim = lim, p.runs = 0.001))
  runs <- 1 + sum(diff(s) != 0)
  mu <- 2*n1*n2/n + 1; v <- 2*n1*n2*(2*n1*n2 - n)/(n^2*(n - 1))
  z <- (runs - mu)/sqrt(v)
  p <- switch(mixing, less = stats::pnorm(z), greater = 1 - stats::pnorm(z), two.sided = 2*min(stats::pnorm(z), 1 - stats::pnorm(z)))
  list(sig3lim = lim, p.runs = round(p, 3))
}

.need <- function(pkg, why) if(!requireNamespace(pkg, quietly = TRUE)) stop(why, ": install.packages('", pkg, "')", call. = FALSE)
