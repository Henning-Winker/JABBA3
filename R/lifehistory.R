# Life-history conversions: Leslie-matrix r and generation time, r <-> Fmsy for the production function shape, and
# Selex EB/SB curves from length-based selectivity. Inputs follow FLicc (lhpar, selpar on the length scale); the Leslie
# matrix follows JABBA::jbleslie() and FLSRTMB::productivity(), the selectivity functions FLicc::sel_curve_flicc()
# (H. Winker, FLicc), ported here so that FLicc is not required.

# ---- input helpers ---------------------------------------------------------------------------------------------------

# one sex from an FLicc-style lhpar: FLPar (optionally with a sex/unit dimension "F", "M"), named numeric or list
.lh_one <- function(lhpar, sex = "F"){
  if(inherits(lhpar, "FLPar")){
    dn <- dimnames(lhpar); sd <- setdiff(names(dn), c("params", "iter"))
    a <- methods::as(lhpar, "array")
    if(length(sd)){
      s <- dn[[sd[1]]]; j <- if(sex %in% s) match(sex, s) else 1
      v <- apply(a, c(1, which(names(dn) == sd[1])), function(z) stats::median(z, na.rm = TRUE))[, j]
    } else v <- apply(a, 1, function(z) stats::median(z, na.rm = TRUE))
    stats::setNames(as.numeric(v), dn$params)
  } else if(is.list(lhpar) && !is.null(lhpar[[sex]])) .lh_one(lhpar[[sex]], sex) else unlist(lhpar)
}
.sexes <- function(lhpar){
  if(inherits(lhpar, "FLPar")){ dn <- dimnames(lhpar); sd <- setdiff(names(dn), c("params", "iter"))
    if(length(sd)) return(dn[[sd[1]]]) }
  if(is.list(lhpar) && all(c("F", "M") %in% names(lhpar))) return(c("F", "M"))
  "F"
}
.lhget <- function(p, nm, default = NA_real_) if(nm %in% names(p) && is.finite(p[[nm]])) as.numeric(p[[nm]]) else default

# age-based biology from lhpar: von Bertalanffy length, weight W = a L^b, maturity logistic at length (L50, L95),
# constant M (or Mk * k); ages amin:amax
.lh_age <- function(p, amin = 0, amax = NULL, t0 = NULL){
  linf <- .lhget(p, "linf"); k <- .lhget(p, "k"); M <- .lhget(p, "M")
  if(!is.finite(M) && is.finite(.lhget(p, "Mk"))) M <- .lhget(p, "Mk")*k
  if(!is.finite(linf) || !is.finite(k) || !is.finite(M)) stop("lhpar needs linf, k and M (or Mk)", call. = FALSE)
  t0 <- t0 %||% .lhget(p, "t0", -0.5)
  L50 <- .lhget(p, "L50", 0.66*linf); L95 <- .lhget(p, "L95", L50 + 0.2*linf)
  a <- .lhget(p, "a", 0.01); b <- .lhget(p, "b", 3)
  amax <- amax %||% .lhget(p, "amax", ceiling(-log(0.01)/M))            # age at which 1% of a cohort survives
  age <- amin:amax
  len <- linf*(1 - exp(-k*(age - t0))); len[len < 0.01*linf] <- 0.01*linf
  list(age = age, len = len, wt = a*len^b, mat = 1/(1 + exp(-log(19)*(len - L50)/(L95 - L50))), M = rep(M, length(age)),
       linf = linf, k = k, t0 = t0, L50 = L50, L95 = L95, a = a, b = b, h = .lhget(p, "s", .lhget(p, "h", 0.7)))
}

# ---- Leslie core (JABBA::jbleslie with steepness) --------------------------------------------------------------------
.leslie_core <- function(age, M, wt, mat, h, plus_group = TRUE){
  n <- length(age); lx <- numeric(n); lx[1] <- 1
  if(n > 1) for(a in 2:n) lx[a] <- lx[a - 1]*exp(-M[a - 1])
  if(plus_group) lx[n] <- lx[n]/max(1e-12, 1 - exp(-M[n]))
  spr0a <- lx*mat*wt; spr0 <- sum(spr0a)
  rs <- 4*h/(spr0*(1 - h))                                             # recruits per spawner at the origin (B-H)
  L <- matrix(0, n, n); L[1, ] <- rs*wt*mat
  if(n > 1) for(a in 2:n) L[a, a - 1] <- exp(-M[a - 1])
  if(plus_group) L[n, n] <- exp(-M[n])
  lambda <- max(Re(eigen(L, only.values = TRUE)$values))
  c(r = log(lambda), gt = sum(age*spr0a)/spr0, lambda = lambda, spr0 = spr0)
}

#' Intrinsic rate of increase and generation time from a Leslie matrix
#'
#' The maximum intrinsic rate of increase \eqn{r = \log\lambda} (dominant eigenvalue of the Leslie matrix with recruitment
#' at the slope at the origin of the Beverton-Holt stock-recruitment relationship, \eqn{4h/(\phi_0(1-h))}) and the
#' generation time GT (mean age of spawners in the unfished population), as \code{JABBA::jbleslie()} and
#' \code{FLSRTMB::productivity()}.
#'
#' Inputs:
#' \describe{
#'   \item{FLicc \code{lhpar}}{an \code{FLPar} (or named vector or list) with \code{linf}, \code{k}, \code{M} (or
#'     \code{Mk}), \code{L50}, \code{L95}, \code{a}, \code{b} and optionally \code{t0} (default -0.5), \code{s} (steepness,
#'     default 0.7) and \code{amax}; biology at age from von Bertalanffy growth, \eqn{W = aL^b} and logistic maturity at
#'     length. Sex-structured: an \code{FLPar} with a second dimension (e.g. \code{sex = c("F", "M")}), or
#'     \code{list(F = , M = )}; reproduction from \code{sex}.}
#'   \item{\code{FLStock}}{\code{m}, \code{mat}, \code{stock.wt} by age (median over years and iterations; sex from the
#'     \code{unit} dimension if it has \code{"F"}).}
#'   \item{\code{FLQuants}}{\code{m}, \code{mat}, \code{wt} (or \code{stock.wt}) by age with \code{unit = c("F", "M")}
#'     for sexes, or by length (first dimension \code{len}) together with \code{lhpar} for growth (length at age).}
#'   \item{list or data.frame}{\code{age}, \code{M}, \code{mat}, \code{wt} (as \code{JABBA::jbio()}).}
#' }
#' @param x lhpar, FLStock, FLQuants, or list/data.frame of biology at age
#' @param h steepness (default from \code{lhpar$s}, else 0.7)
#' @param sex sex that determines reproduction (default \code{"F"})
#' @param lhpar growth parameters for length-based FLQuants
#' @param amin,amax age range for lhpar (default 0 to the age at which 1\% survive)
#' @param plus_group treat the last age as a plus group
#' @return named numeric: r, gt, lambda, spr0 (unfished spawning biomass per recruit)
#' @examples
#' lh <- c(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.2, a = 0.0047, b = 3.15, s = 0.75)
#' jb3_leslie(lh)
#' @export
jb3_leslie <- function(x, h = NULL, sex = "F", lhpar = NULL, amin = 0, amax = NULL, plus_group = TRUE){
  if(inherits(x, "FLStock")){
    q <- list(m = FLCore::m(x), mat = FLCore::mat(x), wt = FLCore::stock.wt(x))
    return(.leslie_flq(q, h %||% 0.7, sex, plus_group, ages = as.numeric(dimnames(q$m)[[1]])))
  }
  if(inherits(x, "FLQuants")){
    q <- as.list(x); if(is.null(q$wt)) q$wt <- q$stock.wt
    if(names(dimnames(q$m))[1] == "len"){
      if(is.null(lhpar)) stop("length-based FLQuants need lhpar (growth) for the length at age", call. = FALSE)
      g <- .lh_age(.lh_one(lhpar, sex), amin, amax)
      at <- function(f){ v <- .flq_vec(f, sex); L <- as.numeric(dimnames(f)[[1]]); stats::approx(L, v, g$len, rule = 2)$y }
      return(.leslie_core(g$age, at(q$m), at(q$wt), at(q$mat), h %||% g$h, plus_group))
    }
    return(.leslie_flq(q, h %||% 0.7, sex, plus_group, ages = as.numeric(dimnames(q$m)[[1]])))
  }
  if(is.data.frame(x) || (is.list(x) && !is.null(x$age))){
    return(.leslie_core(as.numeric(x$age), rep_len(as.numeric(x$M), length(x$age)), as.numeric(x$wt %||% x$Wa),
                        as.numeric(x$mat %||% x$Mat), h %||% x$h %||% 0.7, plus_group))
  }
  g <- .lh_age(.lh_one(x, sex), amin, amax)
  .leslie_core(g$age, g$M, g$wt, g$mat, h %||% g$h, plus_group)
}

.flq_vec <- function(f, sex = "F"){
  a <- methods::as(f, "array"); u <- dimnames(f)$unit
  if(!is.null(u) && sex %in% u) a <- a[, , match(sex, u), , , , drop = FALSE] else a <- a[, , 1, , , , drop = FALSE]
  apply(a, 1, stats::median, na.rm = TRUE)
}
.leslie_flq <- function(q, h, sex, plus_group, ages){
  .leslie_core(ages, .flq_vec(q$m, sex), .flq_vec(q$wt, sex), .flq_vec(q$mat, sex), h, plus_group)
}

#' Convert between Fmsy and the surplus production r for a given shape
#'
#' In the Pella-Tomlinson form used by JABBA3, \eqn{F_{MSY} = r/m} with m the shape implied by \eqn{B_{MSY}/K}, so the
#' r prior that reproduces a given harvest rate at MSY is \eqn{r = m F_{MSY}} (as \code{JABBA::r_from_fmsy_bmsyk()}).
#' Works for \eqn{B_{MSY}/K < 1/e} (m < 1). A Leslie r is the per-capita growth at low biomass, which equals
#' \eqn{r/(m-1)} only for the Schaefer model; matching \eqn{F_{MSY}} (from an age-structured equilibrium, or
#' \eqn{r_{Leslie}/2} as the Schaefer equivalent) keeps the r prior consistent across shapes.
#' @param Fmsy harvest rate at MSY
#' @param r intrinsic rate of the production function
#' @param BmsyK Bmsy/K (default 0.4); ignored when \code{shape} is "schaefer" or "fox"
#' @param shape "pella" (from \code{BmsyK}), "schaefer" or "fox"
#' @return r (\code{jb3_fmsy2r}) or Fmsy (\code{jb3_r2fmsy})
#' @examples
#' jb3_fmsy2r(0.2, BmsyK = 0.4)       # r = 1.19 x 0.2
#' jb3_fmsy2r(0.2, shape = "schaefer")
#' jb3_r2fmsy(0.3, BmsyK = 0.36)
#' @export
jb3_fmsy2r <- function(Fmsy, BmsyK = 0.4, shape = c("pella", "schaefer", "fox")){
  Fmsy*.shape_m(BmsyK, match.arg(shape))
}
#' @rdname jb3_fmsy2r
#' @export
jb3_r2fmsy <- function(r, BmsyK = 0.4, shape = c("pella", "schaefer", "fox")){
  r/.shape_m(BmsyK, match.arg(shape))
}
.shape_m <- function(BmsyK, shape) switch(shape, schaefer = 2, fox = 1.001, pella = sapply(BmsyK, .m_from_bmsyk))

# ---- selectivity (FLicc::sel_curve_flicc) -------------------------------------------------------------------------------
.sel_fun_from <- function(p){
  n <- names(p)
  if(all(c("SL50", "SL95", "SLdesc", "CVdesc", "Smin") %in% n)) "jabba5" else if(all(c("SL50", "SL95") %in% n)) "logistic" else if(all(c("mode", "lsd", "rsd") %in% n)) "dsnormal" else
    if(all(c("mode1", "sd1", "mode2", "lsd2", "rsd2", "h") %in% n)) "bidsnormal" else
      if(all(c("mode1", "sd1", "mode2", "sd2", "h") %in% n)) "binormal" else if(all(c("mode", "sd") %in% n)) "normal" else
        stop("selpar names not recognised (FLicc: SL50/SL95, mode/lsd/rsd, mode/sd, mode1/sd1/mode2/sd2/h, ",
             "mode1/sd1/mode2/lsd2/rsd2/h; JABBA 5-parameter: SL50/SL95/SLdesc/CVdesc/Smin)", call. = FALSE)
}
.sel_len <- function(L, sel_fun, p){
  p <- as.list(p)
  switch(sel_fun,
    jabba5 = jb3_sel5(L, unlist(p)),
    logistic = 1/(1 + exp(-log(19)*(L - p$SL50)/(p$SL95 - p$SL50))),
    dsnormal = { d <- L - p$mode; exp(-d^2/(2*ifelse(L < p$mode, p$lsd, p$rsd)^2)) },
    normal = exp(-(L - p$mode)^2/(2*p$sd^2)),
    binormal = { s <- exp(-(L - p$mode1)^2/(2*p$sd1^2)) + p$h*exp(-(L - p$mode2)^2/(2*p$sd2^2)); s/max(s) },
    bidsnormal = { d2 <- L - p$mode2
      s <- exp(-(L - p$mode1)^2/(2*p$sd1^2)) + p$h*exp(-d2^2/(2*ifelse(L < p$mode2, p$lsd2, p$rsd2)^2)); s/max(s) })
}
# selectivity at age: integrate selectivity at length over a normal length-at-age distribution (CV of length at age)
.sel_age <- function(len, sel_fun, p, CVL = 0.1, linf){
  if(CVL <= 0) return(.sel_len(len, sel_fun, p))
  z <- seq(-3, 3, length.out = 41); w <- stats::dnorm(z); w <- w/sum(w)
  sapply(len, function(l){ L <- pmax(l*(1 + CVL*z), 0.001*linf); sum(w*.sel_len(L, sel_fun, p)) })
}
#' JABBA 5-parameter selectivity
#'
#' The JABBA-Select selectivity (parameters as \code{JABBA::jselex()}): logistic ascending limb with 50\% and 95\%
#' selection at \code{SL50} and \code{SL95}; beyond \code{SLdesc} a normal descending limb (SD = \code{SLdesc} x
#' \code{CVdesc}) declining from 1 towards a minimum \code{Smin} (\code{Smin = 1}, or a large \code{SLdesc}, gives
#' logistic selectivity). Note: \code{JABBA::jselex()} as coded inverts the descending limb (selectivity drops to
#' \code{Smin} just beyond \code{SLdesc} and rises again); here the descending limb declines continuously. Usable on length or age (\code{x} and the parameters on the same scale).
#' @param x lengths or ages
#' @param pars c(SL50, SL95, SLdesc, CVdesc, Smin), named or in this order
#' @param type "length" or "age" (documentation of the scale; the function is the same)
#' @return selectivity in [0, 1]
#' @examples
#' jb3_sel5(seq(10, 100, 10), c(SL50 = 35, SL95 = 45, SLdesc = 70, CVdesc = 0.15, Smin = 0.05))
#' jb3_sel5(0:15, c(2, 3, 8, 0.2, 0.1), type = "age")
#' @export
jb3_sel5 <- function(x, pars, type = c("length", "age")){
  type <- match.arg(type)
  p <- as.numeric(pars); nm <- names(pars)
  if(!is.null(nm) && all(c("SL50", "SL95", "SLdesc", "CVdesc", "Smin") %in% nm)) p <- as.numeric(pars[c("SL50", "SL95", "SLdesc", "CVdesc", "Smin")])
  if(length(p) < 5) stop("pars must be c(SL50, SL95, SLdesc, CVdesc, Smin)", call. = FALSE)
  up <- 1/(1 + exp(-log(19)*(x - p[1])/(p[2] - p[1])))
  dn <- p[5] + (1 - p[5])*exp(-(x - p[3])^2/(2*(p[3]*p[4])^2))         # 1 at SLdesc, declining towards Smin
  s <- ifelse(x <= p[3], up, dn); pmin(pmax(s, 0), 1)
}

.selpar_list <- function(selpar){
  if(inherits(selpar, "FLPars") || (is.list(selpar) && !is.numeric(selpar) && !inherits(selpar, "FLPar") &&
                                     !all(sapply(selpar, length) == 1))){
    s <- lapply(as.list(selpar), function(z) .lh_one(z))
  } else s <- list(.lh_one(selpar))
  if(is.null(names(s)) || any(names(s) == "")) names(s) <- paste0("fleet", seq_along(s))
  s
}

#' Selex inputs from length-based selectivity: EB/SB curves for the spawning-biomass state
#'
#' Computes, for each fleet (selectivity), the equilibrium ratio of the exploitable stock to spawning biomass as a
#' function of spawning-biomass depletion P = SB/SB0 over a grid of F, and fits the saturating curve
#' \eqn{EB/SB(P) = u_1 + (u_2 - u_1)(1 - e^{-u_3(P - P_1)})/(1 - e^{-u_3(P_2 - P_1)})} used by
#' \code{\link{jb3_model}(ebsb = )} (as \code{JABBA::make_ebsb()} and \code{fit_ebsb()}). The exploitable stock is
#' measured in weight (\code{units = "biomass"}, EB/SB, for catches and indices in weight) or in numbers
#' (\code{units = "numbers"}, exploitable numbers per unit SB, EN/SB, for indices in numbers, as \code{I.unit = 0} in
#' JABBA-Select). Because fishing shifts the exploitable stock towards smaller fish, EN/SB rises more steeply with
#' depletion than EB/SB, so an index in numbers declines less than one in weight; mapping it with the weight curve
#' makes the model optimistic. The curve is fitted over \code{fit.range} of SB/SB0 (as JABBA-Select) and is constant
#' outside it. Biology at age comes from an FLicc \code{lhpar} (von Bertalanffy growth, weight and maturity at length,
#' M, steepness \code{s}); selectivity is defined at length with the FLicc functions and parameters (\code{selpar} on
#' the length scale, e.g. from \code{FLicc::selpars_flicc(fit)}) and converted to age by integrating over the
#' length-at-age distribution (coefficient of variation \code{CVL}). Sex-structured \code{lhpar} (dimension
#' \code{sex = c("F", "M")}): both sexes contribute to EB and SB with \code{sex_ratio}.
#'
#' A fleet whose catch is in weight and whose index is in numbers needs two curves from the same selectivity: give the
#' selectivity twice with different names and units, e.g. \code{selpar = list(LL = s, LL.N = s)},
#' \code{units = c("biomass", "numbers")}, then \code{jb3_model(ebsb = sx, ebsb.index = "LL.N", ebsb.catch = "LL")}.
#' @param lhpar FLicc life-history parameters (FLPar, named vector or \code{list(F = , M = )})
#' @param selpar FLicc selectivity parameters on the length scale: an FLPar (one fleet), FLPars or named list (one per
#'   fleet); the selectivity function is identified from the parameter names (logistic: SL50, SL95; dsnormal: mode,
#'   lsd, rsd; normal: mode, sd; binormal; bidsnormal)
#' @param sel.type "length" (default) or "age", per fleet; age-based selectivity uses the JABBA 5-parameter form
#'   (\code{\link{jb3_sel5}()}), which can also be used at length (parameters SL50, SL95, SLdesc, CVdesc, Smin)
#' @param units "biomass" (default) or "numbers", per fleet (recycled): units of the exploitable stock, i.e. of the
#'   index mapped with this curve
#' @param h steepness (default \code{lhpar$s}, else 0.7)
#' @param CVL coefficient of variation of length at age (0: selectivity at the mean length)
#' @param sex_ratio weights of the sexes (default equal)
#' @param amin,amax age range
#' @param Fgrid grid of F (default: from 0 until SB/SB0 < 0.02)
#' @param fit.range range of SB/SB0 over which the curve is fitted (default 0.1-0.9, as JABBA-Select)
#' @return object of class \code{jb3selex}: \code{ebsb} (named list of c(u1, u2, u3, P1, P2), ready for
#'   \code{jb3_model(ebsb = )}), \code{curves} (equilibrium ratio by fleet, column \code{EB_SB}, with \code{units}),
#'   \code{units}, \code{at.len} and \code{at.age} (selectivity and maturity), \code{leslie} (r and GT)
#' @seealso \code{\link{jb3_asem_prior}()} for the matching productivity prior
#' @examples
#' lh <- c(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.2, a = 0.0047, b = 3.15, s = 0.75)
#' sx <- jb3_selex(lh, list(Gillnet = c(mode = 32, lsd = 4, rsd = 15), Trawl = c(SL50 = 22, SL95 = 28)))
#' sx$ebsb
#' # catch in weight and an index in numbers from the same trawl selectivity
#' tr <- c(SL50 = 22, SL95 = 28)
#' sx2 <- jb3_selex(lh, list(Trawl = tr, Trawl.N = tr), units = c("biomass", "numbers"))
#' @export
jb3_selex <- function(lhpar, selpar, sel.type = "length", units = "biomass", h = NULL, CVL = 0.1, sex_ratio = NULL,
                      amin = 0, amax = NULL, Fgrid = NULL, fit.range = c(0.1, 0.9)){
  ss <- .selex_setup(lhpar, selpar, sel.type, CVL, sex_ratio, amin, amax)
  bio <- ss$bio; sela <- ss$sela; sr <- ss$sr; sexes <- ss$sexes; st <- ss$st; sf <- ss$sf; sel <- ss$sel
  h <- h %||% bio[[1]]$h
  if(length(units) > 1 && length(units) != length(sel)) stop("units: one value or one per fleet", call. = FALSE)
  un <- stats::setNames(rep_len(as.character(units), length(sel)), names(sel))
  if(any(!un %in% c("biomass", "numbers"))) stop("units must be \"biomass\" or \"numbers\"", call. = FALSE)
  curves <- do.call(rbind, lapply(names(sel), function(f){
    sb0 <- .eq_pr(0, bio, sela[[f]], sr)[["sb"]]
    Pof <- function(F){ p <- .eq_pr(F, bio, sela[[f]], sr); R <- max(0, (4*h*p[["sb"]] - (1 - h)*sb0)/((5*h - 1)*p[["sb"]]))
      c(P = p[["sb"]]*R/sb0, EB_SB = (if(un[[f]] == "numbers") p[["en"]] else p[["eb"]])/p[["sb"]]) }
    Fg <- Fgrid %||% { Fm <- 0.5; while(Pof(Fm)[["P"]] > 0.02 && Fm < 50) Fm <- 2*Fm; seq(0, Fm, length.out = 150) }
    d <- t(sapply(Fg, Pof)); d <- data.frame(fleet = f, units = un[[f]], F = Fg, d)
    d <- d[d$P > 0.005 & is.finite(d$EB_SB), ]; d[order(d$P), ]
  }))
  ebsb <- lapply(stats::setNames(names(sel), names(sel)), function(f){
    d <- curves[curves$fleet == f & curves$P >= fit.range[1] & curves$P <= fit.range[2], ]
    if(nrow(d) < 5) stop("fewer than 5 equilibrium points in fit.range for ", f, call. = FALSE)
    P1 <- min(d$P); P2 <- max(d$P)
    fit <- try(stats::nls(EB_SB ~ u1 + (u2 - u1)*(1 - exp(-u3*(P - P1)))/(1 - exp(-u3*(P2 - P1))), data = d,
                          start = list(u1 = d$EB_SB[1], u2 = d$EB_SB[nrow(d)], u3 = 3),
                          control = stats::nls.control(maxiter = 200, warnOnly = TRUE)), silent = TRUE)
    if(inherits(fit, "try-error")) stop("EB/SB curve fit failed for ", f, call. = FALSE)
    c(stats::coef(fit), P1 = P1, P2 = P2)
  })
  curves$fit <- unlist(lapply(names(sel), function(f) .rho_curve(curves$P[curves$fleet == f], ebsb[[f]])))
  L <- seq(0, max(sapply(bio, `[[`, "linf"))*1.1, length.out = 200)
  b1 <- bio[[1]]; ageL <- pmax(b1$t0 - log(pmax(1 - L/b1$linf, 1e-6))/b1$k, 0)      # inverse von Bertalanffy
  at.len <- do.call(rbind, c(lapply(names(sel), function(f) data.frame(curve = f, type = "selectivity", len = L,
    value = if(st[[f]] == "age") jb3_sel5(ageL, sel[[f]], "age") else .sel_len(L, sf[[f]], sel[[f]]))),
                             list(data.frame(curve = "maturity", type = "maturity", len = L,
                                             value = 1/(1 + exp(-log(19)*(L - bio[[1]]$L50)/(bio[[1]]$L95 - bio[[1]]$L50)))))))
  at.age <- do.call(rbind, c(lapply(names(sel), function(f) data.frame(curve = f, type = "selectivity", age = bio[[1]]$age, value = sela[[f]][[1]])),
                             list(data.frame(curve = "maturity", type = "maturity", age = bio[[1]]$age, value = bio[[1]]$mat))))
  structure(list(ebsb = ebsb, curves = curves, units = un, at.len = at.len, at.age = at.age, sel.fun = sf, sel.type = st, h = h,
                 sexes = sexes, fit.range = fit.range,
                 leslie = jb3_leslie(lhpar, h = h, sex = sexes[1], amin = amin, amax = amax)), class = "jb3selex")
}

# biology and selectivity at age for jb3_selex() and jb3_asem_prior()
.selex_setup <- function(lhpar, selpar, sel.type = "length", CVL = 0.1, sex_ratio = NULL, amin = 0, amax = NULL){
  sexes <- .sexes(lhpar); sel <- .selpar_list(selpar)
  st <- stats::setNames(rep_len(sel.type, length(sel)), names(sel))
  if(any(st == "age" & sapply(sel, .sel_fun_from) != "jabba5")) stop("age-based selectivity needs the JABBA 5-parameter form (SL50, SL95, SLdesc, CVdesc, Smin)", call. = FALSE)
  bio <- lapply(stats::setNames(sexes, sexes), function(s) .lh_age(.lh_one(lhpar, s), amin, amax))
  na <- min(sapply(bio, function(b) length(b$age))); bio <- lapply(bio, function(b){ for(nm in c("age", "len", "wt", "mat", "M")) b[[nm]] <- b[[nm]][seq_len(na)]; b })
  sr <- sex_ratio %||% rep(1/length(sexes), length(sexes))
  sf <- sapply(sel, .sel_fun_from)
  sela <- lapply(names(sel), function(f) lapply(bio, function(b)
    if(st[[f]] == "age") jb3_sel5(b$age, sel[[f]], "age") else .sel_age(b$len, sf[[f]], sel[[f]], CVL, b$linf)))
  names(sela) <- names(sel)
  list(bio = bio, sela = sela, sr = sr, sexes = sexes, st = st, sf = sf, sel = sel)
}

# equilibrium per recruit at fishing mortality F (selectivity S, a list over sexes): spawning biomass, exploitable
# biomass and numbers (start of year) and yield in weight (Baranov); the last age is a plus group
.eq_pr <- function(F, bio, S, sr, Mmult = 1){
  sb <- 0; eb <- 0; en <- 0; y <- 0
  for(j in seq_along(bio)){ b <- bio[[j]]; Sj <- S[[j]]; na <- length(b$age); M <- b$M*Mmult; Z <- M + F*Sj
    N <- numeric(na); N[1] <- 1
    if(na > 1) for(a in 2:na) N[a] <- N[a - 1]*exp(-Z[a - 1])
    N[na] <- N[na]/max(1e-12, 1 - exp(-Z[na]))
    sb <- sb + sr[j]*sum(N*b$wt*b$mat); eb <- eb + sr[j]*sum(N*b$wt*Sj); en <- en + sr[j]*sum(N*Sj)
    y <- y + sr[j]*sum(F*Sj/Z*(1 - exp(-Z))*N*b$wt) }
  c(sb = sb, eb = eb, en = en, y = y)
}

#' Productivity prior from the age-structured equilibrium, accounting for selectivity
#'
#' Monte Carlo over natural mortality and steepness through the age-structured equilibrium with the fishery's
#' selectivity (as the ASEM prior of JABBA-Select): for each draw, F_MSY is found by maximising the equilibrium yield
#' (Beverton-Holt recruitment), and the harvest rate at MSY and the production-function shape follow from the
#' equilibrium at F_MSY. The prior for r then matches the harvest rate at MSY of the model, which depends on the state
#' and on how catch is removed:
#' \itemize{
#'   \item standard JABBA3 (\code{state = "EB"}): K = EB0, B_MSY/K = EB_MSY/EB0 and r = m MSY/EB_MSY;
#'   \item Selex (\code{state = "SB"}), catch removed through the EB/SB curve (\code{jb3_model(ebsb.removal = TRUE)},
#'     the default; \code{removal = "EB"}): K = SB0, B_MSY/K = SB_MSY/SB0 and r = m MSY/EB_MSY;
#'   \item Selex with catch removed from SB directly (\code{ebsb.removal = FALSE}, as JABBA-Select v1;
#'     \code{removal = "SB"}): r = m MSY/SB_MSY.
#' }
#' Here m is the Pella-Tomlinson shape for B_MSY/K. Unlike \eqn{r_{Leslie}/2}, which ignores selectivity and the
#' stock-recruitment curve beyond its slope and can be wrong by a factor of two, this is the prior to use whenever the
#' biology and selectivity of an age-structured model are available.
#'
#' Use the medians as the model settings and the r prior conditional on the fixed shape:
#' \code{jb3_model(BmsyK = pr$BmsyK)} with \code{jb3_priors(r = pr$r)}. To estimate the shape, use \code{pr$m} as the m
#' prior (independent of r in JABBA3; JABBA-Select v1 used the correlated prior \code{pr$lmvn} on log(H_MSY, m)).
#' @param lhpar,selpar,sel.type,CVL,sex_ratio,amin,amax as in \code{\link{jb3_selex}()}; with several fleets the
#'   selectivity of fleet \code{catch} sets F_MSY
#' @param catch fleet (name or position) whose selectivity defines the fishery
#' @param state "SB" (Selex) or "EB" (standard JABBA3)
#' @param removal for \code{state = "SB"}: "EB" (catch removed through the EB/SB curve, default as in
#'   \code{jb3_model()}) or "SB" (removed from SB directly, \code{ebsb.removal = FALSE})
#' @param M,h natural mortality and steepness: single values (default from \code{lhpar}) or vectors of draws
#' @param M.cv,h.sd lognormal CV of M and normal SD of steepness for the draws when single values are given
#' @param nsim number of draws (1 without uncertainty)
#' @param seed random seed
#' @return object of class \code{jb3asem}: \code{r} (c(median, log.sd) for \code{jb3_priors(r = )}, conditional on the
#'   median shape), \code{BmsyK} (median, for \code{jb3_model(BmsyK = )}), \code{m} (c(median, log.sd, 1) for
#'   \code{jb3_priors(m = )}), \code{lmvn} (mean and covariance of log(H_MSY, m)), \code{draws} (per draw: M, h,
#'   F_MSY, MSY, SB0, SB_MSY, EB0, EB_MSY, B_MSY/K, m, H_MSY, r)
#' @seealso \code{\link{jb3_selex}()}, \code{\link{jb3_fmsy2r}()}, \code{\link{jb3_leslie}()}
#' @examples
#' lh <- c(linf = 147.5, k = 0.126, t0 = -1.89, a = 1.3718e-5, b = 3.0793, M = 0.3, s = 0.64,
#'         L50 = 79.9, L95 = 91.3, amax = 15)
#' sel <- c(SL50 = 57.2, SL95 = 72.4)
#' pr <- jb3_asem_prior(lh, sel, M.cv = 0.3, h.sd = 0.064, nsim = 200)
#' pr
#' # jb3_model(BmsyK = pr$BmsyK, ebsb = jb3_selex(lh, sel)); jb3_priors(r = pr$r)
#' @export
jb3_asem_prior <- function(lhpar, selpar, sel.type = "length", catch = 1, state = c("SB", "EB"), removal = c("EB", "SB"),
                           M = NULL, h = NULL, M.cv = 0, h.sd = 0, nsim = 1000, CVL = 0.1, sex_ratio = NULL, amin = 0,
                           amax = NULL, seed = 1){
  state <- match.arg(state); removal <- match.arg(removal)
  ss <- .selex_setup(lhpar, selpar, sel.type, CVL, sex_ratio, amin, amax)
  bio <- ss$bio; sr <- ss$sr
  S <- ss$sela[[catch]]; if(is.null(S)) stop("catch: fleet ", catch, " not found", call. = FALSE)
  M0 <- bio[[1]]$M[1]; h0 <- bio[[1]]$h
  if(length(M) > 1 || length(h) > 1) nsim <- max(length(M), length(h)) else if(M.cv == 0 && h.sd == 0) nsim <- 1
  set.seed(seed)
  Md <- if(!is.null(M) && length(M) > 1) rep_len(M, nsim) else (M %||% M0)*exp(stats::rnorm(nsim, -0.5*log(1 + M.cv^2), sqrt(log(1 + M.cv^2))))
  hd <- if(!is.null(h) && length(h) > 1) rep_len(h, nsim) else pmin(pmax(stats::rnorm(nsim, h %||% h0, h.sd), 0.21), 0.99)
  one <- function(Mi, hi){
    mm <- Mi/M0; pr0 <- .eq_pr(0, bio, S, sr, mm); sb0 <- pr0[["sb"]]
    eq <- function(F){ p <- .eq_pr(F, bio, S, sr, mm); R <- max(0, (4*hi*p[["sb"]] - (1 - hi)*sb0)/((5*hi - 1)*p[["sb"]]))
      c(R*p, R = R) }
    Fg <- exp(seq(log(0.01*Mi), log(10*Mi), length.out = 200)); Y <- sapply(Fg, function(F) eq(F)[["y"]])
    i <- which.max(Y); lo <- Fg[max(1, i - 1)]; hi_ <- Fg[min(length(Fg), i + 1)]
    Fm <- stats::optimize(function(F) -eq(F)[["y"]], c(lo, hi_))$minimum
    e <- eq(Fm); bk <- if(state == "EB") e[["eb"]]/pr0[["eb"]] else e[["sb"]]/sb0
    H <- e[["y"]]/(if(state == "SB" && removal == "SB") e[["sb"]] else e[["eb"]])
    c(M = Mi, h = hi, Fmsy = Fm, MSY = e[["y"]], SB0 = sb0, SBmsy = e[["sb"]], EB0 = pr0[["eb"]], EBmsy = e[["eb"]],
      BmsyK = bk, Hmsy = H)
  }
  d <- as.data.frame(t(mapply(one, Md, hd)))
  d$m <- sapply(d$BmsyK, .m_from_bmsyk); d$r <- d$Hmsy*d$m
  bk <- stats::median(d$BmsyK); mmed <- .m_from_bmsyk(bk)
  rc <- d$Hmsy*mmed                                                        # r conditional on the fixed median shape
  lsd <- function(x) if(length(x) > 1) stats::sd(log(x)) else 0
  lm <- log(cbind(Hmsy = d$Hmsy, m = d$m))
  structure(list(r = c(median = stats::median(rc), log.sd = lsd(rc)), BmsyK = bk,
                 m = c(median = stats::median(d$m), log.sd = lsd(d$m), on = 1),
                 lmvn = list(mean = colMeans(lm), cov = if(nrow(lm) > 1) stats::cov(lm) else matrix(0, 2, 2)),
                 r.marginal = c(median = stats::median(d$r), log.sd = lsd(d$r)), state = state, removal = removal,
                 draws = d, leslie = jb3_leslie(lhpar, h = h0, sex = ss$sexes[1], amin = amin, amax = amax)), class = "jb3asem")
}

#' @export
print.jb3asem <- function(x, ...){
  d <- x$draws
  cat(sprintf("JABBA3 ASEM productivity prior (%d draws; state %s%s)\n", nrow(d), x$state,
              if(x$state == "SB") paste0(", catch removed from ", x$removal) else ""))
  cat(sprintf("  F_MSY %.3f | H_MSY %.3f | B_MSY/K %.3f (m = %.3f)\n", stats::median(d$Fmsy), stats::median(d$Hmsy), x$BmsyK,
              .m_from_bmsyk(x$BmsyK)))
  cat(sprintf("  r prior (shape fixed at the median): median %.3f, log.sd %.3f\n", x$r[1], x$r[2]))
  cat(sprintf("  for comparison: r_Leslie/2 converted to this shape = %.3f\n", jb3_fmsy2r(x$leslie[["r"]]/2, x$BmsyK)))
  cat("  use: jb3_model(BmsyK = $BmsyK), jb3_priors(r = $r)\n"); invisible(x)
}

#' @export
print.jb3selex <- function(x, ...){
  cat("JABBA3 Selex inputs:", length(x$ebsb), "selectivity curve(s);", if(length(x$sexes) > 1) "sexes" else "sex",
      paste(x$sexes, collapse = ", "), "; steepness", x$h, "\n")
  cat(sprintf("  Leslie r = %.3f, generation time = %.1f years\n", x$leslie[["r"]], x$leslie[["gt"]]))
  u <- x$units %||% stats::setNames(rep("biomass", length(x$ebsb)), names(x$ebsb))
  print(cbind(as.data.frame(signif(do.call(rbind, x$ebsb), 4)), units = u))
  cat("  pass the object to jb3_model(ebsb = ); with several curves set ebsb.index and ebsb.catch\n"); invisible(x)
}

#' Selectivity and EB/SB curves
#'
#' Left: selectivity by fleet with the maturity ogive (dashed), at length (default) or at age. Right: the equilibrium
#' ratio of exploitable to spawning biomass against spawning-biomass depletion (points) with the fitted Selex curve
#' (lines), as used by \code{jb3_model(ebsb = )}.
#' @param x output of \code{\link{jb3_selex}()}
#' @param at "length" or "age"
#' @return a patchwork (or list of ggplots)
#' @export
jbplot_selex <- function(x, at = c("length", "age")){
  at <- match.arg(at)
  d <- if(at == "length") transform(x$at.len, xx = len) else transform(x$at.age, xx = age)
  fl <- names(x$ebsb); cols <- stats::setNames(jb3cols(length(fl)), fl)
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_line(data = d[d$type == "maturity", ], ggplot2::aes(.data$xx, .data$value), linetype = 2, linewidth = 0.7) +
    ggplot2::geom_line(data = d[d$type == "selectivity", ], ggplot2::aes(.data$xx, .data$value, colour = .data$curve), linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = if(at == "length") "Length" else "Age", y = "Selectivity, maturity (dashed)") + theme_jabba3()
  cv <- x$curves; x$units <- x$units %||% stats::setNames(rep("biomass", length(x$ebsb)), names(x$ebsb))
  x$fit.range <- x$fit.range %||% c(0, 1); if(is.null(cv$units)) cv$units <- "biomass"
  rel <- any(x$units == "numbers") && any(x$units == "biomass")
  if(rel) for(f in unique(cv$fleet)){ w <- cv$fleet == f; s <- .rho_curve(1, x$ebsb[[f]]); cv$EB_SB[w] <- cv$EB_SB[w]/s; cv$fit[w] <- cv$fit[w]/s }
  ylab <- if(rel) "Exploitable / SB (relative to SB0)" else if(all(x$units == "numbers")) "EN/SB" else "EB/SB"
  p2 <- ggplot2::ggplot(cv, ggplot2::aes(.data$P, colour = .data$fleet, linetype = .data$units)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$EB_SB), size = 0.8, alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fit), linewidth = 0.9) + ggplot2::geom_hline(yintercept = 1, linetype = 3) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(biomass = 1, numbers = 2), name = NULL, guide = "none") +
    ggplot2::annotate("rect", xmin = x$fit.range[1], xmax = x$fit.range[2], ymin = -Inf, ymax = Inf, alpha = 0.06) +
    ggplot2::labs(x = expression(SB/SB[0]), y = ylab) + theme_jabba3()
  if(!requireNamespace("patchwork", quietly = TRUE)) return(list(p1, p2))
  patchwork::wrap_plots(p1, p2, ncol = 2) + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")
}
