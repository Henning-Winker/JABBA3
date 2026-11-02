# Seasonal time steps (nseason): data parsing and step-level setup
#
# The model runs on nseason discrete steps per year (dt = 1/nseason, an Euler scheme). The state of step (y, s) is the
# depletion at the start of season s of year y; catch of a step is removed during that step. With nseason = 1 everything
# reduces to the annual JABBA model.

.season_col <- function(x){
  if(!is.data.frame(x)) return(NULL)
  nm <- tolower(names(x))
  if(length(nm) >= 3 && nm[2] == "season") return("season")
  if(nm[1] == "time") return("time")
  NULL
}

# long table (year, frac, season, values...) from a data.frame with a season or time column
.season_long <- function(x, type, nsd){
  x <- as.data.frame(x)
  if(type == "season"){
    y <- as.numeric(x[[1]]); s <- suppressWarnings(as.numeric(x[[2]])); v <- x[, -(1:2), drop = FALSE]
    frac <- ifelse(is.na(s), NA, (s - 1)/nsd)
  } else {
    tm <- as.numeric(x[[1]]); y <- floor(tm + 1e-8); frac <- tm - y; v <- x[, -1, drop = FALSE]
    s <- NA
  }
  list(year = y, frac = frac, season = s, v = as.matrix(v))
}

.season_parse <- function(catch, indices, se, nseason, index.time){
  ct <- .season_col(catch); it <- .season_col(indices)
  if(is.null(ct) && is.null(it)) return(NULL)
  nsd <- nseason
  if(is.null(nsd)){
    if(identical(ct, "time") || identical(it, "time")) stop("nseason is needed with a time column", call. = FALSE)
    nsd <- max(c(if(!is.null(ct)) suppressWarnings(as.numeric(catch[[2]])), if(!is.null(it)) suppressWarnings(as.numeric(indices[[2]]))), na.rm = TRUE)
  }
  nsd <- as.integer(nsd); if(nsd < 1) stop("nseason must be >= 1", call. = FALSE)
  out <- list(nsd = nsd, catch.annual = catch, I.annual = indices, SE.annual = se, catch.season = NULL, iobs = NULL)
  if(!is.null(ct)){
    cl <- .season_long(catch, ct, nsd)
    cs <- if(ct == "time") floor(cl$frac*nsd + 1e-6) + 1 else cl$season
    val <- rowSums(cl$v, na.rm = TRUE)
    yrs <- sort(unique(cl$year))
    CS <- matrix(NA_real_, length(yrs), nsd, dimnames = list(yrs, paste0("s", seq_len(nsd))))
    ann <- stats::setNames(numeric(length(yrs)), yrs)
    for(k in seq_along(val)){
      j <- as.character(cl$year[k])
      if(is.na(cs[k])) ann[j] <- ann[j] + val[k] else {
        if(cs[k] > nsd) stop("catch season > nseason", call. = FALSE)
        CS[j, cs[k]] <- sum(CS[j, cs[k]], val[k], na.rm = TRUE)
      }
    }
    part <- rowSums(!is.na(CS)) > 0 & rowSums(!is.na(CS)) < nsd
    if(any(part)) stop("seasonal catch incomplete in years ", paste(yrs[part], collapse = ", "),
                       " (give all seasons or the annual total with season NA)", call. = FALSE)
    both <- rowSums(!is.na(CS)) == nsd & ann > 0
    if(any(both)) stop("catch given both annually and by season in years ", paste(yrs[both], collapse = ", "), call. = FALSE)
    tot <- ifelse(rowSums(!is.na(CS)) == nsd, rowSums(CS), ann)
    out$catch.annual <- data.frame(year = yrs, catch = as.numeric(tot))
    out$catch.season <- CS
  }
  if(!is.null(it)){
    il <- .season_long(indices, it, nsd)
    sl <- if(is.null(se)) NULL else .season_long(se, it, nsd)
    if(!is.null(sl) && (ncol(sl$v) != ncol(il$v) || length(sl$year) != length(il$year)))
      stop("se must have the same rows and columns as the seasonal indices", call. = FALSE)
    nm <- colnames(il$v)
    frac <- ifelse(is.na(il$frac), 0, il$frac)
    long <- do.call(rbind, lapply(seq_along(nm), function(i) data.frame(year = il$year, frac = frac, index = nm[i],
      obs = as.numeric(il$v[, i]), se = if(is.null(sl)) NA_real_ else as.numeric(sl$v[, i]))))
    long <- long[is.finite(long$obs) & long$obs > 0, ]
    long$se[!is.finite(long$se) | long$se <= 0] <- NA
    out$iobs <- long
    # annual (geometric mean) indices for the annual model
    yrs <- sort(unique(il$year))
    I <- data.frame(year = yrs); S <- data.frame(year = yrs)
    for(i in nm){
      w <- long[long$index == i, ]
      I[[i]] <- vapply(yrs, function(y){ o <- w$obs[w$year == y]; if(length(o)) exp(mean(log(o))) else NA_real_ }, 0)
      S[[i]] <- vapply(yrs, function(y){ o <- w$se[w$year == y]; if(length(o) && any(!is.na(o))) sqrt(mean(o^2, na.rm = TRUE)/sum(!is.na(o))) else NA_real_ }, 0)
    }
    out$I.annual <- I; out$SE.annual <- if(all(is.na(as.matrix(S[, -1])))) NULL else S
  }
  out
}

# attach seasonal fields to a jb3data object
.season_attach <- function(d, sea, index.time){
  d$nseason <- if(is.null(sea)) 1L else sea$nsd
  d$catch.season <- NULL
  if(!is.null(sea$catch.season)){
    CS <- matrix(NA_real_, length(d$years), sea$nsd, dimnames = list(d$years, colnames(sea$catch.season)))
    k <- match(rownames(sea$catch.season), d$years); CS[k[!is.na(k)], ] <- sea$catch.season[!is.na(k), ]
    if(any(!is.na(CS))) d$catch.season <- CS
  }
  nI <- ncol(d$I)
  if(!is.null(sea$iobs)){
    io <- sea$iobs[sea$iobs$year %in% d$years & sea$iobs$index %in% colnames(d$I), ]
  } else if(nI > 0){
    tm <- rep(index.time, length.out = nI)
    if(any(tm < 0 | tm >= 1)) stop("index.time must be in [0, 1)", call. = FALSE)
    io <- do.call(rbind, lapply(seq_len(nI), function(i){ o <- which(!is.na(d$I[, i]))
      data.frame(year = d$years[o], frac = tm[i], index = colnames(d$I)[i], obs = d$I[o, i], se = d$SE[o, i]) }))
  } else io <- data.frame(year = numeric(0), frac = numeric(0), index = character(0), obs = numeric(0), se = numeric(0))
  rownames(io) <- NULL
  d$iobs <- io
  d$seasonal.index <- !is.null(sea$iobs)
  d
}

# step-level data for the engine (called from jb3_setup when model$nseason > 1 or seasonal inputs are present)
# returns: N (steps), Ny, years, stepyear, stepseason, I/SE matrices (steps x indices), catch shares, allocation flags
.season_steps <- function(data, model, years){
  ns <- as.integer(model$nseason %||% 1L); Ny <- length(years); N <- Ny*ns
  stepyear <- rep(seq_len(Ny), each = ns); stepseason <- rep(seq_len(ns), Ny)
  nm <- colnames(data$I); nI <- length(nm)
  I <- SE <- matrix(NA_real_, N, max(nI, 0), dimnames = list(paste(rep(years, each = ns), seq_len(ns), sep = "."), nm))
  io <- data$iobs
  if(nI > 0 && nrow(io)){
    if(ns == 1 && isTRUE(data$seasonal.index)){
      I <- data$I; SE <- data$SE           # annual model on seasonal indices: annual geometric means
      message("nseason = 1 with seasonal indices: annual geometric means used")
    } else {
      yi <- match(io$year, years); sh <- round(io$frac*ns)          # nearest season start
      st <- (yi - 1)*ns + sh + 1
      ok <- !is.na(st) & st <= N
      if(any(!ok)) warning(sum(!ok), " index observations after the last step dropped", call. = FALSE)
      io <- io[ok, ]; st <- st[ok]
      for(k in seq_len(nrow(io))){
        i <- match(io$index[k], nm)
        if(!is.na(I[st[k], i])){
          warning("several observations of index ", nm[i], " in step ", rownames(I)[st[k]], ": geometric mean", call. = FALSE)
          I[st[k], i] <- sqrt(I[st[k], i]*io$obs[k])
        } else { I[st[k], i] <- io$obs[k]; SE[st[k], i] <- io$se[k] }
      }
    }
  }
  # catch: annual totals and within-year shares
  alloc <- model$catch.alloc %||% "F"
  share <- rep(1/ns, N); allocF <- rep(0L, Ny)
  if(is.numeric(alloc)){
    if(length(alloc) != ns || any(alloc < 0)) stop("catch.alloc shares: one non-negative value per season", call. = FALSE)
    share <- rep(alloc/sum(alloc), Ny)
  } else if(alloc %in% c("F", "Fstate") && ns > 1) allocF[] <- if(alloc == "F") 1L else 2L
  CS <- data$catch.season
  if(!is.null(CS) && ns > 1){
    if(ncol(CS) != ns) stop("seasonal catch has ", ncol(CS), " seasons but jb3_model(nseason = ", ns, ")", call. = FALSE)
    for(y in which(rowSums(!is.na(CS)) == ns)){
      tot <- sum(CS[y, ]); share[(y - 1)*ns + seq_len(ns)] <- if(tot > 0) CS[y, ]/tot else rep(1/ns, ns); allocF[y] <- 0L
    }
  }
  w <- model$season.w %||% rep(1, ns)
  if(length(w) != ns || any(w < 0)) stop("season.w: one non-negative weight per season", call. = FALSE)
  w <- w/mean(w)
  list(ns = ns, N = N, Ny = Ny, stepyear = stepyear, stepseason = stepseason, I = I, SE = SE, share = share,
       allocF = allocF, w = w, yr1 = (seq_len(Ny) - 1)*ns + 1)
}

# variance factor: v_step = sd_annual^2 / f, f = ns + 2 sum_{k=1}^{ns-1} (ns - k) rho_s^k, so that the sum of the ns step
# deviations of a year has variance sd_annual^2 (rho_s = rho^(1/ns) for AR(1); rho_s = 0 for iid)
.season_vfac <- function(ns, rhos){
  f <- ns + 0*rhos
  if(ns > 1) for(k in 1:(ns - 1)) f <- f + 2*(ns - k)*rhos^k
  f
}
