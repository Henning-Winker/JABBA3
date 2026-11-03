#' Assessment data for JABBA3
#'
#' Collects catch and abundance indices into the internal data object used by \code{\link{jabba3}()}.
#' Native inputs follow FLR: \code{catch} as an \code{FLQuant} (or \code{FLQuants} by fleet, summed) and
#' \code{indices} as \code{FLIndices} of \code{FLIndexBiomass}, with the log-scale SE of each observation in
#' \code{index.var} (NA or empty = not given). Plain data frames in the JABBA layout are accepted as well
#' (see also \code{\link{jb3_data_jabba}()}).
#'
#' @param catch \code{FLQuant}, \code{FLQuants}, or data.frame(year, catch[, more fleets])
#' @param indices \code{FLIndices}, \code{FLIndexBiomass}, or data.frame(year, index1, index2, ...)
#' @param se optional data.frame(year, se1, se2, ...) of log-scale standard errors (data.frame input only)
#' @param catch.se optional year-specific log-scale SE of the catch: a vector along the catch years or data.frame(year, se)
#'   (e.g. from \code{\link{jb3_catchfill}()}); used by the catch error (\code{jb3_model(catch.error = TRUE)}) in place
#'   of the constant \code{sdc}; NA = \code{sdc}
#' @param stock stock name
#' @param auxiliary optional auxiliary data (JABBA \code{auxiliary}): data.frame(year, series1, ...), FLQuants (e.g. from
#'   \code{jbdat2flr()}, SEs in attribute "se"), or a named list of
#'   series, each a data.frame(year, value[, se]) or a list(data = data.frame(year, value), se = , type = ). Types (see
#'   \code{auxiliary.type}): \code{"effort"} (fishing effort, modelled as F / (q_E x creep)), \code{"f"} (fishing mortality,
#'   harvest rate), \code{"ffmsy"}, \code{"bbmsy"}, \code{"bk"} (status estimates, e.g. from another assessment or a
#'   length-based method) and \code{"z"} (total mortality, as JABBA). See \code{\link{jb3_model}()} for the
#'   catchability, additional SD and effort creep of the auxiliary series.
#' @param auxiliary.se optional log-scale SEs of the auxiliary data, same layout as \code{auxiliary} (data.frame input)
#' @param auxiliary.type type of each auxiliary series (recycled), one of "effort", "f", "ffmsy", "bbmsy", "bk", "z"
#' @param catch.gaps "stop" (default): the model years are the catch years and catch must be given in all of them;
#'   "estimate": the model years span the catch and the indices, and catch in years without catch (gaps, years before the
#'   first or after the last catch year) is estimated with a random walk on log harvest rate (see \code{\link{jb3_model}()},
#'   \code{sdf}), as in SPiCT
#' @details \strong{Alignment.} The series need not cover the same years. The model years are the catch years (leading and
#'   trailing years without catch are trimmed; gaps inside the catch series are an error, see \code{\link{jb3_catchfill}()}).
#'   Each index keeps its own years and is matched by year; index observations outside the catch years are not used (message).
#'   Rows may be in any order; the year column is the column named \code{year} (or the first column). Standard errors are matched
#'   to the indices by year and by column name (or by position when the names differ but the number of columns agrees); SEs may
#'   be given for some years or some indices only. Catch by fleet is summed with missing values counted as no catch, unless all
#'   fleets are missing in a year. \code{FLIndices} may hold \code{FLIndexBiomass} objects of different lengths. The years
#'   used are listed by \code{\link{jb3_align}()}.
#'
#' @param nseason number of seasons of seasonal inputs (needed with a \code{time} column; with a \code{season} column the
#'   default is the largest season). Whether the model runs on seasons is set by \code{jb3_model(nseason = )}.
#' @param index.time timing of annual indices within the year (fraction of the year, recycled over indices; 0 = start of
#'   the year, as JABBA). Used when the model runs on seasons: the index is matched to the season starting nearest to it.
#' @details \strong{Seasonal data.} Catch and indices may carry a \code{season} column (1, ..., nseason) or a
#'   \code{time} column (decimal year as in SPiCT, e.g. 2001.25 for the start of the second quarter) after \code{year}
#'   (data frames: \code{year, season, catch} or \code{time, catch}). Catch in a season covers that season; years without
#'   seasonal values may be given annually (rows with season NA). Seasonal index values are observations of the biomass
#'   at the start of their season (\code{time}: matched to the nearest season start). The annual catch and the annual
#'   (geometric mean) indices are kept in \code{catch}, \code{I} for the annual model; the seasonal values in
#'   \code{catch.season} (years x seasons) and \code{iobs} (long table: year, frac, index, obs, se).
#'
#'   \strong{Catch-only.} \code{indices = NULL} gives a catch-only data set: the model is then informed by the catch, the priors
#'   (including status priors, \code{jb3_priors(status = )}) and any auxiliary data; \code{jabba3()} runs it with NUTS by
#'   default.
#' @return object of class \code{jb3data}: \code{years}, \code{catch}, index matrix \code{I}, SE matrix \code{SE}
#'   (years x indices, zero columns for catch-only), auxiliary matrices \code{A}, \code{A.SE} and types \code{A.type}
#'   (NULL without auxiliary data), \code{stock}
#' @examples
#' data(iccat)
#' d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
#' d
#' # catch-only with hypothetical effort (simulated, not real data)
#' d0 <- jb3_data(iccat$bet$catch, NULL, auxiliary = iccat$bet$effort_hyp[, c("year", "effort")],
#'                auxiliary.type = "effort", stock = "bet")
#' @export
jb3_data <- function(catch, indices = NULL, se = NULL, stock = "stock", catch.se = NULL, auxiliary = NULL, auxiliary.se = NULL,
                     auxiliary.type = "effort", catch.gaps = c("stop", "estimate"), nseason = NULL, index.time = 0){
  catch.gaps <- match.arg(catch.gaps)
  if(inherits(catch, "jb3catchfill")){ if(is.null(catch.se)) catch.se <- data.frame(year = catch$year, se = catch$se.log)
    catch <- data.frame(year = catch$year, catch = catch$catch) }
  if(inherits(indices, "jb3jointindex") && is.null(se)){ se <- data.frame(year = indices$year, joint = indices$log.se)
    indices <- data.frame(year = indices$year, joint = indices$index) }
  # seasonal inputs: a "season" or "time" column in catch and/or indices (see Details)
  sea <- .season_parse(catch, indices, se, nseason, index.time)
  if(!is.null(sea)){ catch <- sea$catch.annual; indices <- sea$I.annual; se <- sea$SE.annual }
  if(!is.null(sea) && catch.gaps == "estimate") stop("catch.gaps = \"estimate\" is not yet available with seasonal inputs", call. = FALSE)
  cy <- .catch_df(catch)
  ok <- is.finite(cy$catch)
  if(!any(ok)) stop("no finite catch values", call. = FALSE)
  years <- seq(min(cy$year[ok]), max(cy$year[ok]))                  # leading/trailing NA years (padding) are trimmed
  ix <- if(is.null(indices)) NULL else .index_df(indices, se)
  if(catch.gaps == "estimate"){
    # model years span the catch and the indices; catch in years without catch is estimated (F random walk, see jb3_model())
    yi <- if(is.null(ix)) numeric(0) else ix$I$year[rowSums(is.finite(as.matrix(ix$I[, -1, drop = FALSE])) & as.matrix(ix$I[, -1, drop = FALSE]) > 0) > 0]
    years <- seq(min(c(years, yi)), max(c(years, yi)))
  }
  C <- cy$catch[match(years, cy$year)]
  if(any(!is.finite(C))){
    if(catch.gaps == "stop") stop("catch missing in years ", paste(years[!is.finite(C)], collapse = ", "),
      " (inside the catch series ", min(years), "-", max(years), "); fill gaps first, e.g. with jb3_catchfill(), or use ",
      "jb3_data(catch.gaps = \"estimate\") to estimate them with an F random walk", call. = FALSE)
    message("catch missing in ", sum(!is.finite(C)), " of ", length(years), " model years (", .yr_ranges(years[!is.finite(C)]),
            "): estimated with an F random walk")
  }
  align <- data.frame(series = "catch", type = "catch", first = min(cy$year), last = max(cy$year), n = sum(ok),
                      used_first = min(years), used_last = max(years), n_used = length(years))
  cse <- NULL
  if(!is.null(catch.se)){
    cse <- if(is.data.frame(catch.se)) as.numeric(catch.se[[2]])[match(years, as.numeric(catch.se[[1]]))] else rep(as.numeric(catch.se), length.out = length(years))
  }
  aux <- .aux_df(auxiliary, auxiliary.se, auxiliary.type, years)
  if(is.null(indices)){                                   # catch-only
    I <- SE <- matrix(NA_real_, length(years), 0, dimnames = list(years, NULL))
    return(.season_attach(structure(list(stock = stock, years = years, catch = C, I = I, SE = SE, catch.fleets = cy$fleets, catch.se = cse,
                          A = aux$A, A.SE = aux$SE, A.type = aux$type, align = align, catch.obs = is.finite(C)), class = "jb3data"), sea, index.time))
  }
  # alignment by year: each index keeps its own years within the model years
  for(k in names(ix$I)[-1]){
    v <- ix$I[[k]]; o <- is.finite(v) & v > 0; yo <- ix$I$year[o]; u <- yo >= min(years) & yo <= max(years)
    align <- rbind(align, data.frame(series = k, type = "index", first = if(any(o)) min(yo) else NA, last = if(any(o)) max(yo) else NA,
      n = sum(o), used_first = if(any(u)) min(yo[u]) else NA, used_last = if(any(u)) max(yo[u]) else NA, n_used = sum(u)))
  }
  drop <- align[align$type == "index" & align$n_used < align$n, ]
  if(nrow(drop)) message("index observations outside the model years (", min(years), "-", max(years), ") are not used: ",
    paste(sprintf("%s %d of %d", drop$series, drop$n - drop$n_used, drop$n), collapse = "; "),
    if(catch.gaps == "stop") " (catch needed in every model year; or jb3_data(catch.gaps = \"estimate\"))" else "")
  I <- as.matrix(ix$I[match(years, ix$I$year), -1, drop = FALSE])
  SE <- as.matrix(ix$SE[match(years, ix$SE$year), -1, drop = FALSE])
  I[!is.finite(I) | I <= 0] <- NA
  SE[!is.finite(SE) | SE <= 0] <- NA
  dimnames(I) <- dimnames(SE) <- list(years, colnames(ix$I)[-1])
  keep <- colSums(!is.na(I)) > 0
  if(!any(keep)) stop("no index observations within the catch years", call. = FALSE)
  if(any(!keep)) warning("indices without observations within the catch years dropped: ", paste(colnames(I)[!keep], collapse = ", "), call. = FALSE)
  rownames(align) <- NULL
  .season_attach(structure(list(stock = stock, years = years, catch = C, I = I[, keep, drop = FALSE], SE = SE[, keep, drop = FALSE],
                 catch.fleets = cy$fleets, catch.se = cse, A = aux$A, A.SE = aux$SE, A.type = aux$type, align = align,
                 catch.obs = is.finite(C)), class = "jb3data"), sea, index.time)
}

.aux_types <- c("effort", "f", "ffmsy", "bbmsy", "bk", "z")

# auxiliary data to matrices along the catch years
.aux_df <- function(aux, se, type, years){
  if(is.null(aux)) return(list(A = NULL, SE = NULL, type = NULL))
  N <- length(years)
  if(inherits(aux, "FLQuants") || inherits(aux, "FLQuant")){ x <- .flqs2df(aux, years); aux <- x$auxiliary; se <- x$auxiliary.se %||% se }
  if(is.data.frame(aux) || is.matrix(aux)){
    aux <- as.data.frame(aux); if(ncol(aux) < 2) stop("auxiliary needs columns year and at least one series", call. = FALSE)
    nm <- names(aux)[-1]; y <- as.numeric(aux[[1]])
    A <- as.matrix(aux[match(years, y), -1, drop = FALSE])
    S <- if(is.null(se)) matrix(NA_real_, N, ncol(A)) else { se <- as.data.frame(se)
      if(ncol(se) != ncol(aux)) stop("auxiliary.se must have the same columns as auxiliary", call. = FALSE)
      as.matrix(se[match(years, as.numeric(se[[1]])), -1, drop = FALSE]) }
    tp <- rep(type, length.out = ncol(A))
  } else if(is.list(aux)){
    nm <- names(aux); if(is.null(nm) || any(nm == "")) nm <- paste0("aux", seq_along(aux))
    tp0 <- rep(type, length.out = length(aux))
    A <- S <- matrix(NA_real_, N, length(aux)); tp <- character(length(aux))
    for(k in seq_along(aux)){
      a <- aux[[k]]; s <- NULL; tk <- tp0[k]
      if(!is.data.frame(a) && is.list(a)){ s <- a$se; tk <- a$type %||% tk; a <- a$data }
      a <- as.data.frame(a); y <- as.numeric(a[[1]])
      A[, k] <- as.numeric(a[[2]])[match(years, y)]
      if(is.null(s) && ncol(a) >= 3) s <- as.numeric(a[[3]])
      if(!is.null(s)){ s <- if(is.data.frame(s)) as.numeric(s[[2]])[match(y, as.numeric(s[[1]]))] else rep(as.numeric(s), length.out = length(y))
        S[, k] <- s[match(years, y)] }
      tp[k] <- tk
    }
  } else stop("auxiliary must be a data.frame or a list of series", call. = FALSE)
  tp <- tolower(tp)
  if(!all(tp %in% .aux_types)) stop("auxiliary.type must be one of ", paste(.aux_types, collapse = ", "), call. = FALSE)
  A[!is.finite(A) | A <= 0] <- NA; S[!is.finite(S) | S <= 0] <- NA
  dimnames(A) <- dimnames(S) <- list(years, nm)
  keep <- colSums(!is.na(A)) > 0
  if(!any(keep)){ warning("no auxiliary observations within the catch years", call. = FALSE); return(list(A = NULL, SE = NULL, type = NULL)) }
  list(A = A[, keep, drop = FALSE], SE = S[, keep, drop = FALSE], type = stats::setNames(tp[keep], nm[keep]))
}

#' Convert classic JABBA inputs
#'
#' The JABBA data frames (\code{catch}, \code{cpue}, \code{se} as used by \code{build_jabba()}) to a
#' \code{jb3data} object.
#' @param catch data.frame(year, catch...)
#' @param cpue data.frame(year, index...)
#' @param se optional data.frame(year, se...)
#' @param assessment stock name
#' @export
jb3_data_jabba <- function(catch, cpue, se = NULL, assessment = "stock") jb3_data(catch, cpue, se, stock = assessment)

.catch_df <- function(catch){
  if(inherits(catch, "FLQuants") || (inherits(catch, "FLQuant"))){
    fl <- if(inherits(catch, "FLQuants")) as.list(catch) else list(catch)
    if(any(vapply(fl, function(x) dim(x)[6] > 1, TRUE))) warning("catch has iterations: the first is used", call. = FALSE)
    yrs <- sort(unique(unlist(lapply(fl, function(x) as.numeric(dimnames(x)$year)))))
    M <- vapply(fl, function(x){ y <- as.numeric(dimnames(x)$year); as.numeric(x[1, , 1, 1, 1, 1])[match(yrs, y)] }, numeric(length(yrs)))
    M <- matrix(M, length(yrs))
    return(list(year = yrs, catch = .sum_na(M), fleets = length(fl)))
  }
  catch <- as.data.frame(catch)
  yc <- .year_col(catch)
  if(ncol(catch) < 2) stop("catch needs columns year and catch", call. = FALSE)
  y <- as.numeric(catch[[yc]]); v <- as.matrix(catch[, -yc, drop = FALSE]); storage.mode(v) <- "double"
  if(any(duplicated(y))) stop("catch: duplicated years ", paste(unique(y[duplicated(y)]), collapse = ", "), call. = FALSE)
  o <- order(y)
  list(year = y[o], catch = .sum_na(v[o, , drop = FALSE]), fleets = ncol(v))
}

# sum over fleets: NA in a fleet counts as no catch, unless all fleets are NA in that year
.sum_na <- function(M){ out <- rowSums(M, na.rm = TRUE); out[rowSums(!is.na(M)) == 0] <- NA; out }

# compact year ranges, e.g. 1950-1959, 1990
.yr_ranges <- function(y){ if(!length(y)) return(""); g <- cumsum(c(1, diff(y) != 1))
  paste(vapply(split(y, g), function(v) if(length(v) > 1) paste0(min(v), "-", max(v)) else as.character(v), ""), collapse = ", ") }

# year column: a column named "year" (any case; also "yr") or else the first column
.year_col <- function(df){ k <- which(tolower(names(df)) %in% c("year", "yr", "years")); if(length(k)) k[1] else 1L }

.index_df <- function(indices, se){
  if(inherits(indices, "FLIndex") || inherits(indices, "FLIndexBiomass")) indices <- list(index = indices)
  if(inherits(indices, "FLIndices") || (is.list(indices) && !is.data.frame(indices) && length(indices) &&
                                        all(vapply(indices, function(x) methods::is(x, "FLIndex") || methods::is(x, "FLIndexBiomass"), TRUE)))){
    nm <- names(indices); if(is.null(nm) || any(nm == "")) nm <- paste0("index", seq_along(indices))
    yrs <- sort(unique(unlist(lapply(indices, function(x) as.numeric(dimnames(FLCore::index(x))$year)))))
    I <- data.frame(year = yrs); S <- data.frame(year = yrs)
    for(k in seq_along(indices)){
      x <- indices[[k]]; iq <- FLCore::index(x); y <- as.numeric(dimnames(iq)$year)
      v <- tryCatch(FLCore::index.var(x), error = function(e) NULL)
      I[[nm[k]]] <- as.numeric(iq[1, , 1, 1, 1, 1])[match(yrs, y)]
      S[[nm[k]]] <- if(is.null(v) || all(!is.finite(v))) NA else as.numeric(v[1, , 1, 1, 1, 1])[match(yrs, y)]
    }
    return(list(I = I, SE = S))
  }
  I <- as.data.frame(indices)
  if(ncol(I) < 2) stop("indices need columns year and at least one index", call. = FALSE)
  yc <- .year_col(I); I <- data.frame(year = as.numeric(I[[yc]]), I[, -yc, drop = FALSE], check.names = FALSE)
  if(any(duplicated(I$year))) stop("indices: duplicated years ", paste(unique(I$year[duplicated(I$year)]), collapse = ", "),
                                   " (one row per year; use one column per index)", call. = FALSE)
  I <- I[order(I$year), , drop = FALSE]
  S <- I; S[, -1] <- NA
  if(!is.null(se)){
    se <- as.data.frame(se); ys <- .year_col(se)
    se <- data.frame(year = as.numeric(se[[ys]]), se[, -ys, drop = FALSE], check.names = FALSE)
    if(any(duplicated(se$year))) stop("se: duplicated years", call. = FALSE)
    inm <- names(I)[-1]; snm <- names(se)[-1]
    # columns matched by name when the se names are index names; otherwise by position (same number of columns)
    map <- if(all(snm %in% inm)) stats::setNames(snm, snm) else
      if(!any(snm %in% inm) && length(snm) == length(inm)) stats::setNames(snm, inm) else      # JABBA convention: by position
      stop("se: columns must be named as the indices (a subset is allowed), or one column per index in the same order", call. = FALSE)
    for(k in names(map)) S[[k]] <- as.numeric(se[[map[[k]]]])[match(I$year, se$year)]
  }
  list(I = I, SE = S)
}

#' Alignment of the input series
#'
#' Year ranges of the catch and index inputs as given, and the part used after aligning all series to the catch years.
#' @param x \code{jb3_data()} object
#' @return data.frame(series, type, first, last, n, used_first, used_last, n_used)
#' @export
jb3_align <- function(x) x$align

#' @export
print.jb3data <- function(x, ...){
  cat("JABBA3 data:", x$stock, "\n  years ", min(x$years), "-", max(x$years), " (", length(x$years), "), catch range ",
      signif(min(x$catch, na.rm = TRUE), 3), "-", signif(max(x$catch, na.rm = TRUE), 3), "\n", sep = "")
  if(any(is.na(x$catch))) cat("  catch missing in ", sum(is.na(x$catch)), " years (", .yr_ranges(x$years[is.na(x$catch)]), "): estimated with an F random walk\n", sep = "")
  if(!ncol(x$I)) cat("  catch-only: no abundance index\n")
  for(i in seq_len(ncol(x$I))){
    o <- which(!is.na(x$I[, i]))
    cat("  index ", colnames(x$I)[i], ": ", length(o), " obs, ", x$years[min(o)], "-", x$years[max(o)],
        if(any(!is.na(x$SE[o, i]))) ", with SE" else "", "\n", sep = "")
  }
  for(i in seq_len(NCOL(x$A))){
    if(is.null(x$A)) break
    o <- which(!is.na(x$A[, i]))
    cat("  auxiliary ", colnames(x$A)[i], " (", x$A.type[i], "): ", length(o), " obs, ", x$years[min(o)], "-", x$years[max(o)],
        if(any(!is.na(x$A.SE[o, i]))) ", with SE" else "", "\n", sep = "")
  }
  a <- x$align
  if(!is.null(a) && any(a$n_used < a$n | (a$type == "catch" & (a$first < a$used_first | a$last > a$used_last))))
    cat("  aligned to the catch years: some input years not used (see jb3_align())\n")
  invisible(x)
}

#' Iccat example stocks (from JABBA)
#'
#' Catch, CPUE and SE for Atlantic bigeye tuna (bet), white marlin (whm), South Atlantic swordfish (swos) and blue
#' marlin (bum), as in the JABBA package.
#'
#' Bigeye also carries two \strong{hypothetical} auxiliary series for the guide chapter on catch-only assessment and
#' auxiliary data. They are \strong{not real data}: both are simulated from the F trajectory of the JABBA3 bet example
#' fit (Fox, both JPN longline indices; script \code{data-raw/bet_hypothetical.R}):
#' \describe{
#'   \item{\code{bet$effort_hyp}}{hypothetical effort, data.frame(year, effort), 1960-2017:
#'     E_t = F_t / (q0 x 1.02^(t - 1960)) x lognormal noise (SD 0.15), i.e. with 2\% annual effort creep}
#'   \item{\code{bet$ffmsy_hyp}}{hypothetical F/Fmsy index, data.frame(year, ffmsy, se), 1990-2009 only (a period of high
#'     exploitation): F/Fmsy_t x lognormal noise (SD 0.3), se = 0.3}
#' }
#' @format list of stocks, each a list with catch, cpue, se (data frames) and info
"iccat"

#' Example data for the data-exploration (JARA) models
#'
#' Anonymised example data: generic series names, multiplicative lognormal noise and random rescaling (absolute catches
#' and biomass by one common factor, so that fleets and areas still add up); years and missing values as in the source
#' data.
#' \describe{
#'   \item{areas}{a list with \code{catch}: catch history from six fleets/areas (fleet1-fleet6) with missing early years
#'     for two of them; \code{idx} and \code{idx.se}: relative survey indices for four areas (area1-area4) with log-scale
#'     SEs; \code{bio} and \code{bio.se}: swept-area survey biomass for the same four areas (thousands) with log-scale
#'     SEs.}
#'   \item{regions}{survey and CPUE indices from three regions: \code{areaA} and \code{areaB}, surveys of two adjacent
#'     areas; \code{C.survey1}, \code{C.survey2} and \code{C.cpue1}-\code{C.cpue5}, surveys and CPUE series of a third
#'     region.}
#' }
#' @format list
"jara_examples"
