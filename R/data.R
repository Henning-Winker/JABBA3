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
#' @details \code{indices = NULL} gives a catch-only data set: the model is then informed by the catch, the priors
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
                     auxiliary.type = "effort"){
  if(inherits(catch, "jb3catchfill")){ if(is.null(catch.se)) catch.se <- data.frame(year = catch$year, se = catch$se.log)
    catch <- data.frame(year = catch$year, catch = catch$catch) }
  if(inherits(indices, "jb3jointindex") && is.null(se)){ se <- data.frame(year = indices$year, joint = indices$log.se)
    indices <- data.frame(year = indices$year, joint = indices$index) }
  cy <- .catch_df(catch)
  years <- seq(min(cy$year), max(cy$year))
  if(any(!years %in% cy$year) || any(!is.finite(cy$catch[match(years, cy$year)])))
    stop("catch must be given for every year from ", min(years), " to ", max(years), call. = FALSE)
  C <- cy$catch[match(years, cy$year)]
  cse <- NULL
  if(!is.null(catch.se)){
    cse <- if(is.data.frame(catch.se)) as.numeric(catch.se[[2]])[match(years, as.numeric(catch.se[[1]]))] else rep(as.numeric(catch.se), length.out = length(years))
  }
  aux <- .aux_df(auxiliary, auxiliary.se, auxiliary.type, years)
  if(is.null(indices)){                                   # catch-only
    I <- SE <- matrix(NA_real_, length(years), 0, dimnames = list(years, NULL))
    return(structure(list(stock = stock, years = years, catch = C, I = I, SE = SE, catch.fleets = cy$fleets, catch.se = cse,
                          A = aux$A, A.SE = aux$SE, A.type = aux$type), class = "jb3data"))
  }
  ix <- .index_df(indices, se)
  out_of <- setdiff(ix$I$year, years)
  if(length(out_of)) warning("index years outside the catch years dropped: ", paste(out_of, collapse = ", "), call. = FALSE)
  I <- as.matrix(ix$I[match(years, ix$I$year), -1, drop = FALSE])
  SE <- as.matrix(ix$SE[match(years, ix$SE$year), -1, drop = FALSE])
  I[!is.finite(I) | I <= 0] <- NA
  SE[!is.finite(SE) | SE <= 0] <- NA
  dimnames(I) <- dimnames(SE) <- list(years, colnames(ix$I)[-1])
  keep <- colSums(!is.na(I)) > 0
  if(!any(keep)) stop("no index observations within the catch years", call. = FALSE)
  if(any(!keep)) warning("indices without observations dropped: ", paste(colnames(I)[!keep], collapse = ", "), call. = FALSE)
  structure(list(stock = stock, years = years, catch = C, I = I[, keep, drop = FALSE], SE = SE[, keep, drop = FALSE],
                 catch.fleets = cy$fleets, catch.se = cse, A = aux$A, A.SE = aux$SE, A.type = aux$type), class = "jb3data")
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
  if(inherits(catch, "FLQuants")) catch <- Reduce(`+`, as.list(catch))
  if(inherits(catch, "FLQuant")){
    y <- as.numeric(dimnames(catch)$year)
    if(dim(catch)[6] > 1) warning("catch has iterations: the first is used", call. = FALSE)
    return(list(year = y, catch = as.numeric(catch[, , , , , 1]), fleets = 1))
  }
  catch <- as.data.frame(catch)
  if(ncol(catch) < 2) stop("catch needs columns year and catch", call. = FALSE)
  list(year = as.numeric(catch[[1]]), catch = rowSums(as.matrix(catch[, -1, drop = FALSE])), fleets = ncol(catch) - 1)
}

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
  I[[1]] <- as.numeric(I[[1]]); names(I)[1] <- "year"
  if(is.null(se)){ S <- I; S[, -1] <- NA } else {
    S <- as.data.frame(se); S[[1]] <- as.numeric(S[[1]]); names(S)[1] <- "year"
    if(ncol(S) != ncol(I)) stop("se must have the same columns as the indices", call. = FALSE)
    S <- S[match(I$year, S$year), ]; S$year <- I$year; names(S) <- names(I)
  }
  list(I = I, SE = S)
}

#' @export
print.jb3data <- function(x, ...){
  cat("JABBA3 data:", x$stock, "\n  years ", min(x$years), "-", max(x$years), " (", length(x$years), "), catch range ",
      signif(min(x$catch), 3), "-", signif(max(x$catch), 3), "\n", sep = "")
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
