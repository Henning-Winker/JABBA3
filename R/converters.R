#' Convert classic JABBA data frames to FLR objects
#'
#' Converts the JABBA input data frames (as used by \code{JABBA::build_jabba()}) to FLR objects: catch to an
#' \code{FLQuant} (or \code{FLQuants}, one per fleet), abundance indices to \code{FLIndices} of \code{FLIndexBiomass}
#' with the log-scale standard errors in \code{index.var}, and auxiliary data (e.g. effort) to \code{FLQuants}.
#' The result can be passed straight to \code{\link{jb3_data}()}.
#'
#' @param catch data.frame(year, catch[, catch fleet 2, ...])
#' @param index data.frame(year, index1, index2, ...) (JABBA \code{cpue})
#' @param index.se optional data.frame(year, se1, se2, ...) of log-scale standard errors (JABBA \code{se}); same columns as
#'   \code{index}
#' @param auxiliary optional data.frame(year, aux1, ...) (JABBA \code{auxiliary}, e.g. effort)
#' @param auxiliary.se optional data.frame of log-scale standard errors for \code{auxiliary}
#' @param units.catch,units.index units recorded in the FLR objects
#' @param fleets return catch by fleet as \code{FLQuants} (default: TRUE if more than one catch column)
#' @return list with \code{catch} (FLQuant or FLQuants), \code{indices} (FLIndices) and \code{auxiliary} (FLQuants with
#'   the SEs in attribute "se"; \code{NULL} when not given)
#' @examples
#' \dontrun{
#' data(iccat)
#' x <- jbdat2flr(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se)
#' x$indices
#' d <- jb3_data(x$catch, x$indices, stock = "bet")
#' }
#' @export
jbdat2flr <- function(catch, index, index.se = NULL, auxiliary = NULL, auxiliary.se = NULL,
                      units.catch = "t", units.index = "NA", fleets = NULL){
  .need("FLCore", "jbdat2flr() needs FLCore")
  list(catch = jbcatch2flq(catch, units = units.catch, fleets = fleets),
       indices = jbidx2FLIndices(index, index.se, units = units.index),
       auxiliary = if(is.null(auxiliary)) NULL else .df2flqs(auxiliary, auxiliary.se, units = "NA"))
}

#' @rdname jbdat2flr
#' @param units units of the FLQuant(s)
#' @export
jbcatch2flq <- function(catch, units = "t", fleets = NULL){
  .need("FLCore", "jbcatch2flq() needs FLCore")
  catch <- as.data.frame(catch); yrs <- as.numeric(catch[[1]]); x <- as.matrix(catch[, -1, drop = FALSE])
  if(is.null(fleets)) fleets <- ncol(x) > 1
  mk <- function(v) FLCore::FLQuant(as.numeric(v), dimnames = list(age = "all", year = yrs), units = units)
  if(!fleets) return(mk(rowSums(x, na.rm = TRUE)))
  FLCore::FLQuants(stats::setNames(lapply(seq_len(ncol(x)), function(j) mk(x[, j])), colnames(x)))
}

#' @rdname jbdat2flr
#' @export
jbidx2FLIndices <- function(index, index.se = NULL, units = "NA"){
  .need("FLCore", "jbidx2FLIndices() needs FLCore")
  index <- as.data.frame(index); yrs <- as.numeric(index[[1]]); nm <- colnames(index)[-1]
  if(!is.null(index.se)){
    index.se <- as.data.frame(index.se)
    if(ncol(index.se) != ncol(index)) stop("index.se must have the same columns as index", call. = FALSE)
    index.se <- index.se[match(yrs, as.numeric(index.se[[1]])), , drop = FALSE]
  }
  FLCore::FLIndices(stats::setNames(lapply(seq_along(nm), function(j){
    v <- as.numeric(index[[j + 1]]); ok <- which(is.finite(v) & v > 0)
    yr <- seq(min(yrs[ok]), max(yrs[ok]))                                   # trim leading/trailing NA years
    vv <- v[match(yr, yrs)]
    idx <- FLCore::FLQuant(vv, dimnames = list(age = "all", year = yr), units = units)
    se <- if(is.null(index.se)) FLCore::FLQuant(NA_real_, dimnames = list(age = "all", year = yr)) else
      FLCore::FLQuant(as.numeric(index.se[[j + 1]])[match(yr, yrs)], dimnames = list(age = "all", year = yr))
    fi <- FLCore::FLIndexBiomass(index = idx)
    FLCore::index.var(fi) <- se; FLCore::name(fi) <- nm[j]
    fi
  }), nm))
}

#' Convert FLR inputs back to classic JABBA data frames
#'
#' The inverse of \code{\link{jbdat2flr}()}: an \code{FLQuant}/\code{FLQuants} catch and \code{FLIndices} to the
#' JABBA data frames \code{catch}, \code{cpue} and \code{se}.
#' @param catch FLQuant or FLQuants
#' @param indices FLIndices (or a single FLIndexBiomass)
#' @return list(catch, cpue, se) of data frames
#' @export
flr2jbdat <- function(catch, indices){
  .need("FLCore", "flr2jbdat() needs FLCore")
  cl <- .catch_df(catch); ix <- .index_df(indices, NULL)
  yrs <- seq(min(cl$year, ix$I$year), max(cl$year, ix$I$year))
  cdf <- data.frame(year = yrs, catch = cl$catch[match(yrs, cl$year)])
  if(inherits(catch, "FLQuants")) cdf <- data.frame(year = yrs, sapply(as.list(catch), function(q)
    as.numeric(q[1, , 1, 1, 1, 1])[match(yrs, as.numeric(dimnames(q)$year))]))
  I <- ix$I[match(yrs, ix$I$year), ]; I$year <- yrs; S <- ix$SE[match(yrs, ix$SE$year), ]; S$year <- yrs
  rownames(I) <- rownames(S) <- NULL
  list(catch = cdf, cpue = I, se = S)
}

.df2flqs <- function(x, se = NULL, units = "NA"){
  x <- as.data.frame(x); yrs <- as.numeric(x[[1]])
  q <- FLCore::FLQuants(stats::setNames(lapply(colnames(x)[-1], function(n)
    FLCore::FLQuant(as.numeric(x[[n]]), dimnames = list(age = "all", year = yrs), units = units)), colnames(x)[-1]))
  if(!is.null(se)){ se <- as.data.frame(se)
    attr(q, "se") <- FLCore::FLQuants(stats::setNames(lapply(colnames(se)[-1], function(n)
      FLCore::FLQuant(as.numeric(se[[n]])[match(yrs, as.numeric(se[[1]]))], dimnames = list(age = "all", year = yrs))), colnames(x)[-1])) }
  q
}
