#' Assessment data for JABBA3
#'
#' Collects catch and abundance indices into the internal data object used by \code{\link{jabba3}()}.
#' Native inputs follow FLR: \code{catch} as an \code{FLQuant} (or \code{FLQuants} by fleet, summed) and
#' \code{indices} as \code{FLIndices} of \code{FLIndexBiomass}, with the log-scale SE of each observation in
#' \code{index.var} (NA or empty = not given). Plain data frames in the JABBA 2 layout are accepted as well
#' (see also \code{\link{jb3_data_jabba}()}).
#'
#' @param catch \code{FLQuant}, \code{FLQuants}, or data.frame(year, catch[, more fleets])
#' @param indices \code{FLIndices}, \code{FLIndexBiomass}, or data.frame(year, index1, index2, ...)
#' @param se optional data.frame(year, se1, se2, ...) of log-scale standard errors (data.frame input only)
#' @param stock stock name
#' @return object of class \code{jb3data}: \code{years}, \code{catch}, index matrix \code{I}, SE matrix \code{SE}
#'   (years x indices), \code{stock}
#' @examples
#' data(iccat)
#' d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
#' d
#' @export
jb3_data <- function(catch, indices, se = NULL, stock = "stock"){
  cy <- .catch_df(catch)
  ix <- .index_df(indices, se)
  years <- seq(min(cy$year), max(cy$year))
  if(any(!years %in% cy$year) || any(!is.finite(cy$catch[match(years, cy$year)])))
    stop("catch must be given for every year from ", min(years), " to ", max(years), call. = FALSE)
  C <- cy$catch[match(years, cy$year)]
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
                 catch.fleets = cy$fleets), class = "jb3data")
}

#' Convert classic JABBA inputs
#'
#' The JABBA 2 data frames (\code{catch}, \code{cpue}, \code{se} as used by \code{build_jabba()}) to a
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
  for(i in seq_len(ncol(x$I))){
    o <- which(!is.na(x$I[, i]))
    cat("  index ", colnames(x$I)[i], ": ", length(o), " obs, ", x$years[min(o)], "-", x$years[max(o)],
        if(any(!is.na(x$SE[o, i]))) ", with SE" else "", "\n", sep = "")
  }
  invisible(x)
}

#' Iccat example stocks (from JABBA)
#'
#' Catch, CPUE and SE for Atlantic bigeye tuna (bet), white marlin (whm), South Atlantic swordfish (swos) and blue
#' marlin (bum), as in the JABBA package.
#' @format list of stocks, each a list with catch, cpue, se (data frames) and info
"iccat"
