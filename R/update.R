#' Refit a JABBA3 model with changed inputs
#'
#' Reruns \code{\link{jabba3}()} from a fit, using the data, priors, model and fit settings stored in the fit, with the
#' changes given. Priors and model options are changed by name, with the same arguments as \code{\link{jb3_priors}()} and
#' \code{\link{jb3_model}()}.
#'
#' Some priors are derived from others when \code{jb3_priors()} is called: the FishLife prior on the process-error SD
#' (\code{sdb = "fishlife"}) and the default \code{kdyn} from \code{r} (and \code{sigR}), the default \code{K} from the
#' catch, and \code{K.max = "r"} from \code{r}. By default (\code{recompute = FALSE}) only the priors named in
#' \code{priors} change and all others keep their values in the fit, so changing \code{r} leaves the \code{sdb} prior as it
#' was; set the dependent priors explicitly to control them (e.g. \code{priors = list(r = c(0.3, 0.5), sdb = c(0.15, 0.3))}).
#' With \code{recompute = TRUE} the priors are rebuilt by \code{jb3_priors()} from its stored arguments with the changes,
#' so all priors that were derived (not set explicitly) are derived again from the new values and the (new) data.
#'
#' Model options are always rebuilt with \code{jb3_model()} (quantities such as the shape parameter follow
#' \code{shape} and \code{BmsyK}). For a fit from a workflow (\code{jabba3(..., workflow = )}) the inputs before the
#' workflow steps are used and the workflow is rerun.
#'
#' @param object a \code{jabba3} fit
#' @param data optional new \code{\link{jb3_data}()} object (default: the data of the fit)
#' @param priors a \code{jb3_priors()} object (used as is), or a named list of \code{jb3_priors()} arguments to change
#' @param model a \code{jb3_model()} object (used as is), or a named list of \code{jb3_model()} arguments to change
#' @param recompute rebuild the priors so that derived priors follow the changes (default \code{FALSE}; see Details)
#' @param ... other \code{jabba3()} arguments to change (e.g. \code{engine}, \code{ndraws}, \code{osa}, \code{run});
#'   unchanged settings are taken from the fit
#' @return a new \code{jabba3} fit
#' @examples
#' \donttest{
#' data(iccat)
#' d <- jb3_data(iccat$bet$catch, iccat$bet$cpue, iccat$bet$se, stock = "bet")
#' fit <- jabba3(d, osa = FALSE)
#' f1 <- update(fit, priors = list(r = c(0.3, 0.5)))                    # sdb prior unchanged
#' f2 <- update(fit, priors = list(r = c(0.3, 0.5)), recompute = TRUE)  # sdb prior from the new r
#' f3 <- update(fit, model = list(lag = 3, rho = 0.75), run = "lag + AR(1)")
#' }
#' @importFrom stats update
#' @export
update.jabba3 <- function(object, data = NULL, priors = NULL, model = NULL, recompute = FALSE, ...){
  dat <- data %||% object$data
  pr0 <- object$inputs$priors %||% object$priors
  mo0 <- object$inputs$model %||% object$model
  set <- object$settings
  if(is.null(set)){
    warning("fit has no stored settings (fitted with an earlier JABBA3 version): jabba3() defaults are used for engine, ",
            "ndraws, osa, workflow etc. unless given", call. = FALSE)
    set <- list()
  }
  pr <- .update_priors(pr0, priors, dat, recompute)
  mo <- .update_model(mo0, model)
  args <- .replace_args(set, list(...))
  args$run <- args$run %||% object$run
  do.call(jabba3, c(list(data = dat, priors = pr, model = mo), args))
}

# TRUE if an object built by f() was changed by hand after it was built
.edited <- function(x, rebuilt){
  strip <- function(o){ attr(o, "args") <- NULL; unclass(o) }
  !isTRUE(all.equal(strip(x), strip(rebuilt)))
}

.update_priors <- function(pr0, priors, data, recompute){
  if(is.null(priors) && !recompute) return(pr0)
  if(inherits(priors, "jb3priors")) return(priors)
  if(!is.null(priors) && (!is.list(priors) || is.null(names(priors)) || any(names(priors) == "")))
    stop("priors must be a jb3_priors() object or a named list of jb3_priors() arguments", call. = FALSE)
  bad <- setdiff(names(priors), setdiff(names(formals(jb3_priors)), "data"))
  if(length(bad)) stop("unknown jb3_priors() argument(s): ", paste(bad, collapse = ", "), call. = FALSE)
  a0 <- attr(pr0, "args")
  if(recompute){
    if(is.null(a0)){
      warning("the priors of the fit carry no jb3_priors() arguments (earlier JABBA3 version): recompute = TRUE uses ",
              "jb3_priors() defaults for all priors not given", call. = FALSE)
      a0 <- list()
    } else if(.edited(pr0, do.call(jb3_priors, c(list(data = data), a0))))
      warning("the priors of the fit were edited after jb3_priors(); recompute = TRUE rebuilds them from the jb3_priors() ",
              "arguments, so these edits are lost", call. = FALSE)
    return(do.call(jb3_priors, c(list(data = data), .replace_args(a0, priors %||% list()))))
  }
  # change only the named priors; all others keep their values
  new <- do.call(jb3_priors, c(list(data = data), .replace_args(a0 %||% list(), priors)))
  fields <- unique(unlist(lapply(names(priors), function(n) switch(n, sdb = c("sdb", "sdb.type"), psi = "bkinit", n))))
  out <- pr0
  for(f in fields) out[f] <- list(new[[f]])
  attr(out, "args") <- .replace_args(a0 %||% list(), priors)
  out
}

.update_model <- function(mo0, model){
  if(is.null(model)) return(mo0)
  if(inherits(model, "jb3model")) return(model)
  if(!is.list(model) || is.null(names(model)) || any(names(model) == ""))
    stop("model must be a jb3_model() object or a named list of jb3_model() arguments", call. = FALSE)
  bad <- setdiff(names(model), names(formals(jb3_model)))
  if(length(bad)) stop("unknown jb3_model() argument(s): ", paste(bad, collapse = ", "), call. = FALSE)
  a0 <- attr(mo0, "args")
  if(is.null(a0) || .edited(mo0, do.call(jb3_model, a0))){
    if(!is.null(a0)) warning("the model of the fit was edited after jb3_model(): the changes are applied to the stored ",
                             "model directly (derived options are not updated)", call. = FALSE)
    out <- .replace_args(unclass(mo0), model); attr(out, "args") <- NULL
    return(structure(out, class = "jb3model"))
  }
  do.call(jb3_model, .replace_args(a0, model))
}

# replace (not merge) named elements; NULL values are kept (jabba3()/jb3_priors() defaults)
.replace_args <- function(a, b){
  for(n in names(b)) a[n] <- list(b[[n]])
  a
}
