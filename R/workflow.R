#' BioProcDyn workflows: Tier 1, Tier 2 and Selex
#'
#' Runs the workflows of the BioProcDyn protocol, automatically (\code{jabba3(..., workflow = )}) or step by step
#' (\code{jb3_tier1()}, \code{jb3_tier2()} on any fit).
#' \itemize{
#'   \item \code{"tier1"}: production lag at the age at maturity (\code{amat}).
#'   \item \code{"tier2"}: step 1 = constant sdb (FishLife sdb prior); step 2 = dynamic sdb_t = sdb x (B/Bmsy)^-k with the
#'     level tuned on step 1, sdb = sdb1 / sqrt(mean(x^-2k)) over the step-1 B/Bmsy (bounded to 0.2-3), k estimated with
#'     the life-history prior, plus the lag.
#'   \item \code{"selex"}, \code{"selex_tier1"}, \code{"selex_tier2"}: the same with the spawning-biomass state and
#'     EB/SB-mapped indices and removals; needs \code{model$ebsb} and priors on the SSB scale.
#' }
#' @param data,priors,model as \code{\link{jabba3}()}
#' @param workflow one of "tier1", "tier2", "selex", "selex_tier1", "selex_tier2"
#' @param amat age at 50\% maturity (lag)
#' @param ... passed to \code{\link{jabba3}()}
#' @return the final fit, with all steps in \code{$steps}
#' @export
jb3_workflow <- function(data, priors, model, workflow, amat = NULL, ...){
  sel <- grepl("^selex", workflow)
  if(sel && is.null(model$ebsb)) stop("Selex workflows need EB/SB curves: jb3_model(ebsb = ...)", call. = FALSE)
  if(!sel) model$ebsb <- NULL
  if(grepl("tier", workflow) && is.null(amat)) stop("workflow '", workflow, "' needs amat (age at maturity)", call. = FALSE)
  base <- if(sel) "Selex" else "Current"
  f0 <- jabba3(data, priors, model, run = base, ...)
  if(workflow %in% c("selex")) { f0$steps <- list(f0); names(f0$steps) <- base; return(f0) }
  f1 <- if(grepl("tier1$", workflow)) jb3_tier1(f0, amat, ...) else jb3_tier2(f0, amat, ...)
  f1$steps <- stats::setNames(list(f0, f1), c(base, f1$run))
  f1
}

#' @rdname jb3_workflow
#' @param fit a \code{jabba3} fit (step 1)
#' @export
jb3_tier1 <- function(fit, amat, ...){
  mod <- fit$model; mod$lag <- as.integer(round(amat))
  jabba3(fit$data, fit$priors, mod, run = paste0(if(is.null(mod$ebsb)) "" else "Selex + ", "Tier 1"), ...)
}

#' @rdname jb3_workflow
#' @export
jb3_tier2 <- function(fit, amat = NULL, ...){
  if(fit$model$sdb.form != "const") stop("Tier 2 needs a step-1 fit with constant sdb", call. = FALSE)
  k <- fit$priors$kdyn[["mean"]]
  x <- pmin(pmax(fit$report$BBmsy, 0.2), 3)
  lev <- fit$report$sdb/sqrt(mean(x^(-2*k)))
  mod <- fit$model; mod$sdb.form <- "dyn"; mod$sdb.fixed <- lev
  if(!is.null(amat)) mod$lag <- as.integer(round(amat))
  f <- jabba3(fit$data, fit$priors, mod, run = paste0(if(is.null(mod$ebsb)) "" else "Selex + ", "Tier 2"), ...)
  f$tier2 <- c(sdb.step1 = fit$report$sdb, k.prior = k, sdb.level = lev)
  f
}
