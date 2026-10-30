#' JABBA3: The Next Generation of JABBA, Now Without the JAGS
#'
#' Bayesian state-space surplus production models in 'RTMB', with NUTS (tmbstan) for full posteriors and
#' 'FLR'-based inputs and outputs. See \code{DESIGN.md} in the package source for the planned architecture.
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom stats median quantile
"_PACKAGE"

#' Engines and status of this JABBA3 installation
#'
#' Reports the package version and which estimation engines are available.
#' @return invisibly, a list with the version and engine availability
#' @export
jabba3_info <- function(){
  eng <- c(rtmb = requireNamespace("RTMB", quietly = TRUE), nuts = requireNamespace("tmbstan", quietly = TRUE))
  flr <- requireNamespace("FLCore", quietly = TRUE)
  message("JABBA3 ", utils::packageVersion("JABBA3"), " -- hyperdrive ", if(eng[["rtmb"]]) "online" else "offline (install RTMB)",
          "\n  engines: rtmb ", if(eng[["rtmb"]]) "yes" else "no", ", nuts ", if(eng[["nuts"]]) "yes" else "no (install tmbstan)",
          "\n  FLR: ", if(flr) "yes" else "no (install FLCore for FLR inputs/outputs)")
  invisible(list(version = utils::packageVersion("JABBA3"), engines = eng, flr = flr))
}
