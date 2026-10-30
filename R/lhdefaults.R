#' Life-history guidance for the AR(1) persistence of process error
#'
#' The persistence of process deviations in exploitable biomass increases with generation time (GT). Across age-structured
#' simulations for 486 bony-fish species of the FAO SOFIA list (FishLife life histories, F-contrast trajectories) the lag-1
#' autocorrelation followed \eqn{\rho \approx \exp(-c/GT)} with c = 2.25 (R^2 = 0.81; c = 2.14 with size-dependent
#' selectivity, c = 2.35 in nine ICES case-study operating models). The default \code{rho = 0.75} corresponds to
#' GT of about 6-10 years; use this function to adjust it for shorter- or longer-lived stocks.
#' GT can be computed from the same Leslie matrix as the r prior (e.g. \code{JABBA::jbleslie()} or
#' \code{FLSRTMB::productivity()}).
#' @param gt generation time (years)
#' @param c decay constant (default 2.2)
#' @param digits rounding
#' @return suggested rho
#' @examples
#' jb3_rho_gt(c(2, 4, 8, 15))
#' jb3_rho_table()
#' @export
jb3_rho_gt <- function(gt, c = 2.2, digits = 2) round(exp(-c/gt), digits)

#' @rdname jb3_rho_gt
#' @return \code{jb3_rho_table()}: guidance table by generation-time class with the FishLife calibration medians
#' @export
jb3_rho_table <- function(){
  data.frame(`Generation time (years)` = c("<= 3", "3-6", "6-10", "> 10"),
             `Case-study stocks (GT)` = c("very short-lived species",
                                          "Bay of Biscay sardine (3.2), northern shrimp (3.8)",
                                          "Baltic herring (6.6), North Sea plaice (6.7), southern hake (8.4)",
                                          "Icelandic cod (11.1), western horse mackerel (11.4), Biscay sea bass (11.9), North Sea anglerfish (12.6)"),
             `Calibration median` = c(0.30, 0.63, 0.78, 0.80),
             `Suggested rho` = c("0.3-0.4", "0.6", "0.75 (default)", "0.8"),
             check.names = FALSE)
}
