suppressMessages(library(JABBA3, lib.loc = "NULL"))
src <- readLines("sim_seasonal.R"); i0 <- grep("^sim_om <-", src); i1 <- grep("^runs <-", src) - 1
eval(parse(text = src[i0:i1])); `%||%` <- function(a, b) if(is.null(a)) b else a
sc <- list(pulse_unifF = list(w = c(0,3,1,0), fs = rep(.25,4)), flat_unifF = list(w = rep(1,4), fs = rep(.25,4)),
           pulse_seasF = list(w = c(0,3,1,0), fs = c(.1,.4,.4,.1)))
for(k in names(sc)){
  r <- t(sapply(1:40, function(i){ om <- sim_om(i, w = sc[[k]]$w, fshare = sc[[k]]$fs, hmax = 1.5)
    c(minBB = min(om$truth$BBmsy), annual = fitone(om, 1)[1], F = fitone(om, 4, "F")[1], F_pulse = fitone(om, 4, "F", c(0,3,1,0))[1]) }))
  cat(k, ": median terminal error ", paste(colnames(r), round(apply(r, 2, median, na.rm = TRUE), 3), collapse = "  "), "\n")
}
