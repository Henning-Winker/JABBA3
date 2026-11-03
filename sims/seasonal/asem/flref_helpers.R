# from FLRef (R/pesim.R): vbbrp() and asem2spm(), to avoid the FLasher/mse dependency
vbbrp <- function(x) {
  
  sel <- landings.sel(x) + discards.sel(x)
  sel <- sel / max(sel)
  
  vb <- quantSums(stock.n(x) %*% stock.wt(x) %*% sel)
  units(vb) <- units(stock(x))
  
  return(vb)
}

asem2spm <- function(object,
                     quant = c("vb", "ssb"),
                     fmsy = NULL,
                     rel = FALSE,
                     spcurve = FALSE) {
  
  quant <- match.arg(quant)
  
  pbrp <- object
  brp <- object
  
  if (is.null(fmsy)) {
    fmsy <- refpts(object)["msy", "harvest"]
  }
  
  fbar(brp) <- FLQuant(c(0, fmsy))
  
  MSY <- an(catch(brp)[, 2])
  
  if (quant == "vb") {
    bio <- vbbrp(brp)
  } else if (quant == "ssb") {
    bio <- ssb(brp)
  }
  
  Hmsy <- an(MSY / bio[, 2])
  BmsyK <- bio[, 2] / bio[, 1]
  
  mi <- seq(0.0001, 5, 0.0001)
  
  m <- mi[
    abs(mi^(-1 / (mi - 1)) - an(BmsyK)) ==
      min(abs(mi^(-1 / (mi - 1)) - an(BmsyK)))
  ][1]
  
  r <- Hmsy * (m - 1) / (1 - 1 / m)
  
  pars <- FLPar(
    r = r,
    m = m,
    BmsyK = BmsyK,
    Hmsy = Hmsy,
    MSY = MSY,
    Bmsy = bio[, 2],
    B0 = bio[, 1],
    Fmsy = fmsy
  )
  
  if (spcurve) {
    
    if (quant == "vb") {
      out <- data.frame(
        Biomass = c(vbbrp(pbrp)),
        SP = c(catch(pbrp)),
        quant = quant
      )
    } else if (quant == "ssb") {
      out <- data.frame(
        Biomass = c(ssb(pbrp)),
        SP = c(catch(pbrp)),
        quant = quant
      )
    }
    
    if (rel) {
      out$Biomass <- out$Biomass / an(pars["B0"])
      out$SP <- out$SP / an(pars["MSY"])
    }
    
    return(list(pars = pars, curve = out))
  }
  
  return(pars)
}
