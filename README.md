# JABBA3

**The Next Generation of JABBA, Now Without the JAGS**

JABBA3 is the successor of [JABBA](https://github.com/jabbamodel/JABBA) (Just Another Bayesian Biomass
Assessment): Bayesian state-space surplus production models, rebuilt in [RTMB](https://cran.r-project.org/package=RTMB)
and making the jump to hyperspace.

- **Fast:** posterior-mode fits with Laplace-integrated states in well under a second, for simulation testing,
  sensitivity runs and management strategy evaluation.
- **Fully Bayesian when you want it:** NUTS through [tmbstan](https://cran.r-project.org/package=tmbstan) on the same
  model, for final runs and credible intervals.
- **FLR inside:** inputs and outputs follow [FLR](https://flr-project.org) conventions (FLQuant, FLIndices, FLPar), so
  fits plot with ggplotFL, feed `FLRef::plotAdvice()` and plug into `mse` as an estimator.
- **No JAGS required.** Classic JABBA inputs (catch, cpue and se data frames) convert with a helper, so old scripts
  survive the trip.

JABBA 2.x (the JAGS version) stays available and maintained as the reference.

> **Status:** early development. The architecture and open decisions are in [`DESIGN.md`](DESIGN.md).
> The estimation core is being ported from the validated RTMB engine for JABBA 2
> (BioProcDyn, JABBA patches 0013-0015).

## Installation

```r
# install.packages("remotes")
remotes::install_github("Henning-Winker/JABBA3")
JABBA3::jabba3_info()
```

## Lineage

| Version | Engine | Notes |
|---|---|---|
| JABBA 1.x | JAGS | Winker et al. (2018) |
| JABBA 2.x | JAGS (RTMB engine as an extension) | process-error options, JABBA-Select |
| **JABBA3** | RTMB + NUTS | FLR interface, ggplot2 graphics (from JAGG) |
