# =====================================================================================================================
# JABBA3 user guide as an R script
# ---------------------------------------------------------------------------------------------------------------------
# Generated from inst/guide/chapters/*.Rmd with knitr::purl(documentation = 2): the guide text is kept as #' comments,
# chapter and section headings as "#' #" lines, and every code chunk is marked "## ---- <chunk name> ----".
# Chunks that the guide shows but does not run (eval = FALSE: installation, JABBA/JAGS comparisons, long runs) are
# commented out with "## ".
#
# Requirements: JABBA3 (remotes::install_github("Henning-Winker/JABBA3")), ggplot2, patchwork; NUTS sections need
# tmbstan and rstan; FLR sections need FLCore, ggplotFL and FLRef; converter sections need JABBA.
# Run section by section (RStudio: Ctrl+Enter, or Ctrl+Alt+T for the current chunk); objects are reused across
# chapters (e.g. bet, d, pr, fit), so run from the top.
# =====================================================================================================================

## ----setup--------------------------------------------------------------------
library(JABBA3)
library(ggplot2)

#' 
#' # JABBA3: surplus production models with RTMB speed
#' 
#' JABBA3 fits **Bayesian state-space surplus production models**: biomass
#' dynamics driven by a Pella-Tomlinson production function, fishery removals and
#' process error, fitted to one or more abundance indices. It is the successor of
#' JABBA (Just Another Bayesian Biomass Assessment; Winker et al. 2018), rebuilt
#' in RTMB (Kristensen et al. 2016) and with inputs and outputs in FLR (Kell et
#' al. 2007):
#' 
#' - **The same model family as JABBA.** Schaefer, Fox and Pella-Tomlinson
#'   production, the JABBA priors (r, K, initial depletion, shape), process error
#'   after the catch, additive index variance and catch with error. With the same
#'   priors JABBA3 reproduces the JABBA RTMB engine to within $10^{-4}$ in
#'   $\log B/B_{MSY}$.
#' - **Fast.** The log-depletion states are integrated by the Laplace
#'   approximation and the fixed effects estimated at the posterior mode; a fit
#'   takes well under a second, so retrospectives, hindcasts, sensitivity runs and
#'   simulation testing are cheap.
#' - **Fully Bayesian when needed.** The same objective can be sampled with NUTS
#'   (tmbstan; Monnahan and Kristensen 2018) for final runs and credible
#'   intervals.
#' - **Process error from BioProcDyn.** Production lag at the age at maturity,
#'   AR(1) process deviations, depletion-dependent process error and a
#'   life-history (FishLife) prior on its scale.
#' - **Diagnostics.** One-step-ahead (OSA) residuals by particle filter, runs
#'   tests, the SPiCT residual tests, retrospective analysis (Mohn's rho) and
#'   hindcast cross-validation (MASE).
#' - **FLR inside.** Catch as `FLQuant`, indices as `FLIndices`; status, reference
#'   points and parameters out as `FLQuants` and `FLPar`.
#' - **Converters.** Classic JABBA data frames, JABBA `build_jabba()` objects
#'   and SPiCT input lists convert to JABBA3 inputs.
#' 
#' > **JABBA and JABBA3.** Since version 2.x, JABBA itself has an RTMB engine in
#' > addition to JAGS, which few users are aware of. JABBA3 is the standalone RTMB
#' > successor. In this guide "JABBA" means the original package (JAGS, with its
#' > RTMB engine) and "JABBA3" the new one.
#' 
#' How JABBA options map to JABBA3:
#' 
#' | JABBA | JABBA3 |
#' |---|---|
#' | `build_jabba(catch, cpue, se, ...)` | `jb3_data()` + `jb3_priors()` + `jb3_model()`, or `jbinput2jabba3(jbinput)` |
#' | `fit_jabba()` (JAGS) | `jabba3(engine = "rtmb")` (mode) or `jabba3(engine = "nuts")` |
#' | `r.prior`, `K.prior`, `psi.prior` | `jb3_priors(r = , K = , bkinit = )`, each `c(median, log.sd, on)` |
#' | `sigma.proc.prior` (inverse-gamma) | `jb3_priors(sdb = "fishlife")` (lognormal, life-history centre) |
#' | `model.type`, `BmsyK` | `jb3_model(shape = , BmsyK = )` |
#' | `sigma.est`, `fixed.obsE` | estimated `sdi` + input SEs (`obs.error = "additive"`) |
#' | `add.catch.CV` | `jb3_model(catch.error = TRUE, sdc = 0.1)` |
#' | `jbplot_*()` | `jbplot_*()` (ggplot2) |
#' 
#' ## About this guide
#' 
#' **Part I** (sections 3-9) goes through an assessment from the data to stock
#' status: data, priors, model options, fitting, convergence and residual
#' diagnostics, with **checks** at each step (what to look at, what is a red flag,
#' what to do about it). **Part II** (sections 10-15) covers NUTS, the BioProcDyn
#' process-error workflows, projections, FLR outputs, converters from JABBA and
#' SPiCT, and the model equations. Section 16 summarises the checks and the
#' functions. The examples use the built-in ICCAT data (`data(iccat)`): Atlantic
#' bigeye tuna, white marlin, South Atlantic swordfish and blue marlin.
#' 
#' 
#' # Installation
#' 
#' ## Required
#' 
#' JABBA3 needs `RTMB`, `ggplot2` and `Matrix`; these are installed with the package:
#' 
## ----install, eval = FALSE----------------------------------------------------
## install.packages("remotes")
## remotes::install_github("Henning-Winker/JABBA3")

#' 
#' ## Optional
#' 
#' | Package | Needed for | Install |
#' |---|---|---|
#' | `tmbstan`, `rstan` | NUTS (`engine = "nuts"`): full posteriors for final runs, catch-only and prior predictive runs | `install.packages(c("tmbstan", "rstan"))` |
#' | `patchwork` | multi-panel diagnostics and summary plots | `install.packages("patchwork")` |
#' | `FLCore`, `ggplotFL`, `FLRef` | FLR inputs and outputs (`jb3_FLStockR()`, `jbdat2flr()`) | see below |
#' | `JABBA` | converting JABBA (JAGS) inputs (`jbinput2jabba3()`) | `remotes::install_github("JABBAmodel/JABBA")` |
#' | `MASS` | Kobe contours | `install.packages("MASS")` |
#' 
## ----install-optional, eval = FALSE-------------------------------------------
## install.packages(c("tmbstan", "rstan", "patchwork", "MASS"))
## install.packages(c("FLCore", "ggplotFL", "FLRef"),
##                  repos = c("https://flr.r-universe.dev", "https://cloud.r-project.org"))
## remotes::install_github("JABBAmodel/JABBA")

#' 
#' > **NUTS needs a C++ toolchain.** `rstan` compiles models, so it needs Rtools on Windows
#' > (<https://cran.r-project.org/bin/windows/Rtools/>, the version matching your R) or the Xcode command-line tools
#' > on macOS (`xcode-select --install`). Check the toolchain with `pkgbuild::has_build_tools(debug = TRUE)` before
#' > installing `rstan`. JABBA3 runs without NUTS; only `engine = "nuts"`, the NUTS chapters and catch-only runs
#' > need it.
#' 
#' ## Check the installation
#' 
## ----load---------------------------------------------------------------------
library(JABBA3)
packageVersion("JABBA3")
jabba3_info()
opt <- c("tmbstan", "rstan", "patchwork", "FLCore", "ggplotFL", "FLRef", "JABBA", "MASS")
data.frame(package = opt, installed = vapply(opt, requireNamespace, logical(1), quietly = TRUE), row.names = NULL)

#' 
#' # Quick start: your first fit
#' 
#' Atlantic bigeye tuna is the standard example of this guide: catch from 1950 and
#' two Japanese longline CPUE indices with log-scale standard errors, in the JABBA
#' data-frame format. Atlantic blue marlin, a more complex assessment with ten
#' indices, introduces the rebuilding protocol (chapter *The rebuilding protocol*).
#' 
## ----qs-data------------------------------------------------------------------
data(iccat)
names(iccat)
bet <- iccat$bet
head(bet$catch, 3)
head(bet$cpue[bet$cpue$year >= 1959, ], 3)

#' 
#' Three steps: data, priors, fit.
#' 
## ----qs-fit-------------------------------------------------------------------
d   <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
pr  <- jb3_priors(d, r = c(0.2, 0.5))
fit <- jabba3(d, pr, jb3_model(shape = "fox"), ndraws = 1000)
fit

#' 
#' The print shows the convergence verdict, the estimates with 95% intervals and
#' the stock status in the last year. `ndraws = 1000` adds draws from the
#' multivariate normal approximation at the mode, used for the Kobe plot and
#' projections. The fitting took `r round(fit$time, 2)` s.
#' 
#' ## Quick output: `plot(fit)`
#' 
#' `plot(fit)` (or `jabbaplot3(fit)`) is the one-page summary of a fit, as SPiCT's
#' `plot()`:
#' 
#' - **top row:** biomass (with the observed indices scaled by $q$ as points), $F$
#'   and catch, with MSY reference points and $B_{thr} = 0.5 B_{MSY}$ and
#'   $B_{lim} = 0.3 B_{MSY}$;
#' - **middle row:** $B/B_{MSY}$ (again with the scaled indices), $F/F_{MSY}$ and the
#'   Kobe plot with the uncertainty of the terminal year;
#' - **bottom row:** the process deviations, the realised surplus production
#'   $B_{t+1} - B_t + C_t$ against biomass with the production curve, and the
#'   production phase plot: catches against biomass on the Kobe colours, with the
#'   MSY band.
#' 
## ----qs-plot, fig.height = 9, fig.width = 10----------------------------------
plot(fit)

#' 
#' `plot(fit, "simple")` gives a 2 x 2 version ($B/B_{MSY}$, $F/F_{MSY}$, catch,
#' Kobe; as SPiCT's `plot2()`), and `panels = ` picks any set of panels, e.g.
#' `plot(fit, panels = c("B/Bmsy", "F/Fmsy", "fits", "residuals"))`.
#' 
## ----qs-simple, fig.height = 6, fig.width = 8---------------------------------
plot(fit, "simple")

#' 
#' ## Individual plots
#' 
## ----qs-advice, fig.height = 7, fig.width = 7---------------------------------
jbplot_advice(fit)

#' 
#' The advice plot follows FLRef's `plotAdvice()` and JAGG: the median with
#' nested 90% and 60% bands and the MSY reference points as labelled lines.
#' 
## ----qs-kobe, fig.height = 5--------------------------------------------------
jbplot_kobe(fit)
jbplot_kobe(fit, style = "contours")   # filled 50/80/95% regions, as jbplot_biplot(); contours = FALSE: points only

#' 
#' `jbplot_all(fit, dir)` writes all standard plots of a fit to PNG files.
#' 
#' # Data
#' 
#' ## Data frames (JABBA format)
#' 
#' `jb3_data()` takes the JABBA data frames directly: `catch`
#' (`year, catch[, fleet 2, ...]`), `cpue` (`year, index1, index2, ...`, `NA` for
#' missing years) and optionally `se` (log-scale standard errors, same columns as
#' `cpue`). Catches by fleet are summed for the fit and kept in `catch.fleets`.
#' 
## ----data-df------------------------------------------------------------------
d <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
d$years
str(d$I)

#' 
## ----data-plot, fig.height = 5------------------------------------------------
jbplot_data(d)

#' 
#' `jbplot_indices()` shows each index scaled to its mean with the 95% interval
#' from the input standard errors, as JABBA's plot of the same name:
#' 
## ----data-indices, fig.height = 4---------------------------------------------
jbplot_indices(d)

#' 
#' ## FLR objects
#' 
#' The native input is catch as an `FLQuant` (or `FLQuants` by fleet) and indices
#' as `FLIndices` of `FLIndexBiomass`, with the log-scale standard errors in
#' `index.var`. `jbdat2flr()` converts the JABBA data frames, and `flr2jbdat()`
#' converts back:
#' 
## ----data-flr-----------------------------------------------------------------
library(FLCore)
x <- jbdat2flr(bet$catch, bet$cpue, bet$se)
names(x)                    # auxiliary is NULL when not given
x$catch[, ac(2013:2017)]
names(x$indices)
index.var(x$indices[[2]])[, ac(2013:2017)]
d.flr <- jb3_data(x$catch, x$indices, stock = "bet")
all.equal(d.flr$I, d$I, check.attributes = FALSE)

#' 
#' All functions that take the data work the same way with FLR input; `jbplot_data()`
#' also takes the FLR objects directly:
#' 
## ----data-flr-plot, fig.height = 4--------------------------------------------
jbplot_data(x$catch, x$indices, stock = "bet")

#' 
#' Several catch columns give `FLQuants` (one per fleet) unless `fleets = FALSE`.
#' `jbcatch2flq()` and `jbidx2FLIndices()` convert the catch and the indices
#' separately.
#' 
#' ## Check 1: the data
#' 
#' `jbplot_data()` shows the catch and the indices. Look for:
#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | contrast | indices flat while catch changes, or the reverse | r and K are confounded; tighter r prior, or use status priors |
#' | index overlap | indices in conflict where they overlap | separate `sdi`; drop or down-weight one index; sensitivity runs |
#' | standard errors | very small SEs (< 0.05) | the estimated `sdi` adds variance, but check where the SEs come from |
#' | catch | gaps, early catch reconstructed | `catch.error = TRUE` (default); start later with `bkinit` prior |
#' 
#' # Data exploration: index trends and joint indices
#' 
#' Before fitting a production model it pays to look at the abundance indices on
#' their own: do they tell the same story, how do they scale against each other,
#' and do some of them conflict? JABBA3 includes the trend models of JARA, refitted
#' in RTMB, for this: a state-space random walk with drift on log abundance,
#' observed through one or several indices. For index $i$ in year $t$:
#' 
#' $$x_t = x_{t-1} + r + \eta_t,\quad \eta_t \sim N(0, \sigma_p^2);\qquad
#' \log I_{i,t} = \log q_i + x_t + \varepsilon_{i,t},\quad
#' \varepsilon_{i,t} \sim N(0, se_{i,t}^2 + \sigma_{fix}^2 + \sigma_i^2)$$
#' 
#' The common state $x_t$ is the joint trend; each index has its own catchability
#' $q_i$, relative to the first (reference) index with $q_1 = 1$, and its own
#' additional observation SD $\sigma_i$ (`variance.weighting = "model"`;
#' `"equal"` estimates one common value, `"fixed"` none; `fixed.obsE` is
#' $\sigma_{fix}$). The joint index follows the hierarchical analysis of multiple
#' noisy abundance indices of Conn (2010), widely used in US (SEDAR) assessments:
#' JARA's random-walk form gives practically the same joint index as Conn's
#' independent year effects (`state = "iid"`).
#' 
#' The fits carry the same diagnostics as a `jabba3()` fit (fits, runs tests,
#' joint residuals, posterior predictive checks), with one-step-ahead (OSA)
#' residuals by default and `type = "ordinary"` for the log residuals used in JARA
#' reports.
#' 
#' **Log-scale standard errors.** All models use observation errors on the log
#' scale. From a design-based SD of the untransformed index, use
#' $CV = sd/mean$ and $\sigma = \sqrt{\log(CV^2 + 1)}$ (close to the CV for
#' CV < 0.4); from a 95% interval, $\sigma = (\log(uci) - \log(lci))/(2 \times 1.96)$.
#' 
#' ## Concept: hierarchical analysis of multiple indices
#' 
#' The joint-index models rest on the hierarchical framework of Conn (2010). Each
#' index $i$ is taken to measure the same relative abundance $n_t$, but through its
#' own catchability $c_i = \log q_i$, with two kinds of error: process error
#' $\sigma_{p,i}$, from year-to-year variation in catchability or in the spatial
#' distribution of the stock relative to the fishery or survey, and sampling error
#' $\sigma_{s,it}$, known from the input CV of each value:
#' 
#' $$\log U_{it} \sim N\left(n_t + c_i,\; \sigma_{p,i}^2 + \sigma_{s,it}^2\right)$$
#' 
#' Only differences between the $c_i$ are identifiable, so one constraint is needed:
#' a reference index with $c_1 = 0$ ($q_1 = 1$; JABBA3) or an anchor on the mean of
#' the $c_i$ or $n_t$. Conn treats the $n_t$ as independent year effects; JABBA3
#' (following JARA) places a random walk with drift on $n_t$ (`state = "rw"`), which
#' carries the trend through years with few or no observations, and fits Conn's form
#' with `state = "iid"`. The estimated index-specific $\sigma_{p,i}$ correspond to
#' JABBA3's additional observation SD $\sigma_i$ (`variance.weighting = "model"`).
#' The outputs follow directly: each $q_i$ is the scale of index $i$ *relative to
#' the reference index*, not an absolute catchability; the joint index is the
#' estimated $n_t$ with its uncertainty; and residual patterns that persist for one
#' index against the common trend point to conflicts between indices rather than to
#' noise.
#' 
#' ## Aligning indices: bigeye tuna
#' 
#' The two Japanese longline CPUE series of bigeye cover consecutive periods
#' (1959-1978 and 1979-2017) and do not overlap. A relative trend model aligns
#' them on one consistent trend, each index scaled by its catchability:
#' 
## ----trend-bet----------------------------------------------------------------
tr.bet <- jb3_trend(bet$cpue, bet$se, run = "bigeye")
tr.bet

#' 
## ----trend-bet-plot, fig.height = 3.5-----------------------------------------
jbplot_trend(tr.bet)

#' 
#' `jb3_qtable()` gives the catchability of each index relative to the reference
#' index (the first one, q = 1 by definition) with its 95% interval, and
#' `jbplot_q()` its distribution, with the reference as a vertical line at 1:
#' 
## ----trend-bet-q, fig.height = 3----------------------------------------------
knitr::kable(jb3_qtable(tr.bet, digits = 3), row.names = FALSE)
jbplot_q(tr.bet)

#' 
#' The late series has a catchability about 2.4 times that of the early series.
#' Because the series do not overlap, the scaling rests on the continuity of the
#' trend across 1978/79, and the interval is correspondingly wide. In the
#' production model, the two series keep their own q; the trend model shows how
#' much the data alone say about their relative scale.
#' 
#' ## Why catchabilities differ: changes in vessel, gear and method
#' 
#' Differences in q between series are rarely a nuisance parameter only. In
#' fishery-independent surveys, abrupt changes in catchability follow changes of
#' research vessel or trawl gear. A comparison of design-based swept-area
#' estimates and geostatistical delta-GLMM indices (Thorson et al. 2015) for eight
#' species of a demersal trawl survey (Winker et al. 2017) found clear vessel
#' effects for several species, seen as abrupt steps in the index at the change of
#' research vessel, and a gear effect (lower catchability with the new gear) for one
#' species. Both estimators showed the same steps, so the change was in q, not in
#' the estimator. Such a change should be modelled as a change in q, not read as a
#' change in abundance: split the index at the known change and let the trend model
#' (or the production model) estimate a q for each period.
#' 
#' `jb3_splitindex()` does the split. A hypothetical example: a survey index whose
#' new vessel, from 2009, catches 40% less than the old one (q ratio 0.6).
#' 
## ----trend-split--------------------------------------------------------------
a  <- jara_examples$areas
s  <- a$idx[c("year", "area3")]; se <- a$idx.se[c("year", "area3")]
new.v <- 0.6*s$area3                          # what the new vessel would catch
obs <- s; obs$area3[s$year >= 2009] <- new.v[s$year >= 2009]
# (a) no overlap: the old vessel stops in 2008
x0 <- jb3_splitindex(obs, se, "area3", at = 2009, labels = c("vessel.old", "vessel.new"))
# (b) three calibration years (2009-2011) fished with both vessels
obs1 <- obs; cy <- 2009:2011; obs1$area3[s$year %in% cy] <- s$area3[s$year %in% cy]
x1 <- jb3_splitindex(obs1, se, "area3", at = 2009, overlap = cy, new = new.v[s$year %in% cy],
                     labels = c("vessel.old", "vessel.new"))
q.split <- rbind(cbind(design = "no overlap", jb3_qtable(jb3_trend(x0$I, x0$se), digits = 3)[2, ]),
                 cbind(design = "3 calibration years", jb3_qtable(jb3_trend(x1$I, x1$se), digits = 3)[2, ]))
knitr::kable(q.split[, c("design", "index", "q", "lci", "uci")], row.names = FALSE,
             caption = "q of the new vessel relative to the old (true value 0.6)")

#' 
#' Without overlap the trend model cannot tell a lower q from a decline in the
#' stock: the estimated ratio is close to 1 with a wide interval, and the step
#' appears as a decline in the common trend. Three calibration years recover the
#' ratio (about 0.64, interval excluding 1). The same holds for the bigeye
#' longline series above, whose relative q rests on the continuity of the trend
#' alone. Other indices help to identify a change in q only if they follow the same
#' trend; indices from areas with different trends (chapter *Other data features*)
#' do not.
#' 
#' ## Observation and process noise of indices
#' 
#' Survey indices from a design-based swept-area estimator assume that density is
#' homogeneous within each stratum. For patchy species, or species found in only
#' part of the survey area, the catch rates of randomly placed stations then vary
#' widely, and the index is imprecise and shows spikes that are implausible for the
#' population (Winker et al. 2017). Model-based standardisation, such as the
#' geostatistical delta-GLMM, accounts for the spatial correlation of density and
#' reduces both the observation error and this year-to-year "process noise". In the
#' survey comparison of Winker et al. (2017), the geostatistical indices had lower
#' observation-error CVs for all eight species (e.g. 0.12 against 0.18 for one hake
#' species, 0.27 against 0.39 for a patchy by-catch species), and lower process-error
#' CVs (Francis 2011), with the largest reductions for the patchy by-catch species
#' (0.49 to 0.35) and for the most commonly caught hake (0.53 to 0.46). The mean
#' trends of the two estimators agreed.
#' 
#' For an assessment, process noise in an index is variation in catchability or
#' availability that the sampling error does not capture. It is what the
#' additional SD of each index against the common trend, $\sigma_{p,i}$ (Conn 2010),
#' estimates. `jb3_noise()` sets three measures side by side: the mean input CV, the
#' loess process CV of Francis (2011) from a smoother through the log index, and,
#' for a trend fit, $\sigma_{p,i}$ and the RMSE against the common trend.
#' 
#' ## Index conflicts and a joint index: blue marlin
#' 
#' The 2018 blue marlin assessment used ten CPUE series of different fleets and
#' periods. Before combining them, check whether they agree. A first look:
#' correlations of the log indices over the years each pair shares
#' (`jbplot_indexcor()`; blank where fewer than five years overlap):
#' 
## ----trend-bum-cor, fig.height = 6, fig.width = 7-----------------------------
bum <- iccat$bum
tr.bum <- jb3_trend(bum$cpue, bum$se, run = "blue marlin")
jbplot_indexcor(tr.bum)

#' 
#' Of the 34 pairs with at least five shared years, 14 are negatively
#' correlated, though none significantly, while five positive correlations are
#' significant. Correlations over few shared years are a weak test, so the trend
#' model's residuals are the main check. The joint residual plot shows the
#' residuals of all indices against the common trend, with their RMSE and a loess
#' smoother; conflicting indices show opposite residual patterns in the same
#' years:
#' 
## ----trend-bum-resid, fig.height = 4------------------------------------------
jbplot_residuals(tr.bum)

#' 
#' Runs tests by index (OSA residuals) and posterior predictive checks:
#' 
## ----trend-bum-runs, fig.height = 7, fig.width = 9----------------------------
p.runs <- jbplot_runs(tr.bum)
p.runs
attr(p.runs, "runs")[, c("index", "p.runs", "test")]

#' 
## ----trend-bum-ppc, fig.height = 4--------------------------------------------
p.ppc <- jbplot_ppc(tr.bum)
attr(p.ppc, "ppc")

#' 
#' Four indices fail the runs test against the common trend (`TAI.1st`,
#' `VEN.LL`, `VEN.Rec`, `GHA.GIL`); three of them have among the largest
#' additional observation SDs (`sdi` in the print), i.e. the model already
#' down-weights them. The posterior predictive p-values are all between 0.17 and
#' 0.78, so no index is grossly inconsistent with the joint model. The relative
#' catchabilities follow below. First, the noise of each index:
#' 
## ----trend-bum-noise, fig.height = 3.5----------------------------------------
knitr::kable(jb3_noise(tr.bum), row.names = FALSE)
jbplot_noise(tr.bum)

#' 
#' Three of the four series flagged by the runs tests (`TAI.1st`, `VEN.Rec`,
#' `GHA.GIL`) carry additional SDs of 0.6-1.0 against the common trend, well
#' beyond their input CVs (0.10-0.39). For `VEN.Rec` the loess CV is similarly high
#' (0.60): the series is noisy around its own trend. For `TAI.1st` and `GHA.GIL`
#' the additional SD is two to three times the loess CV: the series are fairly
#' smooth but follow a different trend, which is a conflict with the other indices
#' rather than noise. `VEN.LL`, the fourth, is the reverse case: its input CV (0.57)
#' already exceeds its loess CV and RMSE, so its runs-test failure reflects a
#' residual pattern, not unaccounted noise. The relative catchabilities:
#' 
## ----trend-bum-q, fig.height = 4----------------------------------------------
knitr::kable(jb3_qtable(tr.bum, digits = 3), row.names = FALSE)
jbplot_q(tr.bum)

#' 
#' The joint index combines the ten series into one, the precision-weighted mean
#' of the q-scaled indices by year with its log-scale SE:
#' 
## ----trend-bum-joint, fig.height = 3.5----------------------------------------
bum.joint <- jb3_jointindex(tr.bum, name = "joint")
jbplot_jointindex(tr.bum)
knitr::kable(head(bum.joint, 6), digits = 3, row.names = FALSE)

#' 
#' `bum.joint` is used in the blue marlin case study (chapter *The rebuilding
#' protocol*) to compare an assessment on the joint index with one on the ten
#' separate indices.
#' 
#' **Random walk (JARA) or independent year effects (Conn 2010).** The blue marlin
#' indices overlap enough for Conn's form, which treats the years as independent
#' effects rather than a random walk; the two joint indices are practically the
#' same:
#' 
## ----trend-bum-conn, fig.height = 3.5-----------------------------------------
jc <- jb3_jointindex(jb3_trend(bum$cpue, bum$se, state = "iid", run = "Conn (iid)"))
ok <- is.finite(bum.joint$index) & is.finite(jc$index)
cmp <- data.frame(year = bum.joint$year[ok], rw = bum.joint$index[ok]/mean(bum.joint$index[ok]),
                  iid = jc$index[ok]/mean(jc$index[ok]))
round(cor(log(cmp$rw), log(cmp$iid)), 3)
ggplot2::ggplot(cmp, ggplot2::aes(year)) + ggplot2::geom_line(ggplot2::aes(y = rw, colour = "random walk (JARA)")) +
  ggplot2::geom_point(ggplot2::aes(y = iid, colour = "iid (Conn 2010)")) + ggplot2::labs(y = "Joint index", colour = NULL) +
  theme_jabba3()

#' 
#' For series that do not overlap in time, such as the two bigeye indices, only
#' the random walk can scale them against each other.
#' 
#' Catch reconstruction, swept-area surveys from several areas and checks for
#' conflicts between regions are covered in chapter *Other data features*.
#' 
#' # Priors
#' 
#' All priors are in one list from `jb3_priors()`. Each prior is
#' `c(median, log.sd, on)`: lognormal with the given median and log-scale SD, and
#' an on/off switch (`on = 0` gives a flat prior on the log scale). The defaults
#' follow JABBA (`build_jabba()`), except for the process error.
#' 
## ----priors-------------------------------------------------------------------
pr <- jb3_priors(d, r = c(0.2, 0.5))
pr

#' 
#' | Prior | Default | Note |
#' |---|---|---|
#' | `r` | `c(0.2, 0.5, 1)` | intrinsic growth rate; on by default |
#' | `K` | 8 x max catch, log.sd 0.833 (CV 1) | as JABBA; `c(NA, 0.833, 0)` for flat on log K |
#' | `bkinit` (alias `psi`) | `c(0.9, 0.246, 1)` | initial depletion $B_1/K$ |
#' | `sdb` | `"fishlife"` | process error SD at $B_{MSY}$, see below |
#' | `q` | flat on log q | matrix with one row `c(median, log.sd, on)` per catchability group |
#' | `m` | off | shape fixed from `BmsyK`; `c(median, log.sd, 1)` estimates it |
#' | `rho` | off (`c(0.7, 0.15, 0)`) | normal prior on the AR(1) coefficient, when estimated |
#' | `kdyn` | `max(0, 0.064 - 0.155 log r)`, sd 0.15 | depletion dependence of `sdb` (dynamic SD only) |
#' | `sdi` | off | index SD prior, used with `obs.error = "spict"` |
#' | `status` | `NULL` | data.frame(year, type, median, log.sd); type `"BK"`, `"BBmsy"` or `"FFmsy"` |
#' 
#' ## Process error: the FishLife prior
#' 
#' `sdb = "fishlife"` (default) replaces JABBA's diffuse inverse-gamma prior on
#' the process variance by a lognormal prior with a life-history centre,
#' 
#' $$
#' \sigma_{MSY} \sim LN\left(\log(0.80\,\sigma_R\,r^{0.45}),\ 0.3\right),
#' $$
#' 
#' with $r$ the median of the r prior and $\sigma_R$ the recruitment variability
#' (`sigR = 0.5` when unknown). The relation comes from age-structured
#' simulations of FishLife species (Thorson 2020) with selectivity at maturity
#' (BioProcDyn). The inverse-gamma prior is still available (`sdb = "igamma"`)
#' but its mode is at zero, so it is not recommended for posterior-mode fits; a
#' fixed-scale prior is given as `sdb = c(median, log.sd, 1)`.
#' 
## ----priors-sdb---------------------------------------------------------------
jb3_priors(d, r = c(0.2, 0.5), sigR = 0.8)$sdb
jb3_priors(d, r = c(0.2, 0.5), sdb = c(0.1, 0.3, 1))$sdb

#' 
#' ## r and the shape
#' 
#' The r prior refers to the JABBA parameterisation of the Pella-Tomlinson
#' model, in which $F_{MSY} = r/m$ (section *Model equations*). An r prior
#' derived for the Schaefer model ($m = 2$, $F_{MSY} = r/2$) therefore has to be
#' converted when the shape changes: for the same $F_{MSY}$, the r for shape $m$
#' is $r_m = m\,F_{MSY}$. The JABBA package has the helper
#' `JABBA::r_from_fmsy_bmsyk(Fmsy, BmsyK)`:
#' 
## ----priors-r-----------------------------------------------------------------
if(requireNamespace("JABBA", quietly = TRUE))
  JABBA::r_from_fmsy_bmsyk(Fmsy = 0.1, BmsyK = c(0.5, 0.4, 0.37))

#' 
#' > A Schaefer r of 0.2 ($F_{MSY} = 0.1$) corresponds to r = 0.12 for
#' > $B_{MSY}/K = 0.4$ and r = 0.1 for the Fox model. Leaving the r prior
#' > unchanged while changing the shape changes the implied $F_{MSY}$.
#' 
#' ## Status priors
#' 
#' A prior on depletion, $B/B_{MSY}$ or $F/F_{MSY}$ in any year (e.g. from expert
#' knowledge or a CMSY-type analysis) is added as a data frame:
#' 
## ----priors-status------------------------------------------------------------
pr.st <- jb3_priors(d, r = c(0.2, 0.5),
                    status = data.frame(year = 2010, type = "BBmsy", median = 1, log.sd = 0.3))
pr.st$status

#' 
#' After fitting, `jbplot_bprior()` compares each status prior with the posterior
#' of that quantity in that year, with the posterior-to-prior mean and variance
#' ratios (PPMR, PPVR) as in JABBA:
#' 
## ----priors-bprior, fig.height = 3.5------------------------------------------
fit.st <- jabba3(d, pr.st, jb3_model(shape = "fox"), osa = FALSE)
jbplot_bprior(fit.st)

#' 
#' ## Check 2: prior and posterior
#' 
#' After fitting, `jbplot_priorpost()` compares each prior with the posterior.
#' 
## ----priors-pp, fig.height = 5------------------------------------------------
jbplot_priorpost(fit)

#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | posterior equal to the prior | the data carry no information on that parameter | report it; sensitivity runs over the prior |
#' | posterior at the edge of the prior | prior-data conflict | revisit the prior; check the data and the model |
#' 
#' # Life history: r prior, generation time and Selex inputs
#' 
#' JABBA3 takes life-history information in the format of FLicc (length-based
#' assessment with `fiticc()`): `lhpar`, an `FLPar` with growth, natural
#' mortality, maturity at length and the weight-length relationship, and `selpar`,
#' the selectivity parameters on the length scale (e.g. from
#' `FLicc::selpars_flicc(fit)`). From these, the tools in this chapter give
#' 
#' - the maximum intrinsic rate of increase $r$ and the generation time GT from a
#'   Leslie matrix (`jb3_leslie()`);
#' - the r prior for the chosen production-function shape (`jb3_fmsy2r()`);
#' - guidance for the AR(1) persistence $\rho$ from GT (`jb3_rho_gt()`);
#' - the EB/SB curves for the Selex (spawning-biomass) model from length or age
#'   selectivity (`jb3_selex()`), with `jbplot_selex()`.
#' 
#' FLicc itself is not needed: its selectivity functions are included.
#' 
#' ## Leslie matrix: r and generation time
#' 
#' Life-history parameters of the FLicc alfonsino example (von Bertalanffy
#' $L_\infty$ and $k$, natural mortality $M$, length at 50% maturity, weight at
#' length $W = aL^b$; steepness `s`, default 0.7):
#' 
## ----lh-lhpar-----------------------------------------------------------------
FLPar <- FLCore::FLPar; FLPars <- FLCore::FLPars
lhpar <- FLPar(linf = 55.7, k = 0.08, M = 0.162, L50 = 31.1859,
               a = 0.004721956, b = 3.146168, s = 0.75)
jb3_leslie(lhpar)

#' 
#' $r = \log\lambda$, with $\lambda$ the dominant eigenvalue of the Leslie matrix
#' whose first row is the reproduction at the slope at the origin of the
#' Beverton-Holt stock-recruitment relationship, $4h/(\phi_0(1-h))$ times weight
#' and maturity at age; GT is the mean age of the spawners in the unfished
#' population. Biology at age follows from von Bertalanffy growth (`t0`, default
#' -0.5), $W = aL^b$ and logistic maturity at length ($L_{50}$, $L_{95}$; default
#' $L_{50} + 0.2 L_\infty$), with ages up to where 1% of a cohort survives
#' (`amax`). The calculation reproduces `JABBA::jbleslie()` and
#' `FLSRTMB::productivity()`.
#' 
#' Sex-structured parameters go in a second dimension of the `FLPar` (or
#' `list(F = , M = )`); reproduction is set by `sex` (default females):
#' 
## ----lh-sex-------------------------------------------------------------------
lh2 <- FLPar(array(c(55.7, 0.08, 0.162, 31.2, 0.0047, 3.15, 0.75,
                     50.0, 0.09, 0.180, 29.0, 0.0047, 3.15, 0.75), dim = c(7, 2, 1),
                   dimnames = list(params = c("linf", "k", "M", "L50", "a", "b", "s"),
                                   sex = c("F", "M"), iter = 1)))
lh2[, "F"]
rbind(F = jb3_leslie(lh2, sex = "F"), M = jb3_leslie(lh2, sex = "M"))

#' 
#' `jb3_leslie()` also takes an `FLStock` (median `m`, `mat` and `stock.wt` at age,
#' females from the `unit` dimension if present), `FLQuants` (`m`, `mat`, `wt`) by
#' age with `unit = c("F", "M")`, or by length together with `lhpar` for growth:
#' 
## ----lh-flq-------------------------------------------------------------------
data(ple4, package = "FLCore")
jb3_leslie(ple4, h = 0.7)[c("r", "gt")]

#' 
#' ## The r prior for the production-function shape
#' 
#' The Leslie $r$ is the per-capita growth at low biomass. In the Pella-Tomlinson
#' form of JABBA3 this equals $r/(m-1)$, while $F_{MSY} = r/m$: using the Leslie
#' $r$ directly is only consistent for the Schaefer model. Matching the prior on
#' $F_{MSY}$ (from an age-structured equilibrium, or $r_{Leslie}/2$ as the
#' Schaefer equivalent) keeps it consistent across shapes:
#' 
## ----lh-fmsy2r----------------------------------------------------------------
r.les <- jb3_leslie(lhpar)[["r"]]
Fmsy  <- r.les/2                                   # Schaefer equivalent
c(schaefer = jb3_fmsy2r(Fmsy, shape = "schaefer"),
  pella0.4 = jb3_fmsy2r(Fmsy, BmsyK = 0.4),
  pella0.36 = jb3_fmsy2r(Fmsy, BmsyK = 0.36),
  fox = jb3_fmsy2r(Fmsy, shape = "fox"))
jb3_r2fmsy(0.2, BmsyK = 0.4)                       # and back

#' 
#' ## Persistence from generation time
#' 
#' The generation time from the same Leslie matrix gives the guidance for the AR(1)
#' coefficient (chapter *Model options*):
#' 
## ----lh-rho-------------------------------------------------------------------
gt <- jb3_leslie(lhpar)[["gt"]]
c(GT = gt, rho = jb3_rho_gt(gt))

#' 
#' ## Selex inputs from length selectivity
#' 
#' With the Selex model (`jb3_model(ebsb = )`) the state is spawning biomass and the
#' indices and removals are mapped through the equilibrium ratio of exploitable
#' to spawning biomass, EB/SB, as a function of SB/SB0. `jb3_selex()` computes this
#' curve for each fleet from `lhpar` and length-based selectivity, over a grid of F
#' with Beverton-Holt recruitment, and fits the saturating curve
#' $c(u_1, u_2, u_3, P_1, P_2)$. Selectivity parameters use the FLicc names and
#' functions (logistic: `SL50`, `SL95`; `dsnormal`: `mode`, `lsd`, `rsd`; `normal`;
#' `binormal`; `bidsnormal`), here as an `FLPars` as returned by
#' `FLicc::selpars_flicc(fit)`; selectivity at length is converted to age by
#' integrating over the length-at-age distribution (`CVL = 0.1`).
#' 
## ----lh-selex-----------------------------------------------------------------
selpar <- FLPars(Gillnet = FLPar(mode = 32, lsd = 4, rsd = 15),
                 Longline = FLPar(SL50 = 28, SL95 = 36))
sx <- jb3_selex(lh2, selpar)
sx

#' 
## ----lh-selex-plot, fig.height = 4.5, fig.width = 10--------------------------
jbplot_selex(sx)

#' 
#' The JABBA 5-parameter selectivity (JABBA-Select; `jb3_sel5()`: `SL50`, `SL95`,
#' then a descending limb from `SLdesc` with coefficient of variation `CVdesc`
#' towards a minimum `Smin`) can be used at length, alongside the FLicc functions:
#' 
## ----lh-sel5-len, fig.height = 4.5, fig.width = 10----------------------------
sx5 <- jb3_selex(lhpar, list(Logistic = c(SL50 = 28, SL95 = 36),
                             Dome5 = c(SL50 = 25, SL95 = 32, SLdesc = 42, CVdesc = 0.15, Smin = 0.1)))
jbplot_selex(sx5)

#' 
#' or at age (`sel.type = "age"`; parameters in years):
#' 
## ----lh-sel5-age, fig.height = 4.5, fig.width = 10----------------------------
sxa <- jb3_selex(lhpar, list(Age5 = c(SL50 = 4, SL95 = 6, SLdesc = 15, CVdesc = 0.2, Smin = 0.2)),
                 sel.type = "age")
jbplot_selex(sxa, at = "age")

#' 
#' ## A Selex fit
#' 
#' The output of `jb3_selex()` goes straight into `jb3_model(ebsb = )`; with
#' several fleets, `ebsb.index` and `ebsb.catch` assign the curves to the indices
#' and the removals. For bigeye tuna, with illustrative life-history and longline
#' selectivity parameters:
#' 
## ----lh-selex-fit, fig.height = 4.5, fig.width = 10---------------------------
lh.bet <- FLPar(linf = 190, k = 0.2, M = 0.35, L50 = 105, a = 2.4e-5, b = 2.98, s = 0.8)
sx.bet <- jb3_selex(lh.bet, list(LL = c(SL50 = 80, SL95 = 110)))
c(jb3_leslie(lh.bet)[c("r", "gt")], rho = jb3_rho_gt(jb3_leslie(lh.bet)[["gt"]]))
fsx <- jabba3(d, pr, jb3_model(shape = "fox", ebsb = sx.bet), osa = FALSE)
jbplot_ebsb_status(fsx)

#' 
#' > Check that the selectivity is plausible against the length compositions
#' > (e.g. with an FLicc fit) before using the curve: it determines how the
#' > exploitable-biomass indices map onto the spawning-biomass state.
#' 
#' # Model options
#' 
#' The model structure is set with `jb3_model()`, separately from the priors.
#' 
#' | Setting | Default | Role |
#' |---|---|---|
#' | `shape` | `"pella"` | `"pella"` (shape from `BmsyK`), `"schaefer"` ($B_{MSY}/K = 0.5$) or `"fox"` ($m = 1.001$) |
#' | `BmsyK` | 0.4 | $B_{MSY}/K$ for the Pella-Tomlinson form |
#' | `lag` | 0 | production lag in years (age at maturity) |
#' | `rho` | 0 | AR(1) coefficient of the process deviations |
#' | `rw`, `rw.scale` | `FALSE`, 1 | random walk in the process deviations (the conditional $\rho = 1$ limit) |
#' | `sdb.form` | `"const"` | `"dyn"`: process error SD $\sigma_t = \sigma_{MSY}(B_t/B_{MSY})^{-k}$ |
#' | `sdb.fixed` | `NULL` | fix the process error SD (used by Tier 2) |
#' | `xlim` | `c(0.2, 3)` | bounds of $B/B_{MSY}$ in the dynamic SD |
#' | `biascor` | `FALSE` | lognormal bias correction of the process error |
#' | `start` | `"all"` | `"index"`: no process error before the first index year |
#' | `obs.error` | `"additive"` | `"additive"`: SD $= \sqrt{se^2 + sdi^2}$ (JABBA); `"spict"`: SD $= sdi \times w_t$ |
#' | `se.min`, `se.na` | `NULL`, 0.3 | minimum SE; SE where an index has SEs but a year has none |
#' | `sets.q`, `sets.var` | one per index | catchability and variance groups |
#' | `catch.error`, `sdc` | `TRUE`, 0.1 | catch with lognormal error |
#' | `ebsb`, `ebsb.index`, `ebsb.catch`, `ebsb.removal` | `NULL` | Selex: spawning biomass as the state, with EB/SB curves |
#' | `q.uniform` | `FALSE` | JABBA (JAGS) uniform prior on q |
#' | `P.bound`, `K.bound`, `sdi.bound` | JABBA | soft bounds |
#' 
#' ## Shape
#' 
## ----model-shape--------------------------------------------------------------
fits.shape <- lapply(c(schaefer = "schaefer", pella = "pella", fox = "fox"), function(s)
  jabba3(d, pr, jb3_model(shape = s), osa = FALSE))
t(sapply(fits.shape, function(f) c(Fmsy = f$refpts["Fmsy", "est"], MSY = f$refpts["MSY", "est"],
                                   BmsyK = f$refpts["BmsyK", "est"], B.Bmsy = tail(f$ts$est[f$ts$quantity == "B/Bmsy"], 1))))

#' 
#' > With the same r prior the three shapes imply different $F_{MSY}$; convert
#' > the r prior with the shape (section *Priors*) when comparing shapes.
#' 
#' ## Process error: lag and persistence (recommended)
#' 
#' The BioProcDyn simulations showed that the standard assumption of independent
#' process error with no lag makes surplus production models over-optimistic
#' about rebuilding, and more so for late-maturing stocks. The true process
#' deviations are persistent (lag-1 autocorrelation 0.7-0.8 for most life
#' histories) and a model with independent deviations absorbs this persistence
#' into r and K. The recommended configuration is (worked through on Atlantic blue
#' marlin in chapter *The rebuilding protocol*)
#' 
#' - a **production lag** at the age at 50% maturity, `lag = amat`;
#' - **AR(1) process deviations** with $\rho = 0.75$ (sensitivity 0.5 and 0.9);
#'   $\sigma_{MSY}$ stays the marginal SD, $d_t \sim N(\rho d_{t-1}, \sigma^2(1-\rho^2))$;
#' - the FishLife prior on $\sigma_{MSY}$ (default).
#' 
## ----model-ar1----------------------------------------------------------------
rec <- jabba3(d, pr, jb3_model(shape = "fox", lag = 3, rho = 0.75))
rec

#' 
## ----model-ar1-plot, fig.height = 6-------------------------------------------
jbplot_trj(list(Current = fit, `Lag + AR(1)` = rec), type = c("B/Bmsy", "F/Fmsy"))

#' 
#' ### Adjusting $\rho$ by generation time
#' 
#' $\rho = 0.75$ is the default. The persistence of process deviations increases
#' with generation time (GT): across age-structured simulations for 486 bony-fish
#' species of the FAO SOFIA list (FishLife life histories) the lag-1
#' autocorrelation followed $\rho \approx \exp(-2.2/GT)$ ($R^2 = 0.81$), and the
#' default corresponds to GT of about 6-10 years. For shorter- or longer-lived
#' stocks, adjust $\rho$ with the guidance table or `jb3_rho_gt()`; GT comes from
#' the same Leslie matrix as the r prior (e.g. `JABBA::jbleslie()` or
#' `FLSRTMB::productivity()`). Performance in simulation was flat between 0.5
#' and 0.9, so the adjustment matters mostly for short-lived species, where
#' persistence is weak.
#' 
## ----model-rho-table----------------------------------------------------------
knitr::kable(jb3_rho_table())
jb3_rho_gt(c(2, 4, 8, 15))

#' 
#' The scale prior on $\sigma_{MSY}$ stays tied to r (median
#' $0.80\,\sigma_R\,r^{0.45}$, unbiased across the 486 species), so it follows the
#' r prior directly in sensitivity runs on r.
#' 
#' $\rho$ can be estimated with a normal prior, by switching `priors$rho` on;
#' with a single index this is less stable than fixing it.
#' 
## ----model-rhoest-------------------------------------------------------------
pr.rho <- jb3_priors(d, r = c(0.2, 0.5), rho = c(0.7, 0.15, 1))
fr <- jabba3(d, pr.rho, jb3_model(shape = "fox", lag = 3), osa = FALSE)
fr$report$rho

#' 
#' ## Depletion-dependent process error
#' 
#' `sdb.form = "dyn"` makes the process error SD increase at low biomass,
#' $\sigma_t = \sigma_{MSY} x_t^{-k}$ with $x_t = B_t/B_{MSY}$ bounded to `xlim`,
#' and k estimated with the life-history prior `kdyn`. In BioProcDyn modelling
#' this inflation did not reduce the bias, so it is an option (and part of the
#' Tier 2 workflow, section *Workflows*), not the default.
#' 
## ----model-dyn----------------------------------------------------------------
fd <- jabba3(d, pr, jb3_model(shape = "fox", lag = 3, sdb.form = "dyn"), osa = FALSE)
fd$pars[c("sdb", "kdyn"), c("est", "lci", "uci")]

#' 
#' ## Index observation error
#' 
#' `obs.error = "additive"` (default, as JABBA) adds an estimated SD `sdi` to
#' the input standard errors, $SD_t = \sqrt{se_t^2 + sdi^2}$.
#' `obs.error = "spict"` uses SPiCT's form, $SD_t = sdi \times w_t$, where the
#' input SEs enter only as relative weights $w_t$ (as `stdevfacI`); a prior on
#' `sdi` is then available through `priors$sdi` (here on, as without it the SD of
#' the second bigeye index collapses towards zero). Indices can share catchability
#' (`sets.q`) or variance (`sets.var`), e.g. `sets.q = c(1, 1)` for two segments
#' of one fleet.
#' 
## ----model-obs----------------------------------------------------------------
pr.sdi <- jb3_priors(d, r = c(0.2, 0.5), sdi = c(0.1, 0.5, 1))
fsp <- jabba3(d, pr.sdi, jb3_model(shape = "fox", obs.error = "spict"), osa = FALSE)
fsp$pars[grep("sdi", rownames(fsp$pars)), c("est", "lci", "uci")]

#' 
#' ## Catch error
#' 
#' With `catch.error = TRUE` (default) the catch is estimated with lognormal error
#' (SD `sdc = 0.1`) as random effects, as JABBA's `add.catch.CV = TRUE`; this
#' absorbs small inconsistencies between catch and production. `catch.error = FALSE`
#' treats the catch as known.
#' 
#' ## Selex (spawning biomass as the state)
#' 
#' With `ebsb` the state is spawning biomass and the indices and removals are
#' mapped to it through the equilibrium ratio of exploitable to spawning biomass,
#' a function of depletion given as a curve `c(u1, u2, u3, P1, P2)` (or a named
#' list of curves, one per selectivity), assigned to indices with `ebsb.index`
#' and to the removals with `ebsb.catch`. The priors are then on the
#' spawning-biomass scale. The curves come from an age-structured equilibrium
#' analysis (e.g. FLBRP) or from a JABBA Selex setup via `jbinput2jabba3()`.
#' Selex is a special case for stocks where the fishery selects well below
#' maturity and is not covered further here.
#' 
#' `jbplot_ebsb_status()` shows the spawning biomass status (the state) with the
#' exploitable biomass status implied by each curve, EB/EB~MSY~ = (SB/SB~MSY~)
#' $\times$ $\rho$(SB/K)/$\rho$(SB~MSY~/K). With an illustrative curve:
#' 
## ----model-selex, fig.height = 4----------------------------------------------
fsel <- jabba3(d, pr, jb3_model(shape = "fox", ebsb = c(0.6, 1.2, 3, 0.1, 1)), osa = FALSE)
jbplot_ebsb_status(fsel)

#' 
#' # Fitting
#' 
#' `jabba3()` fits the model.
#' 
#' | Argument | Default | Role |
#' |---|---|---|
#' | `engine` | `"rtmb"` | `"rtmb"`: posterior mode, states Laplace-integrated; `"nuts"`: full posterior with tmbstan |
#' | `ndraws` | 0 | draws from the MVN approximation at the mode (Kobe plot, projections) |
#' | `nuts` | `list()` | tmbstan settings: `chains`, `iter`, `warmup`, `seed`, `cores`, `thin`, `control` |
#' | `osa` | `TRUE` | one-step-ahead residuals by particle filter |
#' | `se` | `TRUE` | standard errors by the delta method; `FALSE` for fast point estimates |
#' | `workflow`, `amat` | `"none"` | BioProcDyn workflows (section *Workflows*) |
#' | `run` | `NULL` | run name for plots |
#' | `peels`, `index.peels` | 0 | remove the last years of all data or of the indices (used by retro and hindcast) |
#' | `nstart` | 4 | guarded restarts |
#' 
#' ## Posterior mode
#' 
#' The default: the log-depletion states (and catch deviations) are integrated by
#' the Laplace approximation, and the fixed effects (log K, log r, log q, initial
#' depletion, process and observation SDs) are estimated at the mode of the
#' posterior. Intervals come from the delta method on the log scale.
#' 
## ----fit-mode-----------------------------------------------------------------
fit$time
fit$pars[, c("est", "lci", "uci")]
fit$refpts[, c("est", "lci", "uci")]
head(fit$ts[fit$ts$quantity == "B/Bmsy", ], 3)

#' 
#' The time series in `fit$ts` hold B, F, B/Bmsy, F/Fmsy, B/K, catch and the
#' process deviations with 95% intervals.
#' 
## ----fit-trj, fig.height = 6--------------------------------------------------
jbplot_trj(fit)

#' 
#' `jbplot_trj()` and `jbplot_advice()` share the JAGG options: `nested = TRUE`
#' for 90/60% bands (the default in `jbplot_advice()`), `extra` for further
#' reference points (e.g. $B_{pa}$; an `extra` $B_{lim}$ or $B_{thr}$ replaces the
#' default), `drop.last` to drop the last year of F and catch, `proj` to append
#' projections and `probs` for other bands. $B_{thr} = 0.5 B_{MSY}$ and
#' $B_{lim} = 0.3 B_{MSY}$ (SPiCT benchmark defaults) are drawn by default
#' (`bthr`, `blim`; `NULL` to omit). A list of fits is coloured by run.
#' 
## ----fit-advice-extra, fig.height = 7, fig.width = 7--------------------------
jbplot_advice(fit, extra = c(Bpa = 0.8), drop.last = TRUE)

#' 
#' `jbplot_kobe(fit, year = )` shows the uncertainty in any year (gold diamond),
#' and `jbplot_biplot()` is the 'post-modern' biplot of Quinn and Collie (2005),
#' as JABBA's: $B/B_{MSY}$ against $F/F_{MSY}$ with the zones of a hockey-stick
#' harvest control rule (target 0.8 $F_{MSY}$, threshold 0.2 $B_{MSY}$) and the
#' 50/80/95% density contours of the terminal status.
#' 
## ----fit-kobe-year, fig.height = 4.5------------------------------------------
jbplot_kobe(fit, year = 2000)

#' 
## ----fit-biplot, fig.height = 4.5---------------------------------------------
jbplot_biplot(fit)

#' 
#' ## Catch
#' 
#' With catch error (the default) the estimated catch can depart a little from the
#' observed catch. `jbplot_catch()` shows both, and `jbplot_catcherror()` adds the
#' log ratio of estimated to observed catch (as JABBA's `jbplot_catcherror()`).
#' 
## ----fit-catch, fig.height = 3.5----------------------------------------------
jbplot_catch(fit)

#' 
## ----fit-catcherror, fig.height = 4.5-----------------------------------------
jbplot_catcherror(fit)

#' 
#' ## Fast point estimates
#' 
#' For simulation testing, MSE and large grids, `se = FALSE, osa = FALSE` skips
#' the standard errors and residuals:
#' 
## ----fit-fast-----------------------------------------------------------------
system.time(fq <- jabba3(d, pr, jb3_model(shape = "fox"), se = FALSE, osa = FALSE))
fq$pars[c("K", "r"), "est"]

#' 
#' ## Full posterior with NUTS
#' 
#' `engine = "nuts"` samples all parameters and states of the same model with
#' tmbstan (no Laplace approximation), started from the mode. Posterior medians
#' and quantiles replace the mode summaries in `pars`, `refpts` and `ts`.
#' 
#' NUTS is meant for final runs, so the defaults are rigorous: 3 chains of 4000
#' iterations (half warmup). After sampling, `jabba3()` checks convergence and warns
#' if the rank-normalised $\hat{R}$ exceeds 1.01, the bulk or tail effective
#' sample size is below 100 per chain, there are divergent transitions, or E-BFMI
#' is below 0.3; the full table is in `fn$nuts$diags` (see *NUTS diagnostics*).
#' For bigeye the default run takes about half a minute.
#' 
## ----fit-nuts-----------------------------------------------------------------
fn <- jabba3(d, pr, jb3_model(shape = "fox"), engine = "nuts", nuts = list(cores = 3))
fn

#' 
#' ## Fit to the indices
#' 
## ----fit-fits, fig.height = 4-------------------------------------------------
jbplot_fits(fit)

#' 
#' `jbplot_fits(fit, log = TRUE)` plots on the log scale; the shaded band is the
#' posterior predictive interval including observation error.
#' 
#' ## Production function and process deviations
#' 
## ----fit-prod, fig.height = 4.5-----------------------------------------------
jbplot_spphase(fit)
jbplot_spdyn(fit)

#' 
#' `jbplot_spphase()` (as JABBA's `jbplot_spphase()`) plots the catch against
#' biomass over the production curve on the Kobe colours mapped to production
#' space: below the line $C = F_{MSY} B$ the stock is fished below $F_{MSY}$. The
#' hatched band left of $B_{MSY}$ marks catches below the production curve but
#' above that line. `jbplot_spdyn()` plots the realised surplus production,
#' $B_{t+1} - B_t + C_t$, against the expected curve; departures are the process
#' deviations (and the lag) in production units. `rel = TRUE` gives $B/K$ and
#' $C/MSY$ axes, and `jbplot_production(fit, type = "phase" | "dyn")` calls either.
#' 
## ----fit-prod-rel, fig.height = 4.5-------------------------------------------
jbplot_production(fit, type = "phase", rel = TRUE)

#' 
## ----fit-procdev, fig.height = 3.5--------------------------------------------
jbplot_procdev(fit)

#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | process deviations | long runs of the same sign | persistence: AR(1) with `rho = 0.75` and the lag |
#' | production plot | catch persistently above the production curve while biomass increases | the production function is wrong for the data: lag, shape, index trends |
#' 
#' ## Tables
#' 
#' `jb3_partable()` lists parameters, reference points and terminal status with 95%
#' intervals; with a list of fits and `format = "wide"` it compares runs side by side.
#' `jb3_tstable()` gives the time series of estimated quantities by year.
#' 
## ----tab-pars-----------------------------------------------------------------
knitr::kable(jb3_partable(fit), digits = 4)
knitr::kable(jb3_partable(list(Fox = fit, Lag.AR1 = rec), format = "wide"))
knitr::kable(tail(jb3_tstable(fit, c("Catch", "B/Bmsy", "F/Fmsy")), 5), digits = 3)

#' 
#' # Check 3: convergence
#' 
#' `jb3_convergence()` follows FLicc: the optimiser code, the largest gradient,
#' a positive definite Hessian, finite standard errors and the Newton decrement
#' (the expected further decrease of the objective by one Newton step). The fit
#' is classed `ok`, `warning` or `failed`.
#' 
## ----conv---------------------------------------------------------------------
jb3_convergence(fit)
fit$convergence$tier

#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | max gradient | > 0.001 | refit (`nstart`), check priors and bounds |
#' | Hessian | not positive definite | a parameter is not identified: priors on r or K, fewer q or variance groups |
#' | Newton decrement | > 0.1 | not at the optimum; simplify the model |
#' | estimates at bounds | K at `K.bound`, B/K at `P.bound` | prior-data conflict; revisit priors |
#' 
#' The fit internally polishes the optimum with Newton steps and restarts from
#' jittered values when the gradient is large (`nstart`).
#' 
#' # Check 4: residual diagnostics
#' 
#' ## One-step-ahead residuals
#' 
#' The residuals of a state-space model computed from the smoothed states are
#' shrunk and autocorrelated by construction, so the classic tests are biased.
#' JABBA3 computes **one-step-ahead (OSA) residuals** (Thygesen et al. 2017) by a
#' particle filter (`jb3_osa()`, run by `jabba3(osa = TRUE)`): each observation is
#' predicted from the data up to the previous one. Under a correct model they are
#' independent N(0, 1).
#' 
## ----diag-osa, fig.height = 4-------------------------------------------------
jbplot_residuals(fit)

#' 
#' `jbplot_residuals(fit, type = "ordinary")` shows the classic log residuals.
#' 
#' `jbplot_stdresiduals()` shows the ordinary residuals standardised by the total
#' observation SD, as JABBA's plot of the same name, with SDNR and its 95%
#' critical value:
#' 
## ----diag-std, fig.height = 4-------------------------------------------------
jbplot_stdresiduals(fit)

#' 
#' ## Posterior predictive check
#' 
#' `jbplot_ppc()` compares, draw by draw, the discrepancy of the observed indices
#' from the fitted ones, $D = \sum ((\log I - \log \hat{I})/\sigma)^2$, with that
#' of indices replicated from the observation model (as JABBA's `jbplot_PPC()`).
#' The Bayesian p-value $P(D^{rep} > D^{obs})$ is near 0.5 for a well-specified
#' observation model; values near 1 mean the indices are fitted more closely than
#' the observation error implies (here the process error absorbs part of the
#' index variation), values near 0 that they are fitted less well. It uses the
#' draws of the fit (NUTS or MVN) and generates MVN draws if there are none.
#' 
## ----diag-ppc, fig.height = 3.5-----------------------------------------------
pp <- jbplot_ppc(fit)
attr(pp, "ppc")
pp

#' 
#' ## Runs tests
#' 
#' The runs test (Wald-Wolfowitz) checks for too few runs of residual signs
#' (autocorrelation), with three-sigma limits from the average moving range (as
#' JABBA, ss3diags and JAGG). Green: passed ($p \ge 0.05$); red: failed.
#' 
## ----diag-runs, fig.height = 4------------------------------------------------
jbplot_runs(fit)

#' 
#' The test table is in the plot's attribute `"runs"`
#' (`attr(jbplot_runs(fit), "runs")`); `mixing = "two.sided"` also flags too many
#' runs. Likewise `jbplot_retro()` and `jbplot_hcxval()` return Mohn's rho and
#' MASE in attributes `"rho"` and `"mase"`.
#' 
#' ## The JABBA and SPiCT test sets
#' 
#' `jb3_diags()` gives one table per index: RMSE, SDNR, runs test and the share
#' outside the three-sigma limits (the JABBA set), and the bias t-test, Ljung-Box
#' test (lag 4) and Shapiro-Wilk normality test (the SPiCT set, as
#' `plotspict.diagnostic()`). `flag` is `TRUE` when any test has $p < 0.05$.
#' 
## ----diag-table---------------------------------------------------------------
jb3_diags(fit)

#' 
## ----diag-spict, fig.height = 5.5---------------------------------------------
jbplot_spdiags(fit)

#' 
#' > Residual tests check the error assumptions; they have little power to detect
#' > structural misspecification such as an omitted production lag. Use the
#' > retrospective and hindcast as well.
#' 
#' ## Retrospective analysis
#' 
#' `jb3_retro()` refits with the last 1 to `n` years of all data removed and
#' computes Mohn's rho.
#' 
## ----diag-retro, fig.height = 4-----------------------------------------------
rt <- jb3_retro(fit, n = 5)
rt$rho
jbplot_retro(rt)

#' 
#' ## Hindcast cross-validation
#' 
#' `jb3_hindcast()` removes the last 1 to `n` years of the indices (catch kept),
#' predicts each removed observation one step ahead, and computes the mean
#' absolute scaled error (MASE; Kell et al. 2021; Carvalho et al. 2021), the prediction error
#' relative to a naive prediction (the last observed value). MASE < 1: the model
#' predicts better than the naive guess.
#' 
## ----diag-hc, fig.height = 4--------------------------------------------------
hc <- jb3_hindcast(fit, n = 5)
hc$mase
jbplot_hcxval(hc)

#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | OSA residuals | trends, runs test failed | lag, AR(1), index conflict; separate `sdi` |
#' | Ljung-Box | p < 0.05 | autocorrelated residuals: process error persistence (AR(1)), or index trends in q |
#' | Mohn's rho | beyond -0.15 to 0.20 | retrospective pattern; check recent data and q |
#' | MASE | > 1 | no prediction skill for that index |
#' 
#' # NUTS diagnostics
#' 
#' > **Requirements.** NUTS runs through `tmbstan`, which needs `rstan` and a C++ toolchain (Rtools on Windows,
#' > Xcode command-line tools on macOS); see chapter *Installation*. Without them, `engine = "nuts"` stops with a
#' > message, and everything else in JABBA3 works.
#' 
#' A NUTS fit has to be checked for convergence of the chains before its
#' posterior is used. `jb3_nuts_diags()` gives the rank-normalised split
#' $\hat{R}$, bulk and tail effective sample sizes (ESS; Vehtari et al. 2021) for
#' every sampled parameter, including the states, and the sampler diagnostics by
#' chain: divergent transitions, hits of the maximum tree depth, mean acceptance,
#' step size and the energy Bayesian fraction of missing information (E-BFMI).
#' 
## ----nuts-diags---------------------------------------------------------------
nd <- jb3_nuts_diags(fn)
nd$sampler
head(nd$pars, 8)
nd$ok

#' 
#' `ok` is `TRUE` when $\hat{R} < 1.01$, bulk and tail ESS > 100 per chain, there
#' are no divergent transitions and E-BFMI > 0.3 in all chains.
#' 
#' The plots are on the sampled scale (mostly logs): `logK`, `logr`, `logq[i]`,
#' `lbk` (log initial depletion), `logsdb`, `logsdi[i]`, and the states `logP[t]`
#' and catch deviations `cdev[t]`. By default only the fixed effects are plotted;
#' `pars` selects parameters (`"logq"` matches all `logq[i]`).
#' 
## ----nuts-summary, fig.width = 11, fig.height = 9-----------------------------
jbplot_nutsdiags(fn)

#' 
#' `jbplot_nuts(fit, type)` gives the single plots:
#' 
#' | `type` | Shows | Healthy |
#' |---|---|---|
#' | `"trace"` | traces by chain; divergent transitions as red rug ticks | chains overlap, no trends |
#' | `"rank"` | rank histograms by chain | flat and overlapping |
#' | `"dens"` | densities by chain | the same in all chains |
#' | `"acf"` | autocorrelation by chain | drops to zero within a few lags |
#' | `"pairs"` | pairwise scatter with divergent transitions in red | no funnels or ridges; divergences not clustered |
#' | `"rhat"` | $\hat{R}$ and ESS for all parameters, states included | $\hat{R} < 1.01$, ESS > 100 per chain |
#' | `"energy"` | marginal energy against energy transitions | similar widths (E-BFMI > 0.3) |
#' 
## ----nuts-pairs, fig.height = 6.5, fig.width = 7------------------------------
jbplot_nuts(fn, "pairs", pars = c("logK", "logr", "logq", "logsdb"))

#' 
## ----nuts-trace, fig.height = 3-----------------------------------------------
jbplot_nuts(fn, "trace", pars = c("logK", "logr"), inc_warmup = TRUE)

#' 
#' Each type is also its own function: `jbplot_trace()`, `jbplot_rank()`,
#' `jbplot_mcmcdens()`, `jbplot_mcmcacf()`, `jbplot_pairs()`, `jbplot_rhat()` and
#' `jbplot_energy()`.
#' 
## ----nuts-single, fig.height = 3----------------------------------------------
jbplot_rank(fn, pars = c("logK", "logr", "logsdb"))
jbplot_mcmcdens(fn, pars = c("logK", "logr", "logsdb"))
jbplot_mcmcacf(fn, pars = c("logK", "logr", "logsdb"))
jbplot_trace(fn, pars = c("logK", "logr", "logsdb"))
jbplot_energy(fn)

#' 
## ----nuts-rhat, fig.height = 6------------------------------------------------
jbplot_rhat(fn)

#' 
## ----nuts-pairs2, fig.height = 4.5, fig.width = 5-----------------------------
jbplot_pairs(fn, pars = c("logK", "logr"))

#' 
#' | What | Red flag | Action |
#' |---|---|---|
#' | $\hat{R}$, ESS | $\hat{R} > 1.01$ or ESS < 100 per chain | more iterations (`iter`); check the pairs plot |
#' | divergent transitions | any | `control = list(adapt_delta = 0.99)`; look where they cluster in the pairs plot |
#' | tree depth | many hits of `max_treedepth` | strong correlations (K and r, K and q): tighter priors |
#' | E-BFMI | < 0.3 | poor exploration of the tails, often the process error SD and the states |
#' | K-r-q ridge | strong negative correlation of log K with log r and log q | expected with little contrast; informative r prior |
#' 
#' > In the short run above the states (`logP`) have low ESS: with process error
#' > SDs around 0.1 the states and `logsdb` form a funnel. Run longer for final
#' > results.
#' 
#' # BioProcDyn workflows
#' 
#' The workflows of the BioProcDyn protocol run automatically with
#' `jabba3(..., workflow = )`, or step by step with `jb3_tier1()` and
#' `jb3_tier2()` on any fit. All steps are kept in `$steps`.
#' 
#' | `workflow` | Steps |
#' |---|---|
#' | `"tier1"` | production lag at the age at maturity (`amat`) |
#' | `"tier2"` | step 1: constant process error (FishLife prior); step 2: dynamic $\sigma_t = \sigma_{MSY} x_t^{-k}$ with the level tuned on step 1, $\sigma_{MSY} = \hat\sigma_1 / \sqrt{\mathrm{mean}(x_{1,t}^{-2k})}$, k estimated with the life-history prior, plus the lag |
#' | `"selex"`, `"selex_tier1"`, `"selex_tier2"` | the same with the spawning-biomass state (needs `model$ebsb`) |
#' 
## ----wf-tier2-----------------------------------------------------------------
wf <- jabba3(d, pr, jb3_model(shape = "fox"), workflow = "tier2", amat = 3)
names(wf$steps)

#' 
## ----wf-steps, fig.height = 6-------------------------------------------------
jbplot_steps(wf)

#' 
#' Step by step, from an existing fit:
#' 
## ----wf-step------------------------------------------------------------------
t1 <- jb3_tier1(fit, amat = 3)
t1$model$lag

#' 
#' > The BioProcDyn simulations found that the depletion-dependent process error
#' > of Tier 2, although present in the operating models, did not remove the
#' > optimism; lag + AR(1) (`jb3_model(lag = amat, rho = 0.75)`) did. The tiers
#' > are kept for comparison and sensitivity.
#' 
#' # Projections
#' 
#' `jb3_project()` projects the stock under constant $F/F_{MSY}$ scenarios (`fmsy`,
#' the default: 0, 0.5, 0.75, 1 and 1.25 $F_{MSY}$) or constant catch (`catch`),
#' from the NUTS draws or from MVN draws at the mode (generated if the fit has
#' none), with process error (and its AR(1) persistence and the lag). F-based
#' scenarios apply $F = f \times F_{MSY}$ draw by draw, so they remain defined
#' when a draw is depleted; constant catches are capped at 95% of the biomass.
#' 
## ----proj---------------------------------------------------------------------
pj <- jb3_project(fit, fmsy = c(0, 0.5, 0.75, 1, 1.25), nyears = 10)
pj$scenarios

#' 
#' Process error in the projection years follows the fitted model. With the
#' default `procerr = "ar1"`, each draw continues its AR(1) process from its last
#' estimated deviation, $d_{t+1} = \rho d_t + \varepsilon_t$ with
#' $\varepsilon_t \sim N(0, \sigma^2(1-\rho^2))$, so the persistence carries into
#' the forecast and the marginal SD stays $\sigma$ (with $\rho = 0$ this is the iid
#' sampling of JABBA's `fw_jabba()`). `procerr = "mean"` uses only the conditional
#' mean $\rho^h d_T$ and `"none"` deterministic production. Projection tables
#' report the risk $P(B < B_{lim})$ with $B_{lim} = 0.3 B_{MSY}$ by default.
#' 
## ----proj-plot, fig.height = 6------------------------------------------------
jbplot_projection(pj)

#' 
#' ## Projection tables
#' 
#' `jb3_projtable()` gives a catch options table in the style of an ICES advice
#' sheet: $F/F_{MSY}$, $F$ and the catch in the first projection year, the median
#' $B/B_{MSY}$ in selected years with the change relative to the terminal year,
#' and the probabilities of $B \ge B_{MSY}$, $F \le F_{MSY}$ and of the Kobe
#' green quadrant.
#' 
## ----proj-tab-----------------------------------------------------------------
knitr::kable(jb3_projtable(pj, years = c(2018, 2022, 2027)))

#' 
#' `jb3_kobematrix()` is the Kobe II strategy matrix (ICCAT): the probability of
#' being in the green quadrant by scenario and year (`what = "B"` or `"F"` for
#' one reference point only).
#' 
## ----proj-k2sm----------------------------------------------------------------
knitr::kable(jb3_kobematrix(pj))

#' ## Kobe II strategy matrices
#' 
#' The Kobe II strategy matrix (K2SM) gives, by management option (rows) and
#' projection year (columns), the probability of meeting reference points.
#' `jb3_kobematrix()` computes it from any projection, in the ICCAT format by
#' default; `jb3_kobematrix_iotc()`, `jb3_kobematrix_gfcm()` and
#' `jb3_kobematrix_ices()` give the IOTC, GFCM and ICES formats, and
#' `jbplot_kobematrix()` plots either, shaded in grey as in the report tables
#' (`colours = "grey"`, default) or red to green (`colours = "kobe"`). All
#' scenarios of a projection share the same process-error draws, so intermediate
#' years are identical across rows and the probabilities change monotonically
#' with the TAC.
#' 
#' ### ICCAT: TAC options
#' 
#' The classic ICCAT matrices (a) $P(F < F_{MSY})$, (b) $P(B > B_{MSY})$ and (c)
#' both (the green Kobe quadrant), for constant catches after one intermediate
#' year at the status quo catch (`imp.years = 1`):
#' 
## ----k2sm-iccat---------------------------------------------------------------
pk2 <- jb3_project(fit, catch = seq(50000, 100000, 5000), nyears = 10, imp.years = 1)
k2 <- jb3_kobematrix(pk2, what = "all", percent = TRUE)
k2

#' 
## ----k2sm-iccat-plot, fig.height = 10, fig.width = 8--------------------------
jbplot_kobematrix(pk2)

#' 
#' `what = "green"` (default), `"B"` or `"F"` returns one matrix, and
#' `percent = FALSE` proportions. Projections from several models (e.g.
#' alternative production functions, or runs with and without the lag) are
#' combined into one matrix with equal or given `weights`, as the ICCAT combined
#' matrices of two base cases:
#' 
## ----k2sm-comb, eval = FALSE--------------------------------------------------
## jb3_kobematrix(list(Fox = pk2.fox, Schaefer = pk2.sch), what = "all", percent = TRUE)

#' 
#' ### IOTC: catch relative to the current catch
#' 
#' The IOTC format reports, for catch options relative to a reference catch
#' (default the catch in the last assessment year), the probability (%) of
#' violating the MSY-based reference points, $P(B < B_{MSY})$ and
#' $P(F > F_{MSY})$, and the limit reference points, $P(B < B_{lim})$ and
#' $P(F > F_{lim})$ (defaults $B_{lim} = 0.5 B_{MSY}$, $F_{lim} = 1.3 F_{MSY}$;
#' `blim`, `flim`), 3 and 10 years after the last assessment year (`years`):
#' 
## ----k2sm-iotc----------------------------------------------------------------
cref <- tail(fit$data$catch, 1)
pio <- jb3_project(fit, catch = cref*seq(0.6, 1.2, 0.1), nyears = 10)
kio <- jb3_kobematrix_iotc(pio, ref.catch = cref)
kio

#' 
## ----k2sm-iotc-plot, fig.height = 4.5, fig.width = 9--------------------------
jbplot_kobematrix(kio, colours = "kobe")

#' 
#' ### F options for GFCM and ICES
#' 
#' GFCM and ICES advice is usually based on F options. Named scenarios
#' (`scenarios` in `jb3_project()`) are reported in the order and with the names
#' given; each is one named value, with `Catch`, `Fmsy` (F/Fmsy), `F` (harvest
#' rate), `Fsq` or `Csq` (multipliers of the status quo).
#' 
#' F options relative to $F_{MSY}$ (`Fmsy = 0.8`) apply $F = 0.8 F_{MSY}$ draw by
#' draw, so $P(F < F_{MSY})$ is 0 or 100%. As in the GFCM and ICES forecasts
#' (e.g. SPiCT `add.man.scenario(fabs = )`), F options are therefore set as
#' absolute harvest rates from the point estimates of $F_{MSY}$ and the current
#' F (mean of the last three years), which gives the probability that each option
#' is below the uncertain $F_{MSY}$. The examples use Atlantic blue marlin
#' (chapter *The rebuilding protocol*), which is below $B_{MSY}$ and overfished,
#' with one intermediate year at the status quo catch:
#' 
## ----k2sm-fopt-fit------------------------------------------------------------
bum <- iccat$bum
d.bum  <- jb3_data(bum$catch, bum$cpue, bum$se, stock = "bum")
pr.bum <- jb3_priors(d.bum, r = c(0.098, 0.18), K = c(50000, sqrt(log(1 + 2^2))),
                     bkinit = c(1, sqrt(log(1 + 0.25^2))))
f.bum <- jabba3(d.bum, pr.bum, jb3_model(shape = "pella", BmsyK = 0.36, rho = 0.75, start = "index"),
                ndraws = 1000)

rp <- f.bum$refpts
Fmsy <- rp["Fmsy", "est"]
Fcur <- mean(tail(f.bum$ts$est[f.bum$ts$quantity == "F"], 3))

#' 
#' #### GFCM: Btgt-Ftgt, Bthr and Blim
#' 
#' `jb3_kobematrix_gfcm()` gives the three matrices of the GFCM advice summary
#' sheet: (a) $P(F < F_{tgt}$ and $B > B_{tgt})$, (b) $P(B > B_{thr})$ and (c) the
#' risk $P(B < B_{lim})$, with $B_{tgt} = B_{MSY}$, $F_{tgt} = F_{MSY}$,
#' $B_{thr} = 0.5 B_{MSY}$ and $B_{lim} = 0.3 B_{MSY}$ by default (`btgt`,
#' `ftgt`, `bthr`, `blim`). The F options follow the GFCM summary sheet: 80-120%
#' of $F_{cur}$, $F_{tgt}$ and the precautionary $F_{MSY}$ PA (35th percentile
#' of $F_{MSY}$):
#' 
## ----k2sm-gfcm----------------------------------------------------------------
Fmsy.pa <- exp(qnorm(0.35, log(Fmsy), rp["Fmsy", "se.log"]))
Ffrac <- c(seq(0.8, 0.95, 0.05), 1, 1.1, 1.2)
fgfcm <- c(setNames(lapply(Ffrac, function(m) c(F = m*Fcur)),
                    ifelse(Ffrac == 1, "Fcur", paste0(100*Ffrac, "% Fcur"))),
           list(Ftgt = c(F = Fmsy), "Fmsy PA" = c(F = Fmsy.pa)))
pg <- jb3_project(f.bum, scenarios = fgfcm, nyears = 10, imp.years = 1)
kg <- jb3_kobematrix_gfcm(pg)
kg

#' 
#' Probabilities are shaded darker (greener) the higher they are, risks darker
#' (redder) the higher the risk, on finer classes (`risk.breaks`):
#' 
## ----k2sm-gfcm-plot, fig.height = 9, fig.width = 8----------------------------
jbplot_kobematrix(kg, colours = "kobe")

#' 
#' #### ICES: Bpa-Fmsy and Blim
#' 
#' For ICES-type advice, `jb3_kobematrix_ices()` gives (a)
#' $P(F < F_{MSY}$ and $B > B_{pa})$ and (b) the risk $P(B < B_{lim})$, to compare
#' with the ICES 5% threshold, with the SPM proxies $B_{pa} = 0.5 B_{MSY}$ and
#' $B_{lim} = 0.3 B_{MSY}$ (`bpa`, `blim`), e.g. for $F_{cur}$ and 0.7-1.3
#' $F_{MSY}$:
#' 
## ----k2sm-ices, eval = FALSE--------------------------------------------------
## mult <- seq(0.7, 1.3, 0.1)
## fices <- c(list(Fcur = c(F = Fcur)),
##            setNames(lapply(mult, function(m) c(F = m*Fmsy)), ifelse(mult == 1, "Fmsy", paste0(mult, "Fmsy"))))
## pices <- jb3_project(f.bum, scenarios = fices, nyears = 10, imp.years = 1)
## jb3_kobematrix_ices(pices)

#' 
#' Other combinations follow from `jb3_kobematrix()` with `bref` and `fref`
#' (multiples of $B_{MSY}$ and $F_{MSY}$, optionally named; `NULL` drops one) and
#' `risk = TRUE` for the complement, e.g.
#' `jb3_kobematrix(pg, what = "all", bref = c(Btrigger = 0.5))`.
#' 

#' ## Short-term advice with a hockey-stick HCR
#' 
#' Short-term advice under a harvest control rule (as in the GFCM Type 1
#' management procedure with JABBA's `fw_jabba_ctrl()` and `jb_hockeystick()`)
#' takes three steps: a forecast with intermediate years at the status quo catch
#' followed by $F_{tgt}$, the hockey-stick rule applied to $B/B_{MSY}$ in the first
#' implementation (TAC) year to obtain the advised $F_{adv}$, and the forecast of
#' the management options. The example continues with the single blue marlin
#' model, which is close to $B_{trigger}$.
#' 
#' Forecast controls (`jb3_fwd()`, as `fw_jabba_ctrl()`) are data frames with one
#' row per year and `quant` `"Catch"`, `"Fmsy"` (F/F~MSY~), `"F"`, `"Csq"`,
#' `"Fsq"` or `"HCR"`; here two intermediate years at the status quo catch (mean
#' of the last three years) and two TAC years:
#' 
## ----hcr-ctrl-----------------------------------------------------------------
y0  <- max(f.bum$years)
yrs <- y0 + 1:4
csq <- mean(tail(f.bum$data$catch, 3))
fsq <- mean(tail(f.bum$ts$est[f.bum$ts$quantity == "F/Fmsy"], 3))
ctrl <- function(v, q) data.frame(year = yrs, value = c(csq, csq, v, v), quant = c("Catch", "Catch", q, q))
ctrl(1, "Fmsy")

#' 
#' `jb3_hcr()` evaluates the rule ($F/F_{MSY} = F_{tgt}$ above $B_{trigger}$,
#' decreasing linearly below) on $B/B_{MSY}$ in the TAC year of the $F_{tgt}$
#' forecast, either at the median (`Fadv.det`) or, precautionary and the default,
#' as the mean over all draws of the rule applied to each draw (`Fadv.pa`):
#' 
## ----hcr-fadv-----------------------------------------------------------------
hcr <- jb3_hcr(jb3_fwd(f.bum, list(Ftgt = ctrl(1, "Fmsy"))), ftgt = 1, btrigger = 0.5, blim = 0.3)
hcr

#' 
## ----hcr-plot, fig.height = 4.5, fig.width = 6--------------------------------
jbplot_hcr(hcr)

#' 
#' The median $B/B_{MSY}$ is above $B_{trigger}$, so the deterministic rule gives
#' $F_{tgt}$, but the draws below $B_{trigger}$ reduce the precautionary $F_{adv}$.
#' The management options, with $F_{adv}$, $F_{tgt}$, $F_{sq}$, $C_{sq}$ and
#' $F = 0$ in the TAC years:
#' 
## ----hcr-stf, fig.height = 6, fig.width = 9-----------------------------------
fc <- list(Fadv = ctrl(hcr$Fadv, "Fmsy"), Ftgt = ctrl(1, "Fmsy"), Fsq = ctrl(fsq, "Fmsy"),
           Csq = ctrl(csq, "Catch"), F0 = ctrl(0, "Fmsy"))
stf <- jb3_fwd(f.bum, fc)
jbplot_projection(stf)

#' 
#' `jb3_fwdtable()` summarises the fishing opportunities in the GFCM format (as
#' `FLRef::fwd2type1()`): the catch and $F/F_{tgt}$ in the TAC year, the biomass
#' at the start of the following year, its change, and the probabilities of
#' being above $B_{lim}$, $B_{pa}$ and $B_{tgt}$:
#' 
## ----hcr-table----------------------------------------------------------------
knitr::kable(jb3_fwdtable(stf))

#' 
#' Alternatively, `quant = "HCR"` applies the rule draw by draw in every
#' projection year (`value` = $F_{tgt}$, rule settings in `hcr`), e.g. for a
#' medium-term projection under the rule:
#' 
## ----hcr-drawwise, eval = FALSE-----------------------------------------------
## hcr.ctrl <- data.frame(year = y0 + 1:10, value = c(csq, csq, rep(1, 8)),
##                        quant = c("Catch", "Catch", rep("HCR", 8)))
## pm <- jb3_fwd(f.bum, list(HCR = hcr.ctrl), hcr = list(btrigger = 0.5))
## jbplot_projection(pm)

#' 
#' For a weighted ensemble, pass a `jb3_ensemble()` instead of the single fit;
#' the steps are the same.
#' 
#' Constant catch scenarios remain available, e.g. for TAC options:
#' 
## ----proj-c, fig.height = 6---------------------------------------------------
pc <- jb3_project(rec, catch = c(60000, 70000, 80000), nyears = 10)
knitr::kable(jb3_projtable(pc)[, 1:9])

#' 
#' # Comparing runs and risk
#' 
#' ## Ensembles and sensitivity runs
#' 
#' `jbplot_ensemble()` compares a list of fits (scenarios, sensitivity runs, the
#' steps of a workflow), as JABBA's plot of the same name: trajectories coloured
#' by run and a Kobe plot with each run's terminal uncertainty and its median.
#' With `joint = TRUE` the runs are pooled with equal weights into one ensemble.
#' 
## ----cmp-runs-----------------------------------------------------------------
fsch <- jabba3(d, jb3_priors(d, r = c(0.2, 0.5)), jb3_model(shape = "schaefer"), osa = FALSE)
runs <- list(Fox = fit, Schaefer = fsch, `Lag + AR(1)` = rec)

#' 
## ----cmp-ensemble, fig.width = 10, fig.height = 5-----------------------------
jbplot_ensemble(runs)

#' 
## ----cmp-ensemble-joint, fig.width = 10, fig.height = 5-----------------------
jbplot_ensemble(runs, joint = TRUE)

#' 
#' `jbplot_summary()` gives the overview of JABBA's `jbplot_summary()`: B, F,
#' $B/B_{MSY}$, $F/F_{MSY}$ and $B/K$ by run, and the production curves with
#' $B_{MSY}$ and MSY.
#' 
## ----cmp-summary, fig.width = 10, fig.height = 8------------------------------
jbplot_summary(runs)

#' 
#' `kobe = FALSE` drops the Kobe panel in `jbplot_summary()` and `jbplot_ensemble()`.
#' `jabbaplot3()` takes a list of runs too, with a choice of panels:
#' 
## ----cmp-jabbaplot3, fig.width = 10, fig.height = 6.5-------------------------
jabbaplot3(runs, panels = c("B/Bmsy", "F/Fmsy", "Catch", "kobe", "production", "procdev"))

#' 
#' ## Projections on the advice plot
#' 
#' `jbplot_advice(fit, proj = )` appends projection scenarios (dashed) to the
#' assessment:
#' 
## ----cmp-advice-proj, fig.height = 7, fig.width = 7---------------------------
jbplot_advice(fit, proj = pj)

#' 
#' ## Risk curves
#' 
#' `jbplot_crisk()` plots, as JABBA's `jbplot_Crisk()`, the probability that
#' $B < B_{MSY}$, that $F > F_{MSY}$ and that $B < 0.5 B_{MSY}$ (`bfrac`) against
#' the scenario value (F/F~MSY~ or catch), one line per projection year:
#' 
## ----cmp-crisk, fig.height = 4------------------------------------------------
jbplot_crisk(pj)

#' 
## ----cmp-crisk-catch, fig.height = 4------------------------------------------
jbplot_crisk(pc)

#' 
#' # The rebuilding protocol: Atlantic blue marlin
#' 
#' Bigeye tuna is the standard example of this guide. This chapter is the more
#' complex one: it introduces the BioProcDyn rebuilding protocol on a stock where
#' it matters. The 2018 ICCAT assessment of Atlantic blue marlin (*Makaira nigricans*) used JABBA
#' and Stock Synthesis (SS3), and the two differed markedly in how fast the stock was
#' projected to rebuild. This chapter reproduces the 2018 JABBA base case, applies
#' the BioProcDyn protocol step by step (life-history process error prior,
#' production lag, persistent AR(1) process error) and compares status, diagnostics
#' and forecasts.
#' 
#' ## The 2018 base case
#' 
#' The 2018 base case was a Pella-Tomlinson model with $B_{MSY}/K = 0.36$, an r
#' prior of 0.098 (CV 0.18), an uninformative K prior, initial depletion near
#' unfished, ten CPUE indices with an additional observation error of 0.25, and the
#' process error fixed at $\sigma = 0.07$ and process deviations starting with the
#' first index year. In JABBA3 (the additional observation variance is estimated
#' by index):
#' 
## ----bum-build----------------------------------------------------------------
bum <- iccat$bum
d.bum  <- jb3_data(bum$catch, bum$cpue, bum$se, stock = "bum")
pr.bum <- jb3_priors(d.bum, r = c(0.098, 0.18), K = c(50000, sqrt(log(1 + 2^2))),
                     bkinit = c(1, sqrt(log(1 + 0.25^2))))
m.base <- jb3_model(shape = "pella", BmsyK = 0.36, sdb.fixed = 0.07, start = "index")
bum.base <- jabba3(d.bum, pr.bum, m.base, run = "2018 base")
bum.base

#' 
#' The same fit follows from the original JABBA input with the converter
#' (`jbinput2jabba3()`, chapter *Converting from JABBA and SPiCT*):
#' 
## ----bum-convert, eval = FALSE------------------------------------------------
## jb.bum <- JABBA::build_jabba(catch = bum$catch, cpue = bum$cpue, se = bum$se,
##   assessment = "BUM", scenario = "2018 base", model.type = "Pella", BmsyK = 0.36,
##   K.prior = c(50000, 2), r.prior = c(0.098, 0.18), psi.prior = c(1, 0.25),
##   fixed.obsE = 0.25, proc.dev.all = FALSE, sigma.proc = 0.07, verbose = FALSE)
## x.bum <- jbinput2jabba3(jb.bum, keep.igamma = TRUE)
## bum.base <- jabba3(x.bum$data, x.bum$priors, x.bum$model, run = "2018 base")

#' 
## ----bum-base-plot, fig.height = 9, fig.width = 10, eval = TRUE---------------
plot(bum.base)

#' 
#' The fit is close to the 2018 JAGS assessment. Note the process deviations:
#' with $\sigma$ fixed at a low value they form long runs of the same sign, a
#' lag-1 autocorrelation of
#' 
## ----bum-acf, eval = TRUE-----------------------------------------------------
pd <- bum.base$ts$est[bum.base$ts$quantity == "Process deviation"]
round(acf(pd, plot = FALSE)$acf[2], 2)

#' 
#' This is the persistence that the BioProcDyn simulations identify as the source
#' of rebuilding optimism: a model with independent process error absorbs the
#' persistent part of the production anomalies into r and K.
#' 
#' ## Joint versus separate indices
#' 
#' Chapter *Data exploration* combined the ten indices into one joint index
#' (`bum.joint`) and found some conflicts among them. Does it matter for the
#' assessment whether the production model is fitted to the ten separate indices
#' or to the joint index? Same priors and settings, the joint index with its own
#' log-scale SE:
#' 
## ----bum-joint, eval = TRUE---------------------------------------------------
if(!exists("bum.joint")) bum.joint <- jb3_jointindex(jb3_trend(bum$cpue, bum$se))
d.joint  <- jb3_data(bum$catch, bum.joint, stock = "bum")
pr.joint <- jb3_priors(d.joint, r = c(0.098, 0.18), K = c(50000, sqrt(log(1 + 2^2))),
                       bkinit = c(1, sqrt(log(1 + 0.25^2))))
bum.ji <- jabba3(d.joint, pr.joint, m.base, run = "joint index")
bum.js <- list(`10 indices` = bum.base, `joint index` = bum.ji)
pt <- jb3_partable(bum.js, format = "wide")
knitr::kable(pt[pt$Parameter %in% c("K", "r", "MSY", "B/Bmsy", "F/Fmsy"), ], row.names = FALSE)

#' 
## ----bum-joint-plot, fig.height = 5, fig.width = 10, eval = TRUE--------------
jabbaplot3(bum.js, panels = c("B/Bmsy", "F/Fmsy", "kobe", "spdyn"))

#' 
## ----bum-joint-fit, fig.height = 3.5, fig.width = 9, eval = TRUE--------------
patchwork::wrap_plots(jbplot_fits(bum.ji), jbplot_runs(bum.ji), ncol = 2)

#' 
## ----bum-joint-diag, eval = TRUE----------------------------------------------
dg <- lapply(bum.js, jb3_diags)
knitr::kable(do.call(rbind, lapply(names(dg), function(n) data.frame(run = n, dg[[n]][, c("index", "n", "RMSE", "SDNR", "runs.p", "flag")]))),
             digits = 3, row.names = FALSE, caption = "Residual diagnostics (OSA residuals) by index")

#' 
#' Status, reference points and productivity are practically the same: $B/B_{MSY}$
#' in 2016 of about 0.78 and $F/F_{MSY}$ of 0.80 with either input, and the joint
#' index gives a slightly narrower interval for $B/B_{MSY}$. The ten indices carry
#' one common signal, and the production model extracts the same signal from the
#' joint index. The separate fit flags six of the ten indices on their residual diagnostics;
#' the joint index has a small RMSE (0.20) but its OSA residuals still fail the
#' runs test (p = 0.02), which reflects the low fixed process error of the base case
#' rather than conflict among indices.
#' 
#' The rest of this chapter continues with the ten separate indices, as in the
#' 2018 assessment: each index keeps its own catchability and observation
#' variance, and the diagnostics can still point to individual series. The joint
#' index is a useful check on the combined signal and a simpler alternative where
#' the separate indices are numerous and short.
#' 
#' ## Applying the protocol
#' 
#' Three steps, each added to the previous one:
#' 
#' 1. **Life-history prior on the process error scale** (the JABBA3 default,
#'    $\sigma_{MSY} = 0.80 \sigma_R r^{0.45}$) instead of the fixed $\sigma = 0.07$;
#' 2. **production lag** at the age at maturity (here 3 years);
#' 3. **AR(1) process deviations** with $\rho = 0.75$ (the default; blue marlin falls
#'    in the 6-10 year generation-time class of `jb3_rho_table()`).
#' 
## ----bum-protocol, eval = TRUE------------------------------------------------
# pr.bum already carries the life-history sdb prior (default); the base case fixed sdb in the model
mod <- function(lag = 0, rho = 0) jb3_model(shape = "pella", BmsyK = 0.36, start = "index", lag = lag, rho = rho)
bum.runs <- list(
  `2018 base`     = bum.base,
  `sdb prior`     = jabba3(d.bum, pr.bum, mod(),        run = "sdb prior"),
  `+ lag`         = jabba3(d.bum, pr.bum, mod(3),       run = "+ lag"),
  `+ lag + AR(1)` = jabba3(d.bum, pr.bum, mod(3, 0.75), run = "+ lag + AR(1)"))
pt <- jb3_partable(bum.runs, format = "wide")
knitr::kable(pt[pt$Parameter %in% c("K", "r", "Bmsy", "Fmsy", "MSY", "sigma.proc", "rho", "B/Bmsy", "F/Fmsy"), ], row.names = FALSE)

#' 
## ----bum-compare, fig.height = 7.5, fig.width = 10, eval = TRUE---------------
jabbaplot3(bum.runs, panels = c("B/Bmsy", "F/Fmsy", "Catch", "kobe", "spdyn", "procdev"))

#' 
#' Each step lowers the estimated status. The base case had $B/B_{MSY} \approx 0.8$
#' in 2016; with the life-history prior, the lag and persistent process error the
#' estimate is about 0.45, with $F > F_{MSY}$.
#' 
#' ## Diagnostics
#' 
#' Retrospective patterns and index prediction skill do not distinguish the runs;
#' both are acceptable in all of them. As in the simulations, the standard
#' diagnostics cannot tell whether a model is optimistic about rebuilding.
#' 
## ----bum-diags, eval = TRUE---------------------------------------------------
bum.hc <- lapply(bum.runs, jb3_hindcast, n = 5)
bum.rt <- lapply(bum.runs, jb3_retro, n = 5)
mase <- sapply(bum.hc, function(h){ m <- h$mase; stats::setNames(m$MASE, m$index)[is.finite(m$MASE)] })
mohn <- sapply(bum.rt, function(r) stats::setNames(r$rho$rho, r$rho$quantity))
knitr::kable(round(rbind(mase, mohn[c("B/Bmsy", "F/Fmsy"), , drop = FALSE]), 2),
             caption = "Hindcast MASE by index and Mohn's rho")

#' 
#' Hindcasts for the two indices that extend to the end of the series and
#' dominate the recent trend (US and Brazilian recreational/longline CPUE), for the
#' 2018 base case and the protocol:
#' 
## ----bum-hcxval, fig.height = 7, fig.width = 9, eval = TRUE-------------------
patchwork::wrap_plots(
  jbplot_hcxval(bum.hc$`2018 base`, indices = c("USA.LL", "BRA.Rec"), ncol = 1) + ggplot2::ggtitle("2018 base"),
  jbplot_hcxval(bum.hc$`+ lag + AR(1)`, indices = c("USA.LL", "BRA.Rec"), ncol = 1) + ggplot2::ggtitle("Lag + AR(1) protocol"),
  ncol = 2)

#' 
#' ## Forecasts
#' 
#' The runs differ most in their forecasts. Constant-catch projections of 1500,
#' 2000 and 2500 t, with the AR(1) process continued from the last estimated
#' deviation in the persistent model:
#' 
## ----bum-proj, eval = TRUE----------------------------------------------------
bum.pc <- lapply(bum.runs[c("2018 base", "+ lag + AR(1)")], jb3_project,
                 catch = c(1500, 2000, 2500), nyears = 10)
lapply(bum.pc, function(p) knitr::kable(jb3_kobematrix(p)[, c(1, 2, 4, 6, 11)],
                                        caption = "P(green) by year"))

#' 
## ----bum-proj-plot, fig.height = 4.5, fig.width = 10, eval = TRUE-------------
patchwork::wrap_plots(
  jbplot_projection(bum.pc$`2018 base`, type = "B/Bmsy") + ggplot2::ggtitle("2018 base"),
  jbplot_projection(bum.pc$`+ lag + AR(1)`, type = "B/Bmsy") + ggplot2::ggtitle("Lag + AR(1) protocol"))

#' 
#' At 2000 t the base case projects a 70% probability of being in the green Kobe
#' quadrant by 2026, the protocol 28%: the base case expects the stock to rebuild
#' much faster. In simulation (BioProcDyn), forecasts of this kind from a model
#' with independent process error were optimistic by 25-45% in $B/B_{MSY}$ during
#' rebuilding, and almost two thirds of stocks still below $B_{MSY}$ were classified
#' as rebuilt, whereas the lag and persistent process error brought the truth
#' close to the median of the estimated distribution.
#' 
#' > The comparison with the 2018 SS3 assessment will be added here.
#' 
#' # Other data features
#' 
#' The data-exploration chapter looked at index trends, catchabilities and
#' conflicts among indices of one stock. This chapter covers other common data
#' situations, with two anonymised example datasets:
#' 
## ----od-data------------------------------------------------------------------
data(jara_examples)
names(jara_examples)
names(jara_examples$areas)
names(jara_examples$regions)

#' 
#' - `areas`: a catch history from six fleets or areas, two of them with missing
#'   early years; relative survey indices with log-SEs for four areas; and
#'   swept-area survey biomass for the same four areas;
#' - `regions`: indices from three regions, a survey each in two adjacent areas
#'   (`areaA`, `areaB`) and two surveys and five CPUE series in a third region
#'   (`C.*`).
#' 
#' ## Catch history with missing years
#' 
#' Two fleets report only in the later years (`fleet4`, `fleet5`):
#' 
## ----od-catch-data------------------------------------------------------------
cc <- jara_examples$areas$catch
colSums(!is.na(cc[, -1]))

#' 
#' A census trend model (one state per fleet, summed to a total) interpolates the
#' missing values. Where a fleet's early catches are not known at all, start values
#' anchor the interpolation: here `fleet4` is assumed to have caught 10 t in the
#' first year and `fleet5` the same as `fleet2`. The start values get a larger
#' log-SE (0.2) than the reported catches (0.05); `jb3_catchfill()` takes the SE as
#' a data frame with the same layout as the catch. A time block lets the trend
#' change after a given year, here when the catches peaked:
#' 
## ----od-catch, fig.height = 5-------------------------------------------------
cc[1, "fleet5"] <- cc[1, "fleet2"]; cc[1, "fleet4"] <- 10
sc <- cc; sc[, -1] <- 0.05; sc[1, c("fleet4", "fleet5")] <- 0.2
cf <- jb3_catchfill(cc, se = sc, timeblock = 2016)
jbplot_trend(attr(cf, "fit"))          # fleet fits and the total

#' 
#' The total catch with its year-specific log-SE, larger where fleets with assumed
#' start values are interpolated:
#' 
## ----od-catch-total, fig.height = 3.5-----------------------------------------
jbplot_catchfill(cf)
knitr::kable(as.data.frame(cf)[c(1:3, 15:18, 24:26), c("year", "catch", "se.log", "lci", "uci")],
             digits = 3, row.names = FALSE)

#' 
#' The filled table, with the original series, the missing values replaced by
#' predictions, and the total:
#' 
## ----od-catch-table-----------------------------------------------------------
knitr::kable(head(jb3_fill_table(cf, digits = 1), 6), row.names = FALSE)

#' 
#' The reconstructed catch goes into `jb3_data()` with its year-specific SE,
#' which the catch error of the model uses in place of the constant `sdc`; the
#' `jb3catchfill` object can be passed directly as the catch:
#' 
## ----od-catch-data2-----------------------------------------------------------
d.cf <- jb3_data(cf, jara_examples$areas$idx, jara_examples$areas$idx.se, stock = "areas")
round(range(d.cf$catch.se), 3)

#' 
#' ## Swept-area surveys in several areas
#' 
#' A survey covering four areas gives two kinds of index: relative indices by
#' area (`idx`, with log-SEs) and swept-area biomass estimates (`bio`). The
#' relative indices are combined into a joint index as in the data-exploration
#' chapter:
#' 
## ----od-joint, fig.height = 3.5-----------------------------------------------
a <- jara_examples$areas
ji.a <- jb3_trend(a$idx, a$idx.se, run = "relative")
jdx.a <- jb3_jointindex(ji.a)
jbplot_jointindex(ji.a)

#' 
#' The swept-area estimates are absolute and additive over areas, so a census
#' model fits each area separately and sums the area states to a total biomass
#' with uncertainty, smoothing through years with missing surveys:
#' 
## ----od-census, fig.height = 5------------------------------------------------
bi <- jb3_trend(a$bio, a$bio.se, type = "census", run = "swept-area census")
jbplot_trend(bi)

#' 
#' The census total and the relative joint index, both scaled to their means:
#' 
## ----od-census-cmp, fig.height = 3.5------------------------------------------
tb <- bi$total
sc1 <- mean(jdx.a$index, na.rm = TRUE); sc2 <- mean(tb$est)
cmp2 <- rbind(data.frame(year = jdx.a$year, index = jdx.a$index/sc1, lci = jdx.a$lci/sc1, uci = jdx.a$uci/sc1,
                         series = "relative joint index"),
              data.frame(year = tb$year, index = tb$est/sc2, lci = tb$lci/sc2, uci = tb$uci/sc2,
                         series = "census total biomass"))
round(cor(log(jdx.a$index), log(tb$est[match(jdx.a$year, tb$year)]), use = "complete.obs"), 2)
ggplot2::ggplot(cmp2, ggplot2::aes(year, index, colour = series, fill = series)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.15, colour = NA) + ggplot2::geom_line() +
  ggplot2::labs(y = "Index (scaled to mean)", colour = NULL, fill = NULL) + theme_jabba3()

#' 
#' The two agree closely. They answer different questions, however: the joint
#' index weights areas by the precision of their indices, the census total by
#' their biomass, so an area with much of the biomass but a noisy index counts
#' more in the census. Where the trends differ, the census total is the better
#' measure of the whole stock; it can also go into the production model as an index of
#' its own, with the log-SE of the census total.
#' 
#' ## Do indices among areas conflict? Stock structure
#' 
#' When indices come from several regions, the question is whether they follow a
#' common trend, in which case the regions can be assessed together, or diverge,
#' which may point to separate stocks or at least to regions with different
#' dynamics. The `regions` data have surveys of two adjacent areas and seven
#' indices of a third region C:
#' 
## ----od-reg-data--------------------------------------------------------------
rg <- jara_examples$regions
colSums(!is.na(rg[, -1]))

#' 
#' **Within regions.** First the joint index of each region: A and B together
#' (adjacent areas), and C:
#' 
## ----od-reg-fits--------------------------------------------------------------
tr.AB <- jb3_trend(rg[, c("year", "areaA", "areaB")], run = "A + B")
tr.C  <- jb3_trend(rg[, c("year", grep("^C\\.", names(rg), value = TRUE))], run = "C")
tr.C

#' 
## ----od-reg-cor-c, fig.height = 4.5, fig.width = 5.5--------------------------
jbplot_indexcor(tr.C)

#' 
#' The five CPUE series of region C are strongly correlated (r 0.77-0.94); the two
#' surveys are noisier and agree less with each other and with the CPUE
#' (the large `sdi` values). Within region C, the indices share one trend.
#' 
#' **Between regions.** The joint indices of the two regions on one plot, and
#' their correlation:
#' 
## ----od-reg-cmp, fig.height = 3.5---------------------------------------------
jAB <- jb3_jointindex(tr.AB, name = "A + B"); jC <- jb3_jointindex(tr.C, name = "C")
reg <- rbind(data.frame(jAB, region = "A + B"), data.frame(jC, region = "C"))
reg <- do.call(rbind, lapply(split(reg, reg$region), function(z){ m <- mean(z$index, na.rm = TRUE)
z$index <- z$index/m; z$lci <- z$lci/m; z$uci <- z$uci/m; z }))
ggplot2::ggplot(reg, ggplot2::aes(year, index, colour = region, fill = region)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.15, colour = NA) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
  ggplot2::labs(y = "Joint index (scaled to mean)", colour = NULL, fill = NULL) + theme_jabba3()
J <- merge(data.frame(year = jAB$year, AB = jAB$index), data.frame(year = jC$year, C = jC$index))
ok <- is.finite(J$AB) & is.finite(J$C)
cor.test(log(J$AB[ok]), log(J$C[ok]))[c("estimate", "p.value")]

#' 
#' **All indices in one model.** If the regions were one stock, a single trend
#' model with all nine indices should fit them without systematic residuals. The
#' pairwise correlations and residual diagnostics of the joint fit:
#' 
## ----od-reg-all, fig.height = 6, fig.width = 7--------------------------------
tr.all <- jb3_trend(rg, run = "all regions")
jbplot_indexcor(tr.all)

#' 
#' The additional observation SD each index needs, fitted within its region and
#' in the model with all regions:
#' 
## ----od-reg-sdi---------------------------------------------------------------
sdi <- function(f) stats::setNames(f$sdi, f$indices)
knitr::kable(data.frame(index = tr.all$indices, within.region = c(sdi(tr.AB), sdi(tr.C))[tr.all$indices],
                        all.regions = sdi(tr.all)), digits = 3, row.names = FALSE)

#' 
## ----od-reg-resid, fig.height = 4---------------------------------------------
jbplot_residuals(tr.all)

#' 
## ----od-reg-runs, fig.height = 6, fig.width = 9-------------------------------
p.reg <- jbplot_runs(tr.all)
p.reg
attr(p.reg, "runs")[, c("index", "p.runs", "test")]

#' 
#' The regional joint indices do not agree: the two regions are uncorrelated
#' (r = 0.06 over 29 shared years), with A + B lowest in the late 1980s and the
#' last years and C highest around 1998-1999 and stable or increasing since 2012.
#' Within regions the indices agree (A and B r = 0.35; the region C CPUE series
#' r 0.77-0.94). Forced onto one trend, the common state follows the numerous
#' and precise region C CPUE series, and the cost falls on the A and B surveys:
#' their additional SD rises from 0.08 to 0.29 (`areaA`) and from 0.23 to 0.39
#' (`areaB`), while the region C CPUE series fit almost as well as on their own,
#' and none of the A or B pairs with region C indices are significantly
#' correlated. The residuals of A and B run negative in the last years, where
#' region C is stable. The abundance data are thus not consistent with a single
#' unit: A and B can be combined, but region C should be assessed separately, or
#' a combined assessment checked against separate ones.
#' 
#' **How to use this.** Agreement among regions (positive correlations, no
#' residual patterns that separate the regions, similar trends in the regional
#' joint indices) supports assessing the regions together with one joint index or
#' with all indices in one model. Conflict between regions, while the indices
#' within each region agree, suggests that the regions have different dynamics:
#' assess them separately, or at least test the combined assessment against
#' separate ones. The trend models cannot by themselves establish stock
#' structure, which needs biological evidence (genetics, tagging, life history),
#' but they show whether the abundance data are consistent with a single unit.
#' 
#' # Special feature: catch-only assessment and auxiliary data
#' 
## ----co-time0, include = FALSE------------------------------------------------
co.t0 <- Sys.time()

#' 
#' Most stocks have no reliable abundance index. This chapter shows how JABBA3
#' handles the data-limited end of the range: catch-only assessments with priors
#' on stock status, auxiliary data such as fishing effort or an index of
#' exploitation, and how each kind of data changes the estimate of stock status.
#' 
#' > **The effort and F/Fmsy series used here are hypothetical.** Bigeye tuna has
#' > no effort or F/Fmsy series in the package. Both series below
#' > (`iccat$bet$effort_hyp` and `iccat$bet$ffmsy_hyp`) were *simulated* from the
#' > F trajectory of the bigeye example fit: effort with 2% annual effort creep
#' > and lognormal noise (SD 0.15), and F/Fmsy with noise (SD 0.3) for 1990-2009
#' > only. They are not real data. They are used only to show the methods, and
#' > their fits are partly circular, because they were simulated from an
#' > assessment that used the same catch.
#' 
#' ## Concept: hierarchical assessment framework
#' 
#' Kell et al. (2025) set out a *hierarchical assessment framework* (HAF) for
#' biomass dynamic models, in which the assessment method follows the data
#' available. There are four tiers:
#' 
#' | Tier | Data | Method |
#' |:--|:--|:--|
#' | Robust quantitative assessment | reliable abundance indices and catch | JABBA |
#' | Preliminary quantitative assessment | partial or short indices, catch | JABBA with scenarios and sensitivity runs |
#' | Auxiliary data | auxiliary indices (length-based F/Fmsy, effort), catch | auxiliary data as priors or as time series of exploitation |
#' | Data-limited | catch only | catch-only model with predictive priors (e.g. current depletion) |
#' 
#' The aim is *risk equivalence*: the risk of misclassifying stock status
#' against $B_{MSY}$ and $F_{MSY}$ should be the same in every tier. It is
#' measured by classification skill, for example the true skill statistic
#' (sensitivity + specificity - 1) for $B < B_{MSY}$ or $F > F_{MSY}$. Each tier
#' starts from *prior predictive modelling* (PPM): the model is run with the
#' priors and the catch but without the data, as a reality check of what the
#' assumptions alone imply. Kell et al. found that auxiliary F/Fmsy indices,
#' especially from periods of high exploitation, substantially improve the
#' classification of status.
#' 
#' JABBA3 covers all tiers with one model:
#' 
#' - `jb3_data(catch, indices = NULL)` gives a catch-only data set;
#' - `jb3_priors(status = jb3_bprior(...))` adds priors on $B/B_{MSY}$, $B/K$ or
#'   $F/F_{MSY}$ in any year, as JABBA's `b.prior`;
#' - `jb3_data(auxiliary = , auxiliary.type = )` adds auxiliary series of type
#'   `"effort"`, `"f"`, `"ffmsy"`, `"bbmsy"`, `"bk"` or `"z"`;
#' - `jabba3(..., prior.only = TRUE)` switches off the likelihood of the indices
#'   and auxiliary data (the catch still drives the dynamics), which gives the
#'   prior predictive run.
#' 
#' Without an index, the posterior mode and its Laplace approximation are
#' unreliable: there is nothing to pin down the states. `jabba3()` therefore
#' samples with NUTS by default when there is no index or when
#' `prior.only = TRUE` (`engine = "rtmb"` can still be set explicitly).
#' 
#' The runs below move up the tiers for bigeye tuna:
#' 
#' 1. catch only (tier 4): the prior predictive run, then the same with a
#'    terminal-year $B/B_{MSY}$ prior;
#' 2. catch + auxiliary data (tier 3): hypothetical effort with 2% creep (and a
#'    no-creep sensitivity), and a hypothetical F/Fmsy index from a period of
#'    high exploitation;
#' 3. catch + index (tier 1): the two Japanese longline CPUE series;
#' 4. catch + auxiliary data + index.
#' 
#' All runs use the Fox model with the r prior of the bigeye example, the default
#' K and initial depletion priors, and NUTS with three chains of 4000 iterations
#' (2000 warm-up; the `jabba3()` default).
#' 
## ----co-setup-----------------------------------------------------------------
bet <- iccat$bet
pr.r <- c(0.2, 0.5)
mod  <- jb3_model(shape = "fox")
nu   <- list(chains = 3, cores = 3)   # 4000 iterations per chain (default), 2000 warm-up

#' 
#' ## Catch only: prior predictive and status priors
#' 
#' With catch only, the information comes from the priors on r, K and initial
#' depletion, and from the catch itself, which has to be sustainable under the
#' estimated production function. The prior predictive run shows what these
#' assumptions imply:
#' 
## ----co-ppm-------------------------------------------------------------------
d0  <- jb3_data(bet$catch, NULL, stock = "bet")
d0
pr0 <- jb3_priors(d0, r = pr.r)
co.ppm <- jabba3(d0, pr0, mod, prior.only = TRUE, nuts = nu, run = "Prior predictive")

#' 
#' Next, a prior on current status. As in JABBA's `b.prior`, `jb3_bprior()`
#' places a lognormal prior on $B/B_{MSY}$ (or $B/K$, $F/F_{MSY}$) in a given year;
#' with `year = NULL` the prior applies to the terminal year. Here
#' $B/B_{MSY} = 0.7$ in 2017 with a CV of 0.2, the kind of prior that might come
#' from expert judgement or a length-based method:
#' 
## ----co-bprior----------------------------------------------------------------
pr.b <- jb3_priors(d0, r = pr.r, status = jb3_bprior(0.7, cv = 0.2, type = "bbmsy"))
pr.b$status
co.b <- jabba3(d0, pr.b, mod, nuts = nu, run = "Catch + B/Bmsy prior")

#' 
#' The JABBA form `status = c(0.5, 0.2, 2010, "bbmsy")` is accepted too, and
#' several priors are combined with `rbind()` (for example
#' `rbind(jb3_bprior(0.7, 0.2), jb3_bprior(1.2, 0.3, 2000, "ffmsy"))`).
#' 
#' `jbplot_priorpost_trj()` compares the two catch-only runs: trajectories, and
#' the marginal distributions of r, K and terminal $B/B_{MSY}$:
#' 
## ----co-pp-plot, fig.height = 6.5, fig.width = 9------------------------------
jbplot_priorpost_trj(list(`Prior predictive` = co.ppm, `Catch + B/Bmsy prior` = co.b))

#' 
#' The prior predictive run is the reality check: with the catch and the
#' default priors alone, bigeye would be lightly exploited ($B/B_{MSY}$ around 2
#' in 2017, $P(B < B_{MSY}) < 0.1$), with a wide interval. This is driven by the
#' default K prior (8 times the maximum catch) and initial depletion near
#' unfished: with a large K, the catch history barely dents the biomass. The
#' prior predictive run makes such implications visible before any data are
#' fitted. The terminal-year status prior pulls $B/B_{MSY}$ down to its median
#' of 0.7 and, through the catch, also the recent trajectory, K and r; but the
#' result is only as good as the prior.
#' 
#' ## Catch + auxiliary data
#' 
#' Auxiliary data enter `jb3_data()` as a data frame (one column per series, with
#' `auxiliary.type`) or as a named list of series. The observation model is
#' lognormal, with the input log-scale SE plus an additional observation SD
#' estimated per series (`jb3_model(aux.sd = NULL)`, the default) or fixed
#' (`aux.sd = 0.15`):
#' 
#' - **effort**: $\log E_t \sim N(\log F_t - \log q_E - (t - t_0)\log(1 + c),\ \sigma^2)$,
#'   where $c$ is the annual effort creep (the increase in catchability per unit
#'   effort). The creep is fixed with `jb3_model(effort.creep = 0.02)` or
#'   estimated with `effort.creep = "estimate"` and a normal prior
#'   `jb3_priors(creep = c(0.02, 0.02))`. $q_E$ is estimated.
#' - **ffmsy, bbmsy, bk**: the status quantity itself with $q = 1$ (absolute
#'   status estimates, e.g. from length-based indicators or another assessment;
#'   set `aux.q = NA` to estimate a scaling instead);
#' - **f** ($q = 1$) and **z** (estimated $q$, as JABBA).
#' 
#' The hypothetical bigeye series:
#' 
## ----co-aux-data, fig.height = 3, fig.width = 9-------------------------------
head(bet$effort_hyp, 3)
head(bet$ffmsy_hyp, 3)
ax <- rbind(data.frame(year = bet$effort_hyp$year, value = bet$effort_hyp$effort, series = "Hypothetical effort"),
            data.frame(year = bet$ffmsy_hyp$year, value = bet$ffmsy_hyp$ffmsy, series = "Hypothetical F/Fmsy"))
ggplot2::ggplot(ax, ggplot2::aes(year, value)) + ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
  ggplot2::facet_wrap(~series, scales = "free_y") + ggplot2::expand_limits(y = 0) +
  ggplot2::labs(x = "Year", y = NULL, title = "Hypothetical auxiliary data (simulated, not real data)") + theme_jabba3()

#' 
#' Catch + hypothetical effort, with the 2% creep it was simulated with, and with
#' no creep as a sensitivity:
#' 
## ----co-effort----------------------------------------------------------------
d.e <- jb3_data(bet$catch, NULL, stock = "bet", auxiliary = bet$effort_hyp, auxiliary.type = "effort")
d.e
co.e  <- jabba3(d.e, pr0, jb3_model(shape = "fox", effort.creep = 0.02), nuts = nu, run = "Catch + effort (creep 2%)")
co.e0 <- jabba3(d.e, pr0, jb3_model(shape = "fox", effort.creep = 0),    nuts = nu, run = "Catch + effort (no creep)")

#' 
#' Catch + hypothetical F/Fmsy index, 1990-2009 only, with its SE of 0.3:
#' 
## ----co-ffmsy-----------------------------------------------------------------
d.f <- jb3_data(bet$catch, NULL, stock = "bet", auxiliary = bet$ffmsy_hyp[, c("year", "ffmsy")],
                auxiliary.se = bet$ffmsy_hyp[, c("year", "se")], auxiliary.type = "ffmsy")
co.f <- jabba3(d.f, pr0, mod, nuts = nu, run = "Catch + F/Fmsy index")

#' 
## ----co-auxfits, fig.height = 3.2, fig.width = 9------------------------------
patchwork::wrap_plots(
  jbplot_auxfits(list(`creep 2%` = co.e, `no creep` = co.e0), title = "Hypothetical effort"),
  jbplot_auxfits(co.f, title = "Hypothetical F/Fmsy index"), ncol = 2, widths = c(1, 1))

#' 
#' Effort informs the *trend* in F, not its scale ($q_E$ is estimated), and the
#' creep assumption is decisive: without creep, the same nominal effort implies
#' a steadily falling F in recent years, and the stock looks lightly exploited
#' ($B/B_{MSY}$ about 2 in 2017). With the 2% creep the series was simulated
#' with, the estimate is close to the index-based assessment, which is partly by
#' construction (the effort was simulated from that assessment's F). An absolute
#' F/Fmsy index informs the *scale* of exploitation: in 1990-2009 the F/Fmsy
#' estimates of the catch + F/Fmsy run are about as precise as those of the
#' index-based run. After 2009 there is no information on F, and the terminal
#' F/Fmsy interval widens again.
#' 
#' ## Catch + index, and catch + auxiliary data + index
#' 
## ----co-index-----------------------------------------------------------------
d.i  <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet")
co.i <- jabba3(d.i, jb3_priors(d.i, r = pr.r), mod, engine = "nuts", nuts = nu, run = "Catch + index")
d.ie <- jb3_data(bet$catch, bet$cpue, bet$se, stock = "bet",
                 auxiliary = bet$effort_hyp, auxiliary.type = "effort")
co.ie <- jabba3(d.ie, jb3_priors(d.ie, r = pr.r), jb3_model(shape = "fox", effort.creep = 0.02),
                engine = "nuts", nuts = nu, run = "Catch + effort + index")

#' 
#' ## Moving up the tiers
#' 
## ----co-compare, fig.height = 7, fig.width = 10-------------------------------
co.runs <- list(`Prior predictive` = co.ppm, `Catch + B/Bmsy prior` = co.b, `Catch + effort` = co.e,
                `Catch + effort (no creep)` = co.e0, `Catch + F/Fmsy index` = co.f, `Catch + index` = co.i,
                `Catch + effort + index` = co.ie)
jabbaplot3(co.runs, panels = c("B/Bmsy", "F/Fmsy", "B/K", "kobe"), ncol = 2)

#' 
## ----co-table-----------------------------------------------------------------
pt <- jb3_partable(co.runs, format = "wide")
knitr::kable(pt[pt$Parameter %in% c("K", "r", "MSY", "qA.effort", "sdA.effort", "sdA.ffmsy", "B/Bmsy", "F/Fmsy"), ],
             row.names = FALSE, caption = "Posterior medians (95% intervals); the effort and F/Fmsy data are hypothetical")

#' 
## ----co-status, fig.height = 3.5, fig.width = 9-------------------------------
stat <- do.call(rbind, lapply(names(co.runs), function(n){ f <- co.runs[[n]]; ly <- as.character(max(f$years))
rbind(data.frame(run = n, quantity = "B/Bmsy 2017", value = f$draws$`B/Bmsy`[, ly]),
      data.frame(run = n, quantity = "F/Fmsy 2017", value = f$draws$`F/Fmsy`[, ly])) }))
stat$run <- factor(stat$run, rev(names(co.runs)))
ggplot2::ggplot(stat, ggplot2::aes(value, run)) + ggplot2::geom_vline(xintercept = 1, linetype = 2) +
  ggplot2::geom_boxplot(outlier.size = 0.3, fill = "grey85") + ggplot2::scale_x_log10() +
  ggplot2::facet_wrap(~quantity, scales = "free_x") + ggplot2::labs(x = NULL, y = NULL) + theme_jabba3()

#' 
#' The probability of being overfished, $P(B < B_{MSY})$, and of overfishing,
#' $P(F > F_{MSY})$, in 2017:
#' 
## ----co-prob------------------------------------------------------------------
knitr::kable(do.call(rbind, lapply(names(co.runs), function(n){ f <- co.runs[[n]]
data.frame(run = n, `P(B<Bmsy)` = mean(f$draws$`B/Bmsy`[, "2017"] < 1), `P(F>Fmsy)` = mean(f$draws$`F/Fmsy`[, "2017"] > 1),
           check.names = FALSE) })), digits = 2)

#' 
#' Prior against posterior for the full assessment (catch + effort + index), with
#' the prior predictive run as the reference:
#' 
## ----co-pp-full, fig.height = 6.5, fig.width = 9------------------------------
jbplot_priorpost_trj(co.ppm, co.ie, labels = c("Prior predictive", "Catch + effort + index"))

#' 
#' Moving up the tiers changes the classification of the stock. The prior
#' predictive run classifies bigeye as not overfished ($P(B < B_{MSY}) < 0.1$);
#' every run with information on the exploitation history, whether a status
#' prior, effort with creep, an F/Fmsy index or the CPUE indices, classifies it
#' as overfished with probability 0.75-0.96, with similar medians of
#' $B/B_{MSY}$ around 0.7-0.8. They differ in precision. The status prior fixes
#' the terminal estimate wherever the prior is placed, so its value lies
#' entirely in how well the prior is founded. The hypothetical F/Fmsy index
#' covers only 20 years of high exploitation, yet moves the classification from
#' "not overfished" to "overfished", as Kell et al. (2025) found for auxiliary
#' F/Fmsy indices from periods of high exploitation; without information after
#' 2009 its terminal interval stays wide. Effort carries the trend in
#' exploitation, but its result depends on the assumed creep, which should be
#' treated as a key uncertainty (fixed scenarios, or estimated with an informed
#' prior). With the CPUE indices, the effort narrows the interval of
#' $B/B_{MSY}$ somewhat, but here the agreement is partly built in, because the
#' hypothetical effort was simulated from the index-based assessment.
#' 
#' ## How much do the data add? Prior-posterior contraction
#' 
#' The prior predictive run is also the reference for a diagnostic of what each
#' tier adds. `jb3_contraction()` compares each quantity of a fit with the prior
#' predictive run (log scale):
#' 
#' - **contraction** = 1 − var(posterior)/var(prior predictive), the share of the
#'   prior uncertainty removed by the data (0: the data add nothing; 1: fully
#'   determined by the data);
#' - **shift** = (mean(posterior) − mean(prior predictive))/sd(prior predictive),
#'   how far the data move the estimate, in prior SDs.
#' 
#' Low contraction means a quantity is prior-driven; a large shift means the data
#' pull it into the tail of the prior (prior-data conflict); high contraction with
#' a moderate shift means it is informed by the data.
#' 
#' > **Source.** The diagnostic is adapted from the *posterior shrinkage* versus
#' > *posterior z-score* plot of Betancourt (2018), part of the principled Bayesian
#' > workflow of Schad et al. (2021). Betancourt computes shrinkage
#' > $= 1 - \mathrm{var(posterior)}/\mathrm{var(prior)}$ and the z-score
#' > $=(\mathrm{mean(posterior)} - \theta_{true})/\mathrm{sd(posterior)}$ over
#' > simulated data sets with a known true value $\theta_{true}$, to check a model
#' > design before it is fitted to real data. A real assessment has no true value,
#' > so JABBA3 replaces the z-score by the shift from the prior predictive mean,
#' > which measures prior-data conflict rather than calibration, and uses the prior
#' > predictive run as the reference distribution; `contraction.pdf`, against the
#' > stated prior density, is the version closest to Betancourt's shrinkage. The
#' > class names and the thresholds (`contraction.min = 0.25`, `shift.max = 2`) are
#' > JABBA3 choices, not from these papers. In simulation testing, where the truth
#' > is known, Betancourt's original form applies directly. By default the prior predictive run is made with the fit's own
#' priors, model and catch (`jb3_contraction(fit)` runs it with NUTS). Here all
#' runs are compared with the common prior predictive run `co.ppm` (catch and the
#' r, K and initial-depletion priors), so that the status prior is treated as one
#' of the tiers:
#' 
## ----co-contraction-----------------------------------------------------------
ct <- jb3_contraction(co.runs[-1], prior = rep(list(co.ppm), length(co.runs) - 1), years = c(1990, 2017))
w <- reshape(as.data.frame(ct)[c("run", "quantity", "contraction")], idvar = "quantity", timevar = "run", direction = "wide")
names(w) <- sub("^contraction\\.", "", names(w))
knitr::kable(w, digits = 2, row.names = FALSE, caption = "Contraction relative to the prior predictive run")

#' 
## ----co-contraction-plot, fig.height = 6.5, fig.width = 10--------------------
jbplot_contraction(ct)

#' 
#' The initial depletion is prior-driven in every run: no tier informs it. All
#' tiers with information on exploitation remove most of the prior uncertainty
#' in terminal status (contraction 0.8-0.9) and move it by about 1.5 prior SDs
#' towards overfished, which says that the default priors (K of 8 times the
#' maximum catch, an unfished start) are optimistic for bigeye. The tiers differ
#' in *where* they inform: the terminal status prior contracts only the terminal
#' year, whereas $B/B_{MSY}$ in 1990 is not informed at all (contraction about −0.4);
#' effort and the indices inform the whole trajectory. The hypothetical F/Fmsy
#' index informs F/Fmsy in 1990 (0.92) but not the terminal status: its posterior
#' for $B/B_{MSY}$ in 2017 is *wider* than the prior predictive (contraction
#' −0.88) and shifted by more than two prior SDs, a flag that the auxiliary data
#' move the stock into a region the rest of the model constrains poorly.
#' 
#' A catch-only fit compared with *its own* prior predictive run, which keeps the
#' status prior, has a contraction of zero by construction, and the diagnostic
#' says so:
#' 
## ----co-contraction-own-------------------------------------------------------
cb <- jb3_contraction(co.b, nuts = nu)
knitr::kable(as.data.frame(cb)[c("quantity", "prior.median", "post.median", "contraction", "contraction.pdf", "class")],
             digits = 3, row.names = FALSE)

#' 
#' `contraction.pdf` compares the posterior of r, K and initial depletion with
#' the stated prior density rather than with the prior predictive run. It shows
#' what the catch and the status prior do to the r and K priors even without any
#' data: the variance of log K is reduced by about 75% and that of log r by
#' about 35%, because only combinations under which the stock survives the catch
#' history and ends near the status prior remain (the CMSY logic). This is a check
#' of the priors (whether much of the prior is incompatible with the catch), not
#' a reason to narrow them: using the catch-conditioned r and K as priors would
#' count the catch twice.
#' 
#' ## Non-centred process deviations
#' 
#' When the data say little about the states, as in a catch-only run, the
#' process deviations and their SD form a "funnel" that NUTS samples poorly in
#' the usual (centred) form. In the non-centred form NUTS samples standard normal
#' innovations $z_t \sim N(0, 1)$ instead, and the deviations follow from
#' 
#' $$d_t = \rho\, d_{t-1} + \sigma_t \sqrt{1 - \rho^2}\, z_t,
#' \qquad \log P_t = \log \bar P_t + d_t,$$
#' 
#' with the production lag, bias correction, random walk and depletion-dependent
#' SD exactly as in the centred form (the model is the same; only the
#' coordinates of the sampler change). It is selected with
#' `jb3_model(param = "noncentred")`; the default is `"centred"`, and the
#' posterior mode always uses the centred form, where the Laplace approximation
#' needs the states themselves. `jb3_nuts_ess()` gives the bulk and tail ESS of
#' the derived states, also per 1000 draws, to compare the two:
#' 
## ----co-ess, eval = FALSE-----------------------------------------------------
## ess <- function(d, pr, p, ...) attr(jb3_nuts_ess(jabba3(d, pr, jb3_model(shape = "fox", param = p), nuts = nu, ...)), "summary")
## rbind(ess(d.i, jb3_priors(d.i, r = pr.r), "centred", engine = "nuts"), ess(d.i, jb3_priors(d.i, r = pr.r), "noncentred", engine = "nuts"),
##       ess(d0, pr.b, "centred"), ess(d0, pr.b, "noncentred"))

#' 
#' | Data | Form | min / median bulk ESS per 1000 draws, B/Bmsy | min tail ESS per 1000, B/Bmsy | Divergent transitions | Run time (s) |
#' |:--|:--|--:|--:|--:|--:|
#' | catch + index | centred | 43 / 54 | 92 | 0 | 33 |
#' | catch + index | non-centred | 113 / 327 | 42 | 1016 | 227 |
#' | catch + B/Bmsy prior | centred | 72 / 134 | 108 | 0 | 23 |
#' | catch + B/Bmsy prior | non-centred | 313 / 813 | 133 | 394 | 280 |
#' 
#' *Three chains of 2000 iterations (3000 post-warm-up draws); F/Fmsy behaves
#' the same way.*
#' 
#' The non-centred form gives a higher bulk ESS per draw, but many divergent
#' transitions, so its draws are not reliable, and it is seven to twelve times
#' slower. The reason is the production model: each state depends on *all*
#' earlier innovations through the nonlinear surplus production and the catch
#' removals, so the non-centred coordinates are strongly and nonlinearly
#' correlated, and the curvature changes sharply where the catch approaches the
#' biomass. The centred form had no divergences in either case and is the
#' default. The non-centred form remains an option for short series or very
#' weakly informed states, always to be checked with `jb3_nuts_diags()`.
#' 
## ----co-time1, include = FALSE------------------------------------------------
co.time <- round(as.numeric(difftime(Sys.time(), co.t0, units = "mins")), 1)

#' 
#' > The NUTS runs of this chapter took `r co.time` minutes in total (three
#' > chains of 4000 iterations each, run in parallel).
#' 
#' # FLR outputs
#' 
#' `jb3_flqs()` returns the time series as `FLQuants` (`stock` = B/Bmsy, `harvest`
#' = F/Fmsy, `biomass`, `F` and `catch`), with the draws as iterations when the
#' fit has draws; `jb3_flpar()` returns parameters, reference points and current
#' status as an `FLPar`.
#' 
## ----flr----------------------------------------------------------------------
fq <- jb3_flqs(fit)
names(fq)
dim(fq$stock)
jb3_flpar(fit)

#' 
#' These plot with ggplotFL and can be used in FLR workflows (e.g. as the output
#' of an estimator in `mse`):
#' 
## ----flr-plot, fig.height = 5, eval = requireNamespace("ggplotFL", quietly = TRUE)----
library(ggplotFL)
plot(fq[c("stock", "harvest", "catch")]) + ggplot2::geom_hline(yintercept = 1, linetype = 2)

#' 
#' ## FLStockR
#' 
#' `jb3_FLStockR()` converts the fit to a biomass-dynamic `FLStockR`, as
#' `FLRef::jabba2FLStockR()`: biomass in `stock` (and `ssb()`), catch, and the
#' harvest rate in `harvest`, with the draws as iterations. The `refpts` slot holds
#' $F_{MSY}$, $B_{MSY}$, MSY, $B_{thr}$ and $B_{lim}$ (by default 0.5 and 0.3
#' $B_{MSY}$, the SPiCT benchmark defaults) and $B_0$. Advice from surplus
#' production models is given in relative terms, so by default the stock is in
#' $B/B_{MSY}$ and $F/F_{MSY}$; `rel = FALSE` gives absolute biomass and harvest rate.
#' 
## ----flr-stockr, fig.height = 7, fig.width = 7, eval = requireNamespace("FLRef", quietly = TRUE)----
stk <- jb3_FLStockR(fit, ndraws = 500)        # relative (B/Bmsy, F/Fmsy) by default
stk@refpts
FLRef::plotAdvice(stk, type = "spm") + ggplot2::facet_wrap(~qname, ncol = 1, scales = "free_y")

#' 
#' Absolute biomass and harvest rate with `rel = FALSE`:
#' 
## ----flr-stockr-abs, fig.height = 7, fig.width = 7, eval = requireNamespace("FLRef", quietly = TRUE)----
FLRef::plotAdvice(jb3_FLStockR(fit, rel = FALSE, ndraws = 500), type = "spm") +
  ggplot2::facet_wrap(~qname, ncol = 1, scales = "free_y")

#' 
#' # Converting from JABBA and SPiCT
#' 
#' ## JABBA: `jbinput2jabba3()`
#' 
#' `jbinput2jabba3()` takes the object returned by `JABBA::build_jabba()` and
#' returns `list(data, priors, model)` for `jabba3()`. It carries over the priors
#' (r, K, psi to `bkinit`, shape, catchability, status priors), the production
#' lag, catch error, catchability and variance groups, `proc.dev.all`, the
#' BioProcDyn process-error options and the Selex curves, and reports each
#' difference:
#' 
#' - JABBA's inverse-gamma prior on the process variance is replaced by the
#'   FishLife prior (`keep.igamma = TRUE` keeps it);
#' - q gets a flat prior on log q (`q.uniform = TRUE` reproduces the JAGS uniform
#'   prior on q);
#' - auxiliary data are not yet supported and are dropped.
#' 
## ----conv-jabba, eval = requireNamespace("JABBA", quietly = TRUE)-------------
jb <- JABBA::build_jabba(catch = bet$catch, cpue = bet$cpue, se = bet$se, model.type = "Fox",
                         r.prior = c(0.2, 0.5), psi.prior = c(0.9, 0.25), verbose = FALSE)
options(scipen = 0)          # build_jabba() sets options(scipen = 999)
x  <- jbinput2jabba3(jb)
f2 <- jabba3(x$data, x$priors, x$model, osa = FALSE)
signif(f2$pars[c("K", "r", "sdb"), c("est", "lci", "uci")], 4)

#' 
#' > With the same process-error prior, the conversion reproduces the JABBA (RTMB engine)
#' > fit of bigeye exactly (K = 1,174,749; $B/B_{MSY}$ = 0.797 in 2017).
#' 
#' ## SPiCT: `spictinp2jabba3()`
#' 
#' `spictinp2jabba3()` takes a SPiCT input list (`inp`, as passed to
#' `spict::fit.spict()`). JABBA3 is an annual model: catches are summed by
#' calendar year and several index observations in a year are averaged
#' (geometric mean). SPiCT's observation-error form is used
#' (`obs.error = "spict"`), with `stdevfacI` as relative weights. The priors
#' `logr`, `logK`, `logbkfrac`, `logsdb`, `logn` and `logsdi` carry over with
#' their on/off switches, and a fixed `n` (phase -1) fixes the shape. The `alpha`
#' and `beta` priors (ratios of observation to process SDs) and the F process
#' have no JABBA3 equivalent and are reported as dropped.
#' 
#' A SPiCT-style input list built from the bigeye data:
#' 
## ----conv-spict---------------------------------------------------------------
ok1 <- is.finite(bet$cpue[, 2]); ok2 <- is.finite(bet$cpue[, 3])
inp <- list(obsC = bet$catch[, 2], timeC = bet$catch[, 1],
            obsI = list(bet$cpue[ok1, 2], bet$cpue[ok2, 3]),
            timeI = list(bet$cpue[ok1, 1] + 0.5, bet$cpue[ok2, 1] + 0.5),
            priors = list(logr = c(log(0.2), 0.5, 1), logn = c(log(2), 2, 0), logbkfrac = c(log(0.8), 0.5, 1)),
            ini = list(logn = log(2)), phases = list(logn = -1))
xs <- spictinp2jabba3(inp, stock = "bet")
fs <- jabba3(xs$data, xs$priors, xs$model, osa = FALSE)
signif(fs$pars[c("K", "r", "sdb"), c("est", "lci", "uci")], 4)

#' 
#' With the spict package installed, its example data convert the same way:
#' 
## ----conv-spict-pol, eval = requireNamespace("spict", quietly = TRUE)---------
## data(pol, package = "spict")
## xp <- spictinp2jabba3(pol$albacore, stock = "albacore")
## fp <- jabba3(xp$data, xp$priors, xp$model)
## fp

#' 
#' # Model equations
#' 
#' **Process equation.** Biomass relative to K, $P_t = B_t/K$, with process error
#' acting after the catch (Winker et al. 2018) and a production lag $L$:
#' 
#' $$
#' \log P_t = \log\left[P_{t-1} + \frac{r}{m-1} P_{t-1-L}\left(1 - P_{t-1-L}^{\,m-1}\right) - \frac{C_{t-1}}{K}\right] + d_t ,
#' $$
#' 
#' where $m$ is the Pella-Tomlinson shape ($m = 2$ Schaefer, $m \to 1$ Fox) and
#' $P_1$ = `bkinit`. Years before the start use the initial depletion for the
#' lagged term.
#' 
#' **Process deviations.** Stationary AR(1), with $\sigma_t$ the marginal SD:
#' 
#' $$
#' d_1 \sim N(0, \sigma_1^2), \qquad d_t \sim N\left(\rho\, d_{t-1},\ \sigma_t^2 (1-\rho^2)\right),
#' $$
#' 
#' $\rho = 0$ giving the independent process error of JABBA and SPiCT. Either
#' $\sigma_t = \sigma_{MSY}$ (constant) or $\sigma_t = \sigma_{MSY} x_t^{-k}$ with
#' $x_t = B_t/B_{MSY}$ bounded to `xlim` (dynamic).
#' 
#' **Reference points.**
#' 
#' $$
#' \frac{B_{MSY}}{K} = m^{-1/(m-1)}, \qquad F_{MSY} = \frac{r}{m}, \qquad MSY = F_{MSY} B_{MSY}.
#' $$
#' 
#' **Observations.** Index $i$ in year $t$:
#' 
#' $$
#' \log I_{i,t} = \log q_i + \log(K P_t) + \varepsilon_{i,t}, \qquad
#' \varepsilon_{i,t} \sim N\left(0,\ se_{i,t}^2 + sdi_i^2\right)
#' $$
#' 
#' (`obs.error = "additive"`), or $\varepsilon_{i,t} \sim N(0, (sdi_i\,w_{i,t})^2)$
#' (`obs.error = "spict"`). Catch: $\log C^{obs}_t \sim N(\log C_t, sdc^2)$ when
#' `catch.error = TRUE`.
#' 
#' **Estimation.** The states $\log P_t$ and catch deviations are random effects,
#' integrated by the Laplace approximation (`engine = "rtmb"`) or sampled with
#' the fixed effects (`engine = "nuts"`).
#' 
#' # Summary of checks and functions
#' 
#' | # | Check | Red flag | Action |
#' |---|---|---|---|
#' | 1 | data (`jbplot_data`) | no contrast; conflicting indices | r prior, status priors; separate `sdi`; sensitivity runs |
#' | 2 | prior and posterior (`jbplot_priorpost`) | posterior = prior; posterior at the edge | report; revisit priors |
#' | 3 | convergence (`jb3_convergence`) | gradient, Hessian, Newton decrement | restarts; priors on r or K; fewer groups |
#' | 4 | residuals (`jbplot_residuals`, `jbplot_runs`, `jb3_diags`) | trends, failed runs or Ljung-Box tests | lag + AR(1); index conflict |
#' | 5 | retrospective (`jb3_retro`) | Mohn's rho beyond -0.15 to 0.20 | recent data, q |
#' | 6 | hindcast (`jb3_hindcast`) | MASE > 1 | index without prediction skill |
#' | 7 | process deviations (`jbplot_procdev`) | long runs of one sign | `rho = 0.75`, `lag = amat` |
#' | 8 | NUTS (`jb3_nuts_diags`, `jbplot_nuts`) | $\hat{R} > 1.01$, low ESS, divergences, E-BFMI < 0.3 | longer chains; `adapt_delta`; priors |
#' 
#' | Topic | Recommendation |
#' |---|---|
#' | Process error | lag at the age at maturity + AR(1) with $\rho = 0.75$ (sensitivity 0.5, 0.9); FishLife prior on $\sigma_{MSY}$ |
#' | Shape | Pella with $B_{MSY}/K = 0.4$ unless known; convert the r prior with the shape |
#' | Engine | posterior mode for exploration, diagnostics, retro and hindcast; NUTS for final runs, and for catch-only and prior predictive runs |
#' | Data-limited | start with the prior predictive run (`prior.only = TRUE`); add status priors or auxiliary data and report how far each moves status |
#' | Simulation testing and MSE | `se = FALSE, osa = FALSE` |
#' 
#' ## Function reference
#' 
#' | Area | Functions |
#' |---|---|
#' | Data | `jb3_data()` (catch-only with `indices = NULL`; auxiliary data), `jb3_data_jabba()` |
#' | Data exploration (JARA) | `jb3_trend()` (relative joint index, census, time blocks), `jb3_jointindex()`, `jb3_catchfill()`, `jb3_fill_table()`, `jb3_qtable()` (relative catchabilities), `jbplot_trend()`, `jbplot_jointindex()`, `jbplot_catchfill()`, `jbplot_q()`, `jbplot_indexcor()` (pairwise index correlations); diagnostics with `jbplot_fits()`, `jbplot_runs()`, `jbplot_residuals()`, `jbplot_ppc()` |
#' | Priors and model | `jb3_priors()`, `jb3_bprior()` (status priors), `jb3_model()` (incl. `effort.creep`, `aux.q`, `aux.sd`, `param`) |
#' | Fitting | `jabba3()` (incl. `prior.only = TRUE`), `jb3_draws()` |
#' | Workflows | `jb3_workflow()`, `jb3_tier1()`, `jb3_tier2()` |
#' | Convergence | `jb3_convergence()` |
#' | Residual diagnostics | `jb3_osa()`, `jb3_diags()`, `jb3_runs()`, `jb3_retro()`, `jb3_hindcast()` |
#' | NUTS diagnostics | `jb3_nuts_diags()`, `jb3_nuts_ess()`, `jbplot_nuts()`, `jbplot_nutsdiags()`, `jbplot_trace()`, `jbplot_rank()`, `jbplot_mcmcdens()`, `jbplot_mcmcacf()`, `jbplot_pairs()`, `jbplot_rhat()`, `jbplot_energy()` |
#' | Projections | `jb3_project()` (F/Fmsy scenarios by default, or constant catch; optional intermediate years), `jb3_fwd()` (forecast controls, as `fw_jabba_ctrl()`) |
#' | Ensembles and HCR | `jb3_ensemble()`, `jb3_weights()`, `jb3_hcr()`, `jbplot_hcr()`, `jb3_fwdtable()` (GFCM fishing opportunities) |
#' | FLR | `jb3_flqs()`, `jb3_flpar()`, `jb3_FLStockR()` |
#' | Life history | `jb3_leslie()`, `jb3_fmsy2r()`, `jb3_r2fmsy()`, `jb3_rho_gt()`, `jb3_rho_table()`, `jb3_selex()`, `jb3_sel5()`, `jbplot_selex()` |
#' | Tables | `jb3_partable()`, `jb3_tstable()`, `jb3_projtable()`, `jb3_kobematrix()`, `jb3_nuts_diags()` |
#' | FLR outputs | `jb3_flqs()`, `jb3_flpar()` |
#' | Converters | `jbdat2flr()`, `jbcatch2flq()`, `jbidx2FLIndices()`, `flr2jbdat()`, `jbinput2jabba3()`, `spictinp2jabba3()` |
#' | Graphics: data | `jbplot_data()`, `jbplot_indices()` |
#' | Graphics: estimates | `jbplot_advice()`, `jbplot_trj()`, `jbplot_kobe()`, `jbplot_biplot()`, `jbplot_procdev()`, `jbplot_catch()`, `jbplot_catcherror()`, `jbplot_spphase()`, `jbplot_spdyn()`, `jbplot_production()`, `jbplot_ebsb_status()` |
#' | Graphics: priors | `jbplot_priorpost()`, `jbplot_priorpost_trj()` (prior predictive vs posterior), `jbplot_bprior()` |
#' | Graphics: diagnostics | `jbplot_fits()`, `jbplot_auxfits()`, `jbplot_residuals()`, `jbplot_stdresiduals()`, `jbplot_runs()`, `jbplot_spdiags()`, `jbplot_ppc()`, `jbplot_retro()`, `jbplot_hcxval()` |
#' | Graphics: runs and projections | `jbplot_steps()`, `jbplot_ensemble()`, `jbplot_summary()`, `jbplot_projection()`, `jbplot_crisk()`, `jbplot_hcr()` |
#' | Graphics: overview | `plot()` / `jabbaplot3()` (3 x 3 panels, single fit or runs) |
#' | Graphics: all | `jbplot_all()` |
#' | Graphics helpers | `theme_jabba3()`, `scale_x_year()`, `jb3cols()` |
#' | Information | `jabba3_info()` |
#' 
#' *Planned: multi-fleet Selex, auxiliary-data lags (JABBA `auxiliary.lag`), a worked
#' simulation-testing example.*
#' 
#' # References
#' 
#' Betancourt, M. (2018). Calibrating model-based inferences and decisions. arXiv:1803.08393. <https://arxiv.org/abs/1803.08393>
#' 
#' Carvalho, F., Winker, H., Courtney, D., Kapur, M., Kell, L., Cardinale, M.,
#' Schirripa, M., Kitakado, T., Yemane, D., Piner, K.R., Maunder, M.N., Taylor, I.,
#' Wetzel, C.R., Doering, K., Johnson, K.F. and Methot, R.D. (2021). A cookbook for
#' using model diagnostics in integrated stock assessments. *Fisheries Research*
#' 240: 105959.
#' 
#' Conn, P.B. (2010). Hierarchical analysis of multiple noisy abundance indices. *Canadian Journal of Fisheries and Aquatic Sciences*, 67, 108-120. doi:10.1139/F09-175
#' 
#' Francis, R.I.C.C. (2011). Data weighting in statistical fisheries stock assessment models. *Canadian Journal of Fisheries and Aquatic Sciences*, 68, 1124-1138.
#' 
#' Kell, L.T., Cardinale, M., Mosqueira, I., Winker, H. and Sharma, R. (2025).
#' Transforming data and capacity-limited stock assessment: achieving risk
#' equivalence with hierarchical assessment frameworks and auxiliary data.
#' *Sustainability* 17: 9383. doi:10.3390/su17219383
#' 
#' Kell, L.T., Mosqueira, I., Grosjean, P., Fromentin, J.-M., Garcia, D., Hillary,
#' R., Jardim, E., Mardle, S., Pastoors, M.A., Poos, J.J., Scott, F. and Scott,
#' R.D. (2007). FLR: an open-source framework for the evaluation and development of
#' management strategies. *ICES Journal of Marine Science* 64(4): 640-646.
#' 
#' Kell, L.T., Sharma, R., Kitakado, T., Winker, H., Mosqueira, I., Cardinale, M.
#' and Fu, D. (2021). Validation of stock assessment methods: is it me or my model
#' talking? *ICES Journal of Marine Science* 78(6): 2244-2255.
#' 
#' Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H. and Bell, B.M. (2016). TMB:
#' automatic differentiation and Laplace approximation. *Journal of Statistical
#' Software* 70(5): 1-21.
#' 
#' Monnahan, C.C. and Kristensen, K. (2018). No-U-turn sampling for fast Bayesian
#' inference in ADMB and TMB: introducing the adnuts and tmbstan R packages. *PLoS
#' ONE* 13(5): e0197954.
#' 
#' Pedersen, M.W. and Berg, C.W. (2017). A stochastic surplus production model in
#' continuous time. *Fish and Fisheries* 18(2): 226-243.
#' 
#' Schad, D.J., Betancourt, M. and Vasishth, S. (2021). Toward a principled Bayesian workflow in cognitive science. *Psychological Methods*, 26, 103-126. Preprint: <https://arxiv.org/abs/1904.12765>
#' 
#' Sherley, R.B. et al. (2020). *Ecology and Evolution* (JARA census model applied to absolute abundance indices). doi:10.1002/ece3.6554
#' 
#' Thorson, J.T. (2020). Predicting recruitment density dependence and intrinsic
#' growth rate for all fishes worldwide using a data-integrated life-history model.
#' *Fish and Fisheries* 21(2): 237-251.
#' 
#' Thorson, J.T., Shelton, A.O., Ward, E.J. and Skaug, H.J. (2015). Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. *ICES Journal of Marine Science*, 72, 1297-1310.
#' 
#' Thygesen, U.H., Albertsen, C.M., Berg, C.W., Kristensen, K. and Nielsen, A.
#' (2017). Validation of ecological state space models using the Laplace
#' approximation. *Environmental and Ecological Statistics* 24(2): 317-339.
#' 
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: an improved $\hat{R}$ for
#' assessing convergence of MCMC. *Bayesian Analysis* 16(2): 667-718.
#' 
#' Winker, H. et al. (2026). Process error in surplus production models:
#' scale, persistence and inflation (BioProcDyn). In preparation.
#' 
#' Winker, H., Carvalho, F. and Kapur, M. (2018). JABBA: Just Another Bayesian
#' Biomass Assessment. *Fisheries Research* 204: 275-288.
#' 
#' Winker, H., Thorson, J.T., Fairweather, T., Leslie, R. and Durholtz, D. (2017). Towards improving precision in South African demersal trawl survey indices using geostatistical GLMMs. MARAM/IWS/2017/Hake/BG8. University of Cape Town.
#'