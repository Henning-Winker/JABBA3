# JABBA3 design notes (draft, 29 Oct 2026)

Working notes for the structure of JABBA3. Items marked **[open]** need a decision before the code that depends on them
is written.

## Goals
1. Same model family as JABBA 2.x (Pella-Tomlinson state-space surplus production with priors), estimated in RTMB, with
   NUTS as the fully Bayesian option. No JAGS.
2. One harmonised set of parameter, argument and output names, aligned with FLicc (and FLicc-SPM).
3. FLR-native inputs and outputs; ggplot2 graphics from JAGG; classic JABBA inputs via a helper.
4. Built-in validation against JABBA 2.x so results are traceable for existing users.

## Decisions so far (29 Oct 2026)
- Follow the FLicc pattern: FLR S4 classes for inputs and outputs (FLQuant / FLQuants / FLIndices in, FLQuants / FLPar out),
  the fit itself a plain list with an S3 tag (`jabba3`) holding engine internals (RTMB objective, optimiser result,
  stanfit, settings).
- Inputs: `catch` and `index` (FLIndexBiomass in FLIndices).
- Error terms named as in SPiCT: `sdb`, `sdi`, `sdc` (and `sdf` if an F process is ever added).
- Priors as one named list; every prior has an on/off switch. Default: run with the r prior.
- No dplyr. Base R in the core; data.table only where it pays (long tables of draws, covariance or MCMC summaries).

## Architecture

### 1. Data
- Native: catch as `FLQuant` (or `FLQuants` by fleet, needed for Selex), indices as `FLIndices` (index, `index.var` for
  SEs, selectivity/maturity references optional).
- Legacy helper: `jb3_data_jabba(catch, cpue, se, ...)` converts the JABBA 2 data frames to the same internal data object.
- Internal: one plain list consumed by the model (no FLR objects inside the objective), so FLR stays optional at fit time.

### 1b. What FLIndexBiomass can store (checked against FLCore 2.6.33)
- `index`: the index (quant "all"). `index.var`: one value per year; use it for the fixed log-scale SD of the index
  (documented convention) with `distribution = "lognormal"`. `index.q`: starting or fixed catchability.
- `sel.pattern` exists but is fixed to quant "all" in FLIndexBiomass (1 x years), so it cannot hold a selectivity curve.
  There is no `catch.sel` slot. Selex EB/SB curves (or selectivity-at-age) therefore need either an age-based `FLIndex`
  (sel.pattern by age) or a separate argument, as in JABBA 2 (`index.ebsb`). **[open]**
- Catch by fleet as `FLQuants`; fleet selectivity likewise from a separate argument (or an FLBRP) for Selex.

### 1c. Error terms (SPiCT naming)
| JABBA3 | SPiCT | JABBA 2 | Meaning |
|---|---|---|---|
| `sdb` | `sdb` | `sigma.proc` | process error SD of log biomass (constant, or at Bmsy when depletion-dependent) |
| `sdi` | `sdi` | `tau` (sigma.est) | estimated additive index observation SD, one per index (or per variance group); total = sqrt(index.var^2 + sdi^2) |
| `sdc` | `sdc` | `catch.cv` | catch observation SD (off = catch known, as default in JABBA 2 sims; on = catch with error) |
| `sdf` | `sdf` | (inactive "fdev") | SD of an F random walk; only if a catch-conditioned-on-F formulation is added later |

### 1d. Index observation error: two forms (`obs.error`)
- `"additive"` (**default, as JABBA 2; SS3-like extra SD**): total SD_t = sqrt(se_t^2 + sdi^2), se_t from `index.var`
  (or the JABBA fixed minimum, 0.1 without SEs, 0.001 added with SEs); `sdi` estimated per index or variance group,
  no prior (JABBA's diffuse gamma on 1/sdi^2 kept for continuity in the JAGS-like setting).
- `"spict"` (option): total SD_t = sdi x w_t, `sdi` estimated per index with no floor; `index.var` enters only as relative
  year weights w_t (normalised to mean 1), as SPiCT's `stdevfacI`; optional lognormal prior on `sdi` via `priors$sdi`
  (off by default, as in the SPiCT benchmarks).

### 2. Specification
- `jb3_priors()`: one named list, each entry `c(median, log.sd, on)` (lognormal; beta allowed for bkinit), every prior with an
  on/off switch. **JABBA stays JABBA:** defaults follow JABBA 2, with one change, the process error prior.

| Prior | Meaning | Default (as JABBA 2 unless noted) |
|---|---|---|
| `r` | intrinsic growth rate | **on** (user value, FishLife or from Fmsy) |
| `K` | carrying capacity | on, weak (JABBA: 8 x max catch, CV 1) |
| `bkinit` | initial depletion B1/K (JABBA 2: `psi`, accepted as alias by the legacy helper; SPiCT: `bkfrac`) | on (JABBA: 0.9, CV 0.25) |
| `m` / `BmsyK` | shape | fixed from BmsyK (0.4 default); on to estimate m |
| `q` | catchability | off (flat on log q) |
| `sdb` | process error SD | **on: LN(0.80 x sigR x r^0.45, 0.3)**, r = median of the r prior, sigR = 0.5 unless given (FishLife relation for sd at Bmsy). Replaces igamma(0.001, 0.001), whose right-skewed posterior collapses the mode towards zero. |
| `sdi` | index SD | estimated, no prior; additive to index.var by default, or SPiCT form (`obs.error = "spict"`) with optional LN prior |
| `sdc` | catch SD | JABBA catch CV when catch error is on |
| `kdyn` | depletion dependence of sdb (Tier 2) | on when the dynamic form is used (life-history centre, sd 0.15) |

  As in the SPiCT benchmarks (r prior on; priors on sdi, sdb, alpha and beta switched off), no ratio priors (alpha = sdi/sdb,
  beta = sdc/sdf) are used.
  To check before release: the FishLife centre is 1.2-3.8x the diffuse-prior JAGS medians on the BioProcDyn case-study OMs
  (e.g. cod 0.24 vs 0.064, sardine 0.44 vs 0.18); the JABBA 2 soft upper bound sigmaproc_bound = 0.2 conflicts with
  centres above 0.2 and should go (or be raised).
- `jb3_model()`: structure. Shape (Schaefer / Fox / Pella / estimated m), production lag, process variance form
  (constant / floor with k / dynamic), Selex EB/SB curves by fleet, catch error, soft bounds.
- Settings separate from priors, mirroring FLicc's `settings` list.

### 3. Engines
- `engine = "rtmb"`: posterior mode of the fixed effects, states (and catch deviations) Laplace-integrated; uncertainty
  from the joint precision (Wald on the log scale) or pseudo-posterior draws.
- `engine = "nuts"`: tmbstan on the same objective, all states sampled.
- Step 1 of the Tier 2 two-step uses an informative lognormal prior on sigma.proc (centre 2 x life-history floor, log sd
  0.3): with it the RTMB mode reproduces the JAGS posterior median (BioProcDyn `output_rtmb/RESULTS_RTMB.md`, section 6).

### 3b. Fitting workflows: Tier 1, Tier 2, Selex (automatic or manual)
- Manual: everything is set in `jb3_model()` (lag, process-variance form and level, k prior, Selex curves) and fitted once.
- Automatic: `workflow =` in the fit call runs the steps internally and returns all fits with the final one on top:
  - `"tier1"`: lag at age at maturity.
  - `"tier2"`: step 1 = constant sdb (FishLife sdb prior), step 2 = depletion-dependent sdb with the level tuned on step 1
    (sd_msy = sdb1 / sqrt(mean(x^-2k)) over the step-1 B/Bmsy) and k estimated, plus the lag.
  - `"selex"` / `"selex_tier2"`: the same on the SSB state with EB/SB-mapped index and removals.
  - `"protocol"`: chooses by selectivity relative to maturity (paper protocol: Lc < Lm50 -> Tier 2 + post-hoc EB-to-SSB
    conversion; at/after maturity or dome -> Selex + Tier 2).
- Each step is also available on its own (e.g. `jb3_tier2(fit1)` takes any step-1 fit), so the steps can be run, inspected
  and adjusted one at a time. With RTMB a full Tier 2 workflow is ~0.5 s, so the automatic route is cheap enough for sims.

### 3c. Optimisation and convergence (as in FLicc)
- nlminb, then Newton polishing (up to 3 steps; Hessian by finite differences of the gradient with random effects; steps
  that do not lower the objective are rejected), then gradient-guarded restarts (accepted only if the gradient falls and
  the objective does not rise), and a final polish.
- `jb3_convergence(fit)`: table of checks with plain-language interpretation (optimiser code, max gradient, Hessian
  positive definite, SEs available, eigenvalues) and an overall tier from the Newton decrement 0.5 g'H^-1 g:
  ok < 0.001, warning <= 0.1, fail above 0.1 or with a non-positive-definite Hessian.
- Simulation outputs carry the tier, so failed fits (e.g. near-unidentified state paths with short-lived stocks and a lag)
  are flagged rather than silently averaged.

### 4. Outputs
- `FLQuants`: stock (B/Bmsy), harvest (F/Fmsy), biomass, catch (observed / estimated), EB for Selex; iterations = draws
  (1 for mode-only fits, with log-scale SEs attached).
- `FLPar` reference points: K, r, m, Bmsy, Fmsy, MSY, B/K and status in the terminal year.
- Diagnostics: index fits, ordinary and OSA residuals (particle filter at the mode), runs tests, retrospective and
  hindcast MASE.
- A light S4 or S3 fit class wrapping these, with `summary()`, `plot()` and `as.FLQuants()` methods. **[open]**

### 4b. Diagnostics routine: `jb3_diags(fit)` = JABBA set + SPiCT set
One routine, run on every fit (fast parts) with the refit-based parts on request.
- **OSA residuals by default** (particle filter at the mode, 0.04-0.09 s per fit; calibrated: mean 0, sd 1, ACF ~0,
  runs-test rejection 2% on simulated data). Ordinary residuals from smoothed states are shrunk and autocorrelated by
  construction, so tests on them are biased; OSA residuals are the primary basis for residual tests.
- **JABBA set** (kept for continuity and ss3diags compatibility): ordinary index residuals, RMSE, SDNR, runs test with
  three-sigma limits (run on OSA residuals by default, on ordinary residuals as "classic"), retrospective with Mohn's rho
  (`retro = n`), hindcast cross-validation with MASE (`hindcast = n`), process-error deviations.
- **SPiCT set** (on OSA residuals, as `plotspict.diagnostic`): per index, bias t-test, Ljung-Box (lag 4), Shapiro-Wilk;
  panels residuals over time / ACF / QQ with the p-values in the titles (flagged when < 0.05); catch too when catch error
  is on.
- **Posterior predictive check only from full posteriors** (NUTS or JAGS). With the RTMB mode the MVN draws would test the
  Laplace approximation rather than the model, so PPC is skipped (an approximate version only on explicit request).
- Convergence verdict (`jb3_convergence`) included, and a one-table summary of all tests with pass/flag.
- Power is limited: on simulated data the OSA residuals did not detect an omitted 3-year production lag that biased
  terminal status by +0.14. Residual tests check the error assumptions; they are not structural diagnostics.

### 5. Graphics
- JAGG functions moved in (advice plot after `FLRef::plotAdvice()`, fits, residuals, runs test, Kobe, retro, hindcast),
  reading from the FLR outputs through one tidy layer.

## Open decisions
- **[open] RTMB or compiled TMB.** RTMB: no compiler, model readable in R, fast enough for surplus production models.
  Compiled TMB: matches FLicc's build setup. Leaning RTMB.
- **[open] Relation to FLicc-SPM.** JABBA3 as the surplus production side of the FLicc family, or a shared production /
  selectivity core used by both.
- Fit class: decided, FLicc pattern (see above). For reference, Laurie Kell's FLABBA prototype (2020) defines an S4 class
  `biodyn` extending FLComp with slots catch (FLQuant), indices and indices.se (FLQuants), stock, diags, priors (array); the
  rest is a copy of JABBA 1.x with base graphics and no dplyr.
- **[open] Version 1 scope.** Proposed: single stock; single- and multi-fleet Selex; lag, dynamic process variance (Tier 2)
  and sigma.proc prior; projections; FLR outputs; graphics; validation. Later: JARA, auxiliary data, catch-only.
- **[open] Authors and licence.** EUPL (>= 1.2) as JABBA/JAGG; author list to complete.

## Validation plan
- Reference: JABBA 2.3.7 (JAGS) with BioProcDyn patches 0001-0015.
- ICCAT example stocks (Schaefer, Fox, Pella, lag): parameters and status within MCMC noise of JABBA 2.
- BioProcDyn case-study OMs (9 stocks x 20 iterations) and the 81-OM grid: bias and coverage tables for current practice,
  Tier 1, Tier 2, Selex and Selex + Tier 2, RTMB and NUTS against JAGS.
- OSA residuals calibrated on simulated surplus production data (N(0,1), white).
- These become the `tests/` suite: a structural change is accepted only if it reproduces them.

## Roadmap
1. Core: port the RTMB objective (BioProcDyn `jabba_patch/0013`), legacy data helper, mode + draws + NUTS. Test against
   JABBA 2 on the ICCAT examples.
2. FLR interface: FLQuant/FLIndices inputs, FLQuants/FLPar outputs, `mse` estimator wrapper.
3. Graphics and diagnostics from JAGG; OSA residuals; hindcast/retro.
4. Validation suite (ICCAT, case-study OMs); vignettes; pkgdown site.
