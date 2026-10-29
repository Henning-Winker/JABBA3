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

### 2. Specification
- `jb3_priors()`: one named list, each entry `c(median, log.sd, on)` on the natural scale with a lognormal prior (as in
  SPiCT, where each prior is c(mean on log scale, sd, use flag)). Proposed entries and defaults:

| Prior | Meaning | Default | SPiCT name (default) |
|---|---|---|---|
| `r` | intrinsic growth rate | **on** (user, FishLife or from Fmsy) | `logr` (off) |
| `K` | carrying capacity | off (soft bounds only) **[open]**; JABBA 2 always had a weak K prior (CV 1) | `logK` (off; start 4 x max catch) |
| `psi` | initial depletion B1/K | on, JABBA 2 default (0.9, CV 0.25) **[open]** | `logbkfrac` (off) |
| `m` / `BmsyK` | shape | off: shape fixed from BmsyK = 0.4; on to estimate | `logn` (on: log 2, sd 2) |
| `q` | catchability (per index) | off | `logq` (off) |
| `sdb` | process error SD | on: LN(2 x life-history floor, 0.3) | `logsdb` (off) |
| `sdi` | additive index SD | off | `logsdi` (off) |
| `sdc` | catch SD | fixed when catch error is on | `logsdc` (off) |
| `sdi/sdb` | ratio (stabiliser) | off; optional | `logalpha` (on: log 1, sd 2) |
| `kdyn` | depletion dependence of sdb (Tier 2) | on when the dynamic form is used (life-history centre, sd 0.15) | (none) |

  SPiCT estimates K freely by default: only `logn`, `logalpha` and `logbeta` (= sdc/sdf) carry default priors, all wide
  (sd 2); everything else is available but switched off. The JABBA3 defaults stay closer to JABBA (r and psi on, shape
  fixed) and add the informative `sdb` prior that makes the posterior mode reliable.
  Name for the depletion-dependence exponent: `kdyn` (not `k`, to avoid confusion with K). **[open]**
- `jb3_model()`: structure. Shape (Schaefer / Fox / Pella / estimated m), production lag, process variance form
  (constant / floor with k / dynamic), Selex EB/SB curves by fleet, catch error, soft bounds.
- Settings separate from priors, mirroring FLicc's `settings` list.

### 3. Engines
- `engine = "rtmb"`: posterior mode of the fixed effects, states (and catch deviations) Laplace-integrated; uncertainty
  from the joint precision (Wald on the log scale) or pseudo-posterior draws.
- `engine = "nuts"`: tmbstan on the same objective, all states sampled.
- Step 1 of the Tier 2 two-step uses an informative lognormal prior on sigma.proc (centre 2 x life-history floor, log sd
  0.3): with it the RTMB mode reproduces the JAGS posterior median (BioProcDyn `output_rtmb/RESULTS_RTMB.md`, section 6).

### 4. Outputs
- `FLQuants`: stock (B/Bmsy), harvest (F/Fmsy), biomass, catch (observed / estimated), EB for Selex; iterations = draws
  (1 for mode-only fits, with log-scale SEs attached).
- `FLPar` reference points: K, r, m, Bmsy, Fmsy, MSY, B/K and status in the terminal year.
- Diagnostics: index fits, ordinary and OSA residuals (particle filter at the mode), runs tests, retrospective and
  hindcast MASE.
- A light S4 or S3 fit class wrapping these, with `summary()`, `plot()` and `as.FLQuants()` methods. **[open]**

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
