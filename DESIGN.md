# JABBA3 design notes (draft, 29 Oct 2026)

Working notes for the structure of JABBA3. Items marked **[open]** need a decision before the code that depends on them
is written.

## Goals
1. Same model family as JABBA 2.x (Pella-Tomlinson state-space surplus production with priors), estimated in RTMB, with
   NUTS as the fully Bayesian option. No JAGS.
2. One harmonised set of parameter, argument and output names, aligned with FLicc (and FLicc-SPM).
3. FLR-native inputs and outputs; ggplot2 graphics from JAGG; classic JABBA inputs via a helper.
4. Built-in validation against JABBA 2.x so results are traceable for existing users.

## Architecture

### 1. Data
- Native: catch as `FLQuant` (or `FLQuants` by fleet, needed for Selex), indices as `FLIndices` (index, `index.var` for
  SEs, selectivity/maturity references optional).
- Legacy helper: `jb3_data_jabba(catch, cpue, se, ...)` converts the JABBA 2 data frames to the same internal data object.
- Internal: one plain list consumed by the model (no FLR objects inside the objective), so FLR stays optional at fit time.

### 2. Specification
- `jb3_priors()`: one named list. K, r (or derived from Fmsy), psi, shape (Bmsy/K), q, sigma.proc (lognormal by default,
  life-history helper), sigma.obs (fixed + estimated), k (depletion dependence), catch CV.
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
- **[open] Fit class.** S4 (FLR style, e.g. building on Laurie Kell's FLABBA prototype `jabba-class.R`) or S3 list (FLicc
  style).
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
