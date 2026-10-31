# Along-Track Speed Budget — reproduction and verification code

MATLAB code accompanying the manuscript:

> **Along-Track Speed Budget Governs Formation Error in Saturated UAV Swarms:
> Closed-Form Corner and Curvature Criteria**
> Guibang Ma, Danni Zhang
> School of Information Engineering, Shenyang University of Chemical Technology,
> Shenyang 110142, China

---

## What this code does

1. **Reproduces** the ARD-PF (attractive–repulsive distance-based potential field) swarm
   simulator and the delay–stability experiments of the source study.
2. **Verifies** every numerical claim in the paper through a single self-checking acceptance
   suite. The suite recomputes each comparison quantity from scratch rather than quoting
   stored values, so no reported number depends on a cached figure.
3. **Generates** all figures of the manuscript and the graphical abstract.

The three parts correspond to the three entry points: `src/acceptance_test.m`,
`figures/make_figures.m`, `figures/make_fig7.m`, `figures/make_graphical_abstract.m`.

---

## Requirements

| | |
|---|---|
| MATLAB | R2019b or later (developed and tested on R2025b) |
| Statistics and Machine Learning Toolbox | required — used only for `corr` in the turn-reversal test |

No other toolboxes are needed.

---

## How to run

```matlab
cd <this folder>
run_all
```

`run_all` executes the whole chain — acceptance suite → manuscript figures → graphical
abstract — and takes roughly **15–25 minutes** on a desktop machine, most of it in the
trajectory simulations.

To run only the acceptance suite:

```matlab
addpath(genpath(pwd));
acceptance_test
```

### Expected result

The suite prints a PASS/FAIL table and writes `acceptance_result.txt`. **All 13 checks
pass** (three reproduction checks P1–P3 plus T1–T10). A reference output from MATLAB R2025b is shipped in `results/`.

---

## Files

### `src/` — model and verification

| File | Purpose |
|---|---|
| `ardpf_run.m` | ARD-PF leader–follower simulator under uniform communication delay. Reproduces the source study's implementation, including its one-step leader-position lag in the attraction term and its repulsion direction convention. Accepts optional `Katt` (attraction gain) and `Kff` (velocity-feedforward gain) arguments for the controller-variant sweep. |
| `controller_variants.m` | re-measures the zero-delay gate under the controller variants of Table 15 (repulsion form, attraction gain, velocity feedforward) |
| `dynamics_sweep.m` | re-measures the zero-delay gate under a first-order speed-loop (rotorcraft) model with time constant `tau_v` |
| `delay_sweep.m` | re-measures the geometry gate across the communication-delay range $\tau\in\{0,\dots,15\}$ s, showing it is delay-invariant for $\tau\le5$ s |
| `circle_headroom.m` | shows that on a circular arc the curvature error vanishes once $v_{\max}\ge\max_i S_i$ (Eq. 13), reducing the error from 11.2 m to the baseline $b_0$ |
| `gate_vs_angle.m` | evaluates the gate for several corner angles and verifies the $M{=}6$ (60°) prediction $\overline{\lvert d\rvert}^\star=17.1$ m against the measured crossing |
| `mixpoly_check.m` | verifies the mixed corner–arc trajectory: the unified map reproduces the full simulation to 0.6%, while the mixed error exceeds either endpoint (corner/rotation coupling) |
| `mixed_check.m` | the acceptance primitive for the mixed trajectory (two $90°$ corners + two $90°$ arcs, $\rho=12$ m): returns full-simulation, unified discrete-map and residual, wired into the acceptance suite as T10 |
| `stability_proof_check.m` | numerical verification of Appendix C: re-derives the per-follower peaks and $A$, and checks each lemma and contraction bound of Theorem C.1 / Proposition C.8 |
| `mixpoly_table.m` | arc-length parameterisation of the mixed corner/arc trajectory |
| `sim_instr.m` | instrumented simulator returning per-follower error histories; pointwise-identical to `ardpf_run` (max deviation ≈ 3 × 10⁻¹⁶) |
| `traj_kappa.m` | single shared implementation of trajectory geometry and curvature |
| `reduced_ode.m` | single shared integrator for the reduced-order single-follower ODE; returns **both** full-run and last-lap statistics so that a model value and its reference can never be compared across windows by accident |
| `A_closedform.m` | corner accumulation factor *A*, from relaxation along a conserved parabola |
| `A_headroom_closedform.m` | headroom generalisation of *A*: the closed form $A(\delta)$ for $v_{\max}>v$ (Appendix D), from relaxation along the conserved curve $q/P^{1+2/\delta}$ |
| `generalized_A.m` | composite-map fixed point for polygons with unequal sides and turns |
| `unified_discrete.m` | unified discrete map covering all eight trajectory families |
| `decompose_independent.m` | three-term decomposition in which every term comes from an independent model |
| `acceptance_test.m` | **verification entry point** — runs all checks |

### `figures/`

| File | Purpose |
|---|---|
| `make_figures.m` | Figures 1–6, 8 and 9 of the manuscript (400 dpi PNG + vector EPS) |
| `make_fig7.m` | Figure 7 (speed-headroom U-curve, 400 dpi PNG + vector EPS) |
| `make_graphical_abstract.m` | graphical abstract (2200 × 1120 px and 4400 × 2240 px) |

### `results/`

Reference outputs from a run on MATLAB R2025b, provided so that readers can compare
without executing the full chain.

---

## Notes on reproducibility

- **Determinism.** At zero delay the delay draw collapses to a constant, so the simulator is
  deterministic and independent of the random seed; results are bit-reproducible.
- **Caching.** `make_figures` caches simulation data in `fig_data.mat` on first run. Delete
  that file to force recomputation.
- **Two comparison windows.** The paper reports both full-run `max_t e_form` and last-lap
  steady state. `reduced_ode` returns both explicitly; mixed use of the two is a documented
  source of apparent disagreement.

---

## License

Released under the MIT License — see `LICENSE`.

## Citation

If you use this code, please cite the manuscript. Machine-readable metadata is in
`CITATION.cff`.
