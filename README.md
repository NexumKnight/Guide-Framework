# GUIDE

**Generalized Uncertainty Inference for Drag Exit-Times**

Closed-form first-exit-time CDFs for LEO station-keeping decisions under thermospheric density uncertainty.

GUIDE replaces Monte Carlo screening for maneuver timing with an analytical cumulative distribution function. It models the fractional thermospheric density perturbation as an Ornstein-Uhlenbeck (OU) process; the resulting semi-major axis (SMA) error is a scalar Gaussian process whose first-exit time through a station-keeping control band has a closed form. That CDF drives maneuver timing and delta-v budgeting at a fraction of the cost of a Monte Carlo screen.

## The pipeline

GUIDE is the decision layer of a three-stage pipeline:

```
   ROPE                    DOPE                     GUIDE
   ----                    ----                     -----
   probabilistic    -->    augmented-state   -->    closed-form first-exit
   density emulator        filter                   decision layer
   (uncertainty model)     (density estimate)       (this package)
```

ROPE supplies the density uncertainty model. DOPE supplies the current density
estimate from tracking data. GUIDE turns those into an operational
station-keeping decision: when to maneuver, how much warning is available, and
what the recentre cost is.

## What GUIDE computes

Given a spacecraft, a density driver, and a control band, GUIDE returns:

- the **first-exit-time CDF** `F(t) = P(band exit by time t)`, in closed form,
- a **maneuver trigger time** for a chosen risk tolerance (for example, schedule
  the burn when the exit probability first reaches 5 percent),
- the **median exit time** and the **actionable lead time** between trigger and
  median,
- the **tangential delta-v** to recentre the SMA error.

The object GUIDE tracks is the *uncertainty* in the SMA error driven by
unmodeled density fluctuations, not the deterministic decay. The deterministic
decay is known and can be planned around; what forces an unscheduled maneuver is
the accumulated density-driven error breaching the band.

## Installation

```bash
git clone https://github.com/NexumKnight/guide-framework.git
cd guide-framework
pip install -e .
```

Runtime dependencies are NumPy and SciPy. To run the test suite, install the
test extra:

```bash
pip install -e ".[test]"
pytest
```

## Quickstart

```python
import numpy as np
from guide import Spacecraft, DensityDriver, ControlBand, GuideProblem, decide

problem = GuideProblem(
    spacecraft=Spacecraft.from_physical(
        cd=2.2, area_m2=5.0, mass_kg=473.0, sma_m=6_871_000.0
    ),
    driver=DensityDriver(rho_nom=5.0e-12, gamma=0.45, tau=3600.0),
    band=ControlBand(half_width=120.0),
)

decision = decide(problem, risk_tolerance=0.05)
print(f"maneuver trigger : {decision.trigger_time_s / 3600:.2f} h")
print(f"median band exit : {decision.median_exit_s / 3600:.2f} h")
print(f"lead time        : {decision.lead_time_s / 3600:.2f} h")
print(f"recentre delta-v : {decision.delta_v_mps * 1000:.2f} mm/s")
```

Two runnable scripts are in `examples/`:

- `examples/quickstart.py`: the high-level `decide()` path on a moderately
  active day.
- `examples/gannon_kanopus.py`: a Gannon-storm screening for a KANOPUS-V class
  satellite, with a Monte Carlo cross-check that confirms the closed form is a
  conservative bound on exit probability.
- `examples/rope_to_guide.py`: calibrating the OU driver from a ROPE forecast
  (runs on a synthetic stand-in; the real ROPE wiring is documented inline).

## ROPE integration

GUIDE's density inputs (`rho_nom`, `gamma`, `tau`) can be calibrated from ROPE
(https://github.com/anirudh53/ROPE_META), the lab's reduced-order probabilistic
density emulator. ROPE returns a mean density field and an Unscented-Transform
uncertainty, both `(H, 72, 36, 45)` over local solar time, latitude, and
altitude, in kg/m^3, and ships a `DensityInterpolator` that samples them along a
satellite track as `{density, sigma}`.

`guide.drivers.rope_adapter` turns that along-track series into a `DensityDriver`:

- `gamma` is the stationary fractional 1-sigma, from `sigma_rho / rho_nom`.
- `rho_nom` is a representative nominal density collapsed from the track.
- `tau` is the correlation time. ROPE's `density_std` is an instantaneous
  spread and does not by itself determine `tau`; the adapter estimates `tau`
  from a realized error path (the residual against accelerometer truth, or the
  ROPE ensemble-member deviations), and falls back to the nominal 3600 s with a
  warning if none is supplied. Calibrating `tau`, and promoting it to a
  time-varying `tau(t)` for storms, is the time-inhomogeneous OU extension.

GUIDE does not import ROPE, TensorFlow, or PyTorch: the adapter operates on the
plain NumPy series you extract from a ROPE run, so GUIDE stays installable
without ROPE's stack. Note that running ROPE itself requires its downloaded
model weights (hosted on Google Drive, placed under `ROPE_META/Models/`); the
drivers go back to 1957, so Halloween 2003, St Patrick's 2015, and Gannon 2024
are all in range.

## Package layout

```
src/guide/
  kernels.py         integrated-OU variance shape B(u) and derivatives (stable)
  ou_process.py      exact (non-Euler) OU simulation and path integration
  dynamics.py        circular-orbit SMA decay, variance scale, delta-v
  first_exit.py      closed-form first-exit-time CDF (the core result)
  monte_carlo.py     OU-driven Monte Carlo with Brownian-bridge correction
  conservatism.py    conservatism diagnostics and manuscript hooks
  config.py          Spacecraft, DensityDriver, ControlBand, GuideProblem
  stationkeeping.py  decide(): the operational decision layer
  drivers/
    rope_adapter.py  calibrate a DensityDriver from ROPE along-track output
docs/
  theory.md          derivations behind every module
  architecture.md    the ROPE -> DOPE -> GUIDE pipeline in detail
examples/            runnable scripts
tests/               full pytest suite
```

## Numerical design notes

- **Variance kernel stability.** The integrated-OU variance shape is
  `B(u) = u - 1 + exp(-u)`, which loses all significant digits below about
  `u = 1e-3` if evaluated naively. The implementation uses `expm1` to avoid the
  catastrophic cancellation.
- **Two-series CDF.** The band survival probability has an eigenfunction (theta)
  series that converges fast at large variance and a method-of-images (Gaussian)
  series that converges fast at small variance. GUIDE uses each in its fast
  regime, which keeps the CDF accurate and strictly monotone across the whole
  range (a single series loses monotonicity through truncation error in its slow
  regime).
- **Conservatism is one-sided and bounded.** The closed form treats the scalar
  process as a time-changed Brownian motion. True integrated OU has positively
  correlated increments (it is smoother than Brownian motion), so it crosses the
  band slightly later on average. The closed form therefore never materially
  under-reports exit probability and over-reports it only by a small, bounded
  margin. The Monte Carlo cross-check in the test suite asserts this direction
  explicitly.

## Validation status

The validated, TLE-driven results (KANOPUS-V 3 decay of 791 m, 61 m / 7.7
percent retrospective RMS, 148 m / 18.0 percent predictive hold-out, validated
across KANOPUS-V 3/5/IK and Swarm-A during the May 2024 G5 Gannon storm) come
from the full external pipeline, not from the self-contained examples in this
repository. The example parameters are illustrative and chosen to exercise the
operational regime.

## Manuscript transcription points

A few results are specific to the accepted manuscript (JGCD / AAS 26-664) and
are intentionally **not** inferred numerically in this code. They ship as
documented hooks that raise `NotImplementedError` with a pointer to the source
equation, so a missing transcription fails loudly instead of silently using a
guess:

- **Proposition 2 conservatism factor `c(u)`** (`guide.conservatism.
  proposition2_conservatism_factor`). Target signature: `c` increasing,
  `c(u*) = 2` at `u* ~ 2.08`, trigger band `u < 1.45`.
- **Proposition 2 crossover `u*`** (`guide.conservatism.proposition2_crossover`).
- **Delta-v cost model.** The package ships the monotone "Formula A"
  (`guide.dynamics.delta_v_for_sma_correction`), which preserves the
  monotonicity used by the Proposition 3 log-convexity argument. Reconcile
  against the exact manuscript cost model before publication.

A transparent, derivable stand-in (`guide.conservatism.transition_diagnostic`,
the ratio `u B'(u) / B(u)`) is provided for plotting until `c(u)` is
transcribed. It is explicitly not the manuscript factor.

## Putting this on GitHub

This repository is initialized with an initial commit. To push it to a new
GitHub repository, create an empty repository on GitHub (no README, license, or
gitignore, so the histories do not conflict), then run:

```bash
git remote add origin https://github.com/NexumKnight/guide-framework.git
git branch -M main
git push -u origin main
```

After the first push, the usual `git add`, `git commit`, `git push` cycle keeps
the work backed up.

## License

MIT. See `LICENSE`.

## Citation

See `CITATION.cff`, or cite the accepted AAS and JGCD papers directly once
published.
