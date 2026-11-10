# Experiment 41B — forensic solver audit and attempted repair

## Verdict first

The mass-budget problem is **not repaired**. The tested constraint-reintegrated solver loses quasi-local mass even though its imposed outer scalar flux vanishes. A second, mathematically derived evolution form advances the Misner–Sharp mass with the Einstein flux equation and closes the outer ledger, but the independently measured Hamiltonian constraint develops a persistent residual. That formulation is also **not accepted as physically valid**. No daughter-universe or horizon-formation inference is made.

What did improve: #40's mass-flux residual had compared one full RK4 mass increment against flux evaluated only at the final stage. It now compares the constraint-based mass increment against the RK4-weighted flux from all four stages. This makes the diagnostic time-integrator consistent. It changes the reported intermediate max flux residual from `1.76018e−3` to `1.66388e−3` and the strong value from `.151851` to `.150390`; it does **not** cure the `−2.1711%` mass drift.

## Independent field-equation audit

The local model uses `c=1`, `κ=8πG`, action

`S=∫√−g [R/(2κ) − (∇φ)²/2 − V(φ)] d⁴x`,

and the canonical stress tensor. In polar-areal, zero-shift coordinates

`ds²=−α²dt²+A²dr²+r²dΩ²`, `f=A⁻²=1−2m/r`,
`Π=(A/α)φ_t`, `Φ=φ_r`,

the orthonormal-normal energy density and radial pressure are

`ρ=f(Π²+Φ²)/2+V`, `p_r=f(Π²+Φ²)/2−V`.

The Einstein constraints and scalar equations are

```
m_r = (κ/2) r² ρ
(ln α)_r = [m/r² + (κ/2) r p_r] / f
φ_t = α Π/A
Π_t = r⁻² ∂r(r² α Φ/A) − α A V_φ
```

With `K_ij=−(1/2)L_n γ_ij`, the radial momentum relation and Misner–Sharp mass flux are

```
K^r_r = −(κ/2) r Π Φ/A
m_t = (κ/2) r² α f^(3/2) Π Φ.
```

Here `m(r)` is the geometric Misner–Sharp mass length `G M_MS/c²`. At a finite `r_max`, it is a quasi-local mass. ADM mass is the asymptotic limit `lim(r→∞)m(r)` for asymptotically flat data; `m(r_max)` approximates it only when the exterior tail is negligible and boundary effects are controlled. Parameters use `μr`, `μt`, `φ/M_Pl` and `V/(μ²M_Pl²)`; no physical value of `μ` is chosen.

The equations above are mutually consistent in the continuum: differentiating `m_r` in time, substituting the Klein–Gordon equation and using `f_r` plus the lapse equation gives `∂r(m_t)` with the stated flux. The sign/factors of the source term match the code convention `Π=A n(φ)`. The potential derivative also matches `V=λφ²(φ−1)²/4+δ(3φ²−2φ³)`.

At the regular center, `m(0)=0` and a smooth spherical scalar has even parity (`Φ(0)=0`); the implementation uses the regular spherical divergence limit. The outer condition currently fixes `φ=Π=0` and normalizes `α(r_max)=1`. Fixing both scalar variables is not an outgoing radiation condition; if a pulse reaches the boundary, this boundary can reflect it and should be replaced by a tested incoming-characteristic/Sommerfeld treatment. In this tested baseline, the outgoing null ray launched at `r=2` reaches only `r=4.0307` by `μt=2`, far from `r_max=16`; the tested ray therefore does not reach that boundary.

## Forensic reproduction of the mass remainder

Original intermediate case: `A=1`, `R0=2.8`, `κ=.1`, `r_max=16`, `dr=.1`, `dt=.0125`, `μt_end=2`.

```
M_MS(0)       = 0.0970447857396 μ⁻¹
M_MS(final)   = 0.0949378821303 μ⁻¹
ΔM            = −0.00210690360926 μ⁻¹
ΔM/M(0)       = −2.17106318%
max |F_boundary| = 0
∫ F_boundary dt  = 0
ledger remainder = −0.00210690360926 μ⁻¹ (−2.17106318%)
```

The boundary flux is evaluated from the Einstein flux expression and the saved boundary values; it is zero because `Π(r_max)=0`. The finite-radius outer mass changes anyway. This is a numerical constraint/evolution inconsistency, not a conservation law for the scalar matter alone and not physical leakage at the outer boundary.

Sensitivity checks were performed without changing the initial physical profile:

| Check | settings | `ΔM/M(0)` |
|---|---|---:|
| baseline | `dr=.1`, `dt=.0125`, `r_max=16` | −2.17106318% |
| half timestep | `dr=.1`, `dt=.00625` | −2.17106319% |
| quarter timestep | `dr=.1`, `dt=.003125` | −2.17106319% |
| larger domain | `r_max=24` | −2.17422464% |
| larger domain | `r_max=32` | −2.17578166% |
| spatial refinement | `dr=.2,.1,.05`, `dt=dr/8` | −2.15376778%, −2.17106318%, −2.17885401% |

The near-invariance under `dt` reduction argues against RK4 temporal error as the main cause. Expanding `r_max` does not cure the drift, and no outgoing ray reaches the boundary in the baseline interval, so physical boundary leakage/reflection is not its cause here. Spatial refinement trends toward a nonzero remainder rather than resolving it away. The local least-squares fit of measured `m_t` to `r² α f^(3/2)ΠΦ` over active bulk cells gives coefficient `.04913` (continuum expectation `κ/2=.05`) and correlation `.99858`; this supports the local flux factor but does not close the global ledger. The remaining error is a failure of the implemented discrete constraint/evolution compatibility. The tests localize it to the bulk numerical formulation, but do not yet isolate a single offending operator/product term. Claiming a more specific root cause would exceed the evidence.

## RK4-consistent flux residual

For the constraint-reintegrated baseline, the scalar and metric fields are evaluated at RK4 stages `1,2,3,4`. The flux diagnostic is now

`R_m = max_r |m_new−m_old − Δt(F₁+2F₂+2F₃+F₄)/6| / Δt`.

This is preferable to comparing a full-step mass change to the endpoint-only flux. It remains a residual between a radially re-solved constraint and an evolution equation; it is not an independently evolved mass. New values are: weak `7.4246e−5`, intermediate `1.66388e−3`, strong `1.50390e−1`. Constraint-reintegrated intermediate/strong values remain invalid for a high-precision formation claim.

## Alternative mass-equation evolution: test, not accepted repair

To test the local cause without hand-editing the mass, a separate path initializes `m(r)` from the Hamiltonian constraint and then evolves it with the Einstein equation `m_t=(κ/2)r²αf^(3/2)ΠΦ` at all RK4 stages. The radial lapse is recomputed, while `m_r−(κ/2)r²ρ` is independently monitored. With a zero-flux fixed boundary this path preserves outer `m` and satisfies its same-stage RK ledger to round-off. That ledger is a scheme identity because the same stage fluxes enter the update, **not** independent evidence of correct physics.

| Profile | `dr` | mass change | max absolute Hamiltonian residual | relative H L2 | max compactness |
|---|---:|---:|---:|---:|---:|
| zero | .1 | 0 | 0 | 0 | 0 |
| weak `A=.35` | .1 | 0 | `5.622e−5` | `8.190e−3` | `.008960` |
| intermediate `A=1` | .1 | 0 | `1.554e−3` | `3.483e−2` | `.045271` |
| strong `A=1.8,R=1,κ=.5` | .025 | 0 | `3.094e−2` | `4.974e−2` | `.432443` |
| intermediate grid | .2 | 0 | `1.846e−3` | `4.304e−2` | `.045461` |
| intermediate grid | .05 | 0 | `1.475e−3` | `3.487e−2` | `.045221` |
| intermediate grid | .025 | 0 | `1.456e−3` | `3.504e−2` | `.045209` |

The nonzero Hamiltonian residual approaches roughly `1.45e−3` under this modest refinement. The strong run reaches `max(2m/r)=.432443<1`, but its Hamiltonian residual is too large; it is neither evidence of a horizon nor a reliable collapse result. The momentum residual in this alternate path is tied algebraically to the mass-flux-derived `K` and is not an independent evolution test.

This comparison establishes the problem sharply: the original formulation maintains the radial Hamiltonian constraint at each slice but fails mass balance; evolving `m` by its Einstein equation conserves that boundary ledger but fails the radial Hamiltonian constraint. A solver satisfying both requires a compatible discretization of the scalar stress-energy flux, radial constraint, lapse and regular-center operators, preferably an SBP/finite-volume formulation with a proven discrete summation-by-parts identity. The current two formulations are diagnostic endpoints, not a selection of whichever result looks favorable.

## Before/after status

| Quantity | original #41 baseline | #41B change/attempt | accepted conclusion |
|---|---:|---:|---|
| intermediate `ΔM/M` | −2.1711% | unchanged in constraint form; zero in mass-equation trial | mass budget unresolved because trial H constraint fails |
| outer integrated flux | 0 | 0 | no physical boundary leakage measured |
| intermediate RK flux residual | `1.76018e−3` endpoint comparison | `1.66388e−3` RK4-stage comparison | diagnostic improved; bulk mismatch persists |
| intermediate relative H L2 | `4.61e−3` | original unchanged; mass-evolved trial `3.48e−2` | neither formulation passes strict accuracy |
| intermediate relative momentum L2 | `5.76e−2` | original unchanged | constraints remain inadequate |
| strong max `2m/r` | `.4324` | no trustworthy horizon; mass-evolved H residual `.03094` | strong trajectory INVALID for physical interpretation |
| H1 parent/daughter time | NOT_COMPARABLE | unchanged | no `a_d`, no linked parent observer |

## Tests and next implementation requirement

The unit tests assert the finite-boundary mass ledger, zero/weak/intermediate/strong cases, RK4 stage accounting, convergence trend of the alternate Hamiltonian residual, and the invalidity of the strong case. Run locally:

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
uv run --no-sync python -m pytest -q
uv run --no-sync python -m compileall -q src tests
```

No horizon-penetrating gauge is implemented. A proper repair next needs: (1) SBP or finite-volume radial operators for `φ`, flux divergence and constraints; (2) exact regular-center parity rows; (3) a constraint-propagation test over at least `dr=.2,.1,.05,.025` and independent CFL refinement; (4) incoming-characteristic outer boundary and boundary-energy ledger; (5) independent matter-energy and Misner–Sharp mass checks. Do not resume strong-field causal-formation claims until both the Hamiltonian constraint and mass ledger converge together.

The simulation time is `μt`. It is not automatically daughter FLRW proper time. Without a physically defined `a_d(t_d)` and parent observer, `dt_d/dT_p=K0a_d` remains **NOT_COMPARABLE**. No numerical parent seconds are assigned.
