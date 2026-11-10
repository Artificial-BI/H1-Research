# Experiment 41D — FV prototype for coupled Einstein–scalar evolution

**Status: PARTIAL ONLY — not a repaired or validated coupled solver.** This note records an opt-in cell-centred finite-volume prototype and its failed acceptance checks. It does not replace the historical #40/#41 calculations, and it does not compute a horizon-penetrating evolution, a daughter universe, or a parent–daughter clock relation.

## Model and continuum equations

The prototype keeps the #40 action and potential, in units where the code's length/time variable is dimensionless `mu` and `c=1`:

`S = integral sqrt(-g) [R/(2 kappa) - (nabla phi)^2/2 - V(phi)] d^4x`, with `kappa=8 pi G`,

`V(phi)=lambda phi^2 (phi-1)^2/4 + delta (3 phi^2-2 phi^3)`.

In polar-areal coordinates, `ds^2=-alpha^2 dt^2+f^-1 dr^2+r^2 dOmega^2`, `f=1-2m/r`, define `Pi=phi_t/(alpha sqrt(f))` and `Phi=partial_r phi`. Then

* `rho = f(Pi^2+Phi^2)/2 + V`, `p_r=f(Pi^2+Phi^2)/2-V`;
* `m_r=(kappa/2) r^2 rho`;
* `partial_r ln(alpha)=[m/r^2+(kappa/2) r p_r]/f`;
* `phi_t=alpha sqrt(f) Pi`;
* `Pi_t=r^-2 partial_r(r^2 alpha sqrt(f) Phi)-alpha/sqrt(f) V_phi`;
* `m_t=(kappa/2) r^2 alpha f^(3/2) Pi Phi`.

The last identity is the polar-areal Misner–Sharp mass balance. With `M(r,t)=m(r,t)`, its finite-domain ledger is `M(R,t2)-M(R,t1)=integral(m_t(R,t) dt)`. For an outgoing wave, the boundary derivative satisfies `Phi(R)=-Pi(R)` and the expected mass flux is negative. The finite-radius quantity is a quasi-local Misner–Sharp mass; it is an ADM mass only in the asymptotically flat limit with suitable falloff.

## Discrete prototype

`src/h1lab/research/experiments/coupled_fv.py` stores cell-centred `phi,Pi` and face-centred `m,f`. Shell coordinate volumes are `DeltaV_i=(r_(i+1/2)^3-r_(i-1/2)^3)/3`; a midpoint shell-average of the Hamiltonian source gives an algebraically solved recurrence for `m_(i+1/2)`, using the same shell average in `f=1-2m/r`. The lapse is integrated inward from `alpha(R)=1`. Scalar flux `alpha sqrt(f) Phi` is face-centred and differenced as a spherical flux divergence. The outgoing outer characteristic is applied as a face value, not by zeroing the last cell. Classical RK4 recomputes constraints at all four stages.

The Hamiltonian defect is independently reconstructed with two-point Gauss integration of a linearly reconstructed scalar and interpolated mass; it is not the midpoint recurrence residual. The time-differenced mass-flux and momentum checks are likewise independent reconstructions. No mass correction is applied. `f<=1e-5` stops the polar-areal chart; no continuation through a trapped surface is attempted.

This is not an SBP/finite-volume proof of the discrete Bianchi identities. Cell reconstruction, lapse quadrature, nonlinear metric factors and RK stages do not yet share an exact compatible product rule. Refinement reduces the defects in some norms, but does not make them sufficiently small or uniformly monotone. The backend therefore remains a diagnostic prototype.

## Fixed numerical evidence

Preset: `A=.35`, `R0=1.5`, `kappa=.1`, `r_max=8`, duration `mu t=8`, `dt/dr=.2`; potential defaults `lambda=1`, `delta=.02`. Global signal-normalized norms are used in the table to avoid dividing by exponentially tiny vacuum-cell sources.

| dr | dt | independent H L2 | local mass-flux L2 | momentum Linf | cumulative outer ledger / M0 | max outward power | min f |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .2 | .04 | 1.22698e-2 | 2.53507e-2 | 1.46021e-1 | 1.70841e-3 | 1.97833e-3 | .994311 |
| .1 | .02 | 3.10474e-3 | 6.76763e-3 | 4.20425e-2 | 1.31061e-4 | 1.90055e-3 | .994312 |
| .05 | .01 | 7.74751e-4 | 1.88917e-3 | 1.09777e-2 | 1.35136e-4 | 1.87118e-3 | .994313 |

Hamiltonian residual decreases by approximately four per halving of `dr` (second-order trend); local flux and momentum residuals also decrease. The integrated mass-ledger remainder improves from 0.171% to 0.0131%, then is 0.0135% at the finest listed grid: it has not demonstrated clean monotone convergence. At `dr=.1`, reducing `dt` from `.02` to `.01` to `.005` changes H L2 only `3.1047e-3, 3.0990e-3, 3.0960e-3`, flux L2 `6.7676e-3, 6.8585e-3, 6.8883e-3`, momentum Linf `4.2042e-2, 4.2257e-2, 4.2292e-2`, and ledger remainder `1.3106e-4, 1.3102e-4, 1.3101e-4`; the observed floor is spatial, not cured by smaller time steps.

A finite-boundary pulse check (`A=.35`, `R0=.5`, `R=4`, `dr=.1`, `dt=.02`, duration 4) has nonzero outward power `2.96230e-3`; outer mass changes from `.0017645480` to `.0009765178`, while the integrated balance remainder is `3.46214e-5` (1.962% of initial mass). This is a real nonzero-flux test, but the remainder is too large to call the gravitating ledger validated.

For comparison, the historical #40/#41B constraint-reintegrated scheme had a medium outer-mass drift of `-2.1711%` with modelled zero boundary flux and effective profile convergence `p≈1.508`. The #41B mass-evolution variant conserved its discrete ledger but had H residual `1.554e-3` at `dr=.1` and `1.456e-3` at `.025`. The new prototype improves the medium ledger discrepancy on its `r_max=8` refinement sequence, but its `dr=.1` independent H L2 is `3.10e-3` and momentum Linf `4.20e-2`; it does **not** satisfy both constraints and mass balance to a demonstrated reliable error budget. These runs are not a valid strong-field calculation: `min f≈.9943`, far from a horizon.

Exact Minkowski control (`phi=Pi=0`, `V(0)=0`) yields `m=0`, `alpha=1`, `f=1`, and zero residuals. It is a code sanity control only. The UI's optional coarse trial uses the same fixed preset at `dr=.2`; it labels the outcome PARTIAL ONLY and does not alter legacy panels or results.

## Scope and next numerical requirement

The code demonstrates a gravitating constrained radial FV *attempt* with matter evolving on a metric recomputed from the Hamiltonian/lapse equations. It does not yet provide a mutually compatible semidiscrete Einstein–matter system. The next solver step must derive a shared face flux/source discretization (or SBP identity) for the discrete Hamiltonian, momentum and mass-balance laws, including the radial geometric sources, and then verify a nonzero outgoing boundary wave on at least three grids with fixed and varied CFL. Any projection must be separately budgeted rather than hiding its mass defect.

The polar-areal chart never crosses `f=0`; no apparent/event horizon or global causal region is established here. `mu t` is the solver coordinate, not a daughter FLRW cosmic time. There is no computed `a_d` and no physically identified parent observer, so `dt_d/dT_p=K0 a_d` is **NOT_COMPARABLE** and no parent seconds are reported. No claim about H1 formation follows from the prototype.

Reproduce with `.venv[локальный путь исключён] -m pytest tests/test_coupled_fv.py` or `uv run pytest tests/test_coupled_fv.py` when the local uv cache is writable. Full project commands are recorded in the task report; no new dependencies are used.
