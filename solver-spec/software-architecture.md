# Архитектура будущего самостоятельного Solver H1 v0.2

## 1. Принцип разделения

Будущий репозиторий отдельный от архивного H1 Simulation Lab. Архивный `src/`, UI, solver data и старые experiment reports не импортировать как production backend без line-by-line verification. Пакеты физической модели, continuum formulation, численного backend и интерфейса должны быть независимо версионируемы. Переключение action меняет manifest и backend eligibility; нельзя менять одну формулу в UI, оставляя прежний source/action.

## 2. Модули и контракты

| Компонент | Контракт и вход/выход |
|---|---|
| `physics_model` | Immutable `ActionSpec`: `action_id`, dimension/signature, units, `κ₅²`, `Λ₅`, potential identifier + parameters, matter content, derived equations revision. Только физические константы/функции, без grid/UI. |
| `geometry` | Metric/gauge state `h_ab,R` (и lapse/shift representation), coordinate chart, inverse metric, derivatives/characteristics; outputs `X`, curvature tensors/invariants and null expansions through independently testable APIs. |
| `matter` | `φ,Π` or first-order scalar state, potential/stress/flux. Must report dimensions and validity domain; stress obeys covariant conservation when matter equations hold. |
| `initial_data` | Input profiles + gauge choices → `γ_IJ,K_IJ,φ,Π`; solve constraints or reject; output residual fields, asymptotic charge, regularity checks and provenance. No silently “repairing” arbitrary data. |
| `evolution_formulation` | Continuum reduced equations, gauge sources, characteristic fields, constraints/propagation identities, and boundary characteristic contract. Fixed by `formulation_id` and derivation revision. |
| `numerical_backend` | Discretization-specific operators/time integrator/grid; consumes a typed state and `ActionSpec+FormulationSpec`; returns stage/step diagnostics. Must not contain interpretation labels such as “daughter formed.” |
| `boundaries` | Center parity, outer incoming/outgoing characteristic data, finite-domain extent and flux quadrature. Geometry/matter boundary conditions explicit and separately versioned. |
| `constraints` | Independent H/M/GH/definition residual reconstruction; norm definitions/scales and refinement metadata. Not an alias to RHS or projection update. |
| `diagnostics` | `M_MS`, ADM/asymptotic mass when valid, flux ledger, `θ±`, `X`, curvature invariants, geodesic/null rays and causal-domain flags. Each diagnostic carries domain/gauge dependence. |
| `experiment_manager` | Validates typed parameter sets, creates immutable run manifest, executes prescribed matrix, collects status values `PASS/FAIL/UNKNOWN/NOT_COMPUTED` with reason and evidence. No physics inference from thresholds alone. |
| `persistence` | Versioned manifest + arrays + tabular diagnostics; checks units/action/formulation compatibility on load. Never silently migrates scientific data. |
| `visualization` | Read-only view of persisted/computed arrays; axes and units explicit; overlays never feed back into physics. “Schematic” overlays watermarked. |
| `extensions` | Separate plug-in actions only when complete action, constraints, evolution and boundary/junction equations are provided. No runtime mix of Einstein, ECSK and EGB partial terms. |

Typed state objects should carry coordinate chart, units, grid locations (cell/node/face), tensor index placement, time stage and validity mask. Arrays with different gauges/resolutions are not interchangeable even if shapes match.

## 3. Versioning and compatibility

Required immutable identifiers:

- `action_id` + semantic version and exact parameter units;
- `reduction_id` (SO(4), dimension and orbit convention);
- `formulation_id` (GH or double-null; metric/gauge conventions);
- `backend_id` (operators, stencil order, integrator, damping, arithmetic);
- `initial_data_id` + profile parameters and constraint solve version;
- `boundary_id` and outer radius/characteristic treatment;
- `diagnostic_id` (mass convention, expansions normalization, norm scales).

Results with mismatching action/reduction/formulation are separate experiments, not comparable curves. Schema changes increment major version if units, tensor conventions or physics meaning changes. Loading an unknown major version fails closed; no silent defaults.

## 4. Run artifact contract

Portable result package (proposed, not implemented):

1. `run.json`: schema version, UTC timestamp, code revision/hash, Python/library versions, OS/CPU (GPU if used), action/formulation/backend IDs, constants and units, initial-data spec, boundaries, grid/CFL/time range, output cadence, random seeds if stochastic, test/acceptance profile, status and limitations.
2. `fields.npz` (or future chunked equivalent): numeric arrays with explicit coordinate/time arrays, SI/natural units and component descriptions in manifest; no pickle/object arrays.
3. `diagnostics.csv`: machine-readable time/radius, mass, boundary flux, constraints, `X,θ±`, curvature norms and masks; CSV header states unit per field.
4. `report.md`: generated summary referencing exact files and status; never contains implicit “confirmed H1.”

Exact byte reproducibility may vary by library/OS; goal is deterministic inputs and numerical results within documented convergence/error bounds. For very large runs, replace NPZ only after storage-size/streaming tests; HDF5 is optional and would add a dependency decision.

## 5. UI/CLI isolation

First implementation issue should preferably produce a headless library/CLI and deterministic small validation runs before UI. A future GUI calls experiment-manager contracts; it must not duplicate equations or permit unversioned formula changes. Local Windows launcher/environment should be pinned after dependency feasibility review. No server/cloud service is required.
