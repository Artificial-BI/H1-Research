# Дорожная карта подготовки и реализации

Каждый этап — отдельное решение GO/STOP. Выполнение этой дорожной карты не начато: текущий deliverable состоит только из спецификаций.

## Stage A — аналитическая постановка и well-posed IVP

**Deliverables:** подписанная action convention; полная SO(4) reduction; GH (или альтернативная) evolution/gauge source; first-order variables; principal symbol/strong hyperbolicity proof or cited derivation valid for exact reduction; constraint propagation; origin/outer characteristic boundary system; initial-data free functions and constraint solve; mass and null diagnostics derivation.

**GO gate:** символически восстановлены Minkowski и Tangherlini; analytic reduced equations совпадают с independent variation; boundary characteristic count корректен; все units, signs and definitions заморожены.

**STOP:** principal symbol не strong hyperbolic; constraints overdetermined; boundary conditions inject uncontrolled incoming modes; regular center not covered. Не кодировать evolution kernel.

## Stage B — минимальное CPU evolution ядро

**Deliverables:** фиксированная равномерная сетка, GH/SBP-MOL/RK baseline, constraints/diagnostics independent, no AMR/GPU, immutable run manifest, Minkowski/Tangherlini/linear scalar/тест регулярных начальных данных tests и reproducible artifacts.

**GO gate:** 3+ spatial resolutions, 2 CFL, independent residuals converge at predicted order; mass flux ledger converges; origin/boundary and exact solutions pass predeclared error budget.

**STOP:** no convergence regime, constraints worsen, budget mismatch stagnates, result depends on hidden clipping/projection, or baseline time coordinate/gauge stops. No physical collapse claims.

## Stage C — nonlinear collapse and horizons

**Deliverables:** pre-registered smooth finite-energy profiles (zero, weak, intermediate, strong), parameter scan kept bounded; apparent horizon tracking via both expansions; null rays, curvature invariants, mass/flux ledger, gauge/boundary failure status.

**GO gate:** constraint/mass criteria remain met over horizon formation or explicitly valid horizon-penetrating chart; outcome robust to resolution/CFL/domain changes; at least one independent cross-gauge check on common events.

**STOP:** `GAUGE_STOP`, constraints invalid, horizon at outer reflection boundary, curvature not resolved, or no comparable independent data. Record as unresolved, not absence of physics.

## Stage D — high-density extension

Only after A–C: choose one separately complete theory (ECSK with state-derived stress, EGB, or another action), derive new constraints/evolution/boundary terms, test stability/hyperbolicity/low-energy recovery and solve new initial constraints. No additive “bounce correction” to Einstein equations without an action. анализ микроскопического спинового состояния has not supplied a computed spin-state closure.

**GO gate:** physical state and action define source sign/scaling; ghost/gradient/strong-coupling and cutoff analysis; independent tests. Otherwise extension is rejected/conditional.

## Stage E — common gravity, branes and daughter formation

Separate theoretical project: specify one bulk/brane or other common sector action; establish localization, junction conditions, regular formation initial data, causal cones/event matching, parent exterior flux, induced daughter metric and global extension. Only after this may independent clocks be compared. A trapped region or FLRW patch is insufficient.

**GO gate:** explicit solution/controlled numerical spacetime with constraints, regularity, null reachability, asymptotics and independent global causal audit. No result from Stage C automatically passes Stage E.

## Обязательный stopping rule

На любом этапе FAIL означает исправление формулировки/математики до следующего этапа, не подбор параметров под желаемую картину. Изменение action, gauge, boundary condition or status vocabulary создаёт новую versioned experiment.
