# Валидация и диагностики Solver H1 v0.2

## 1. Правило приемки

Тест может подтвердить только: (a) continuum identity; (b) код конкретной action/formulation; (c) численную сходимость в данном domain. Он не доказывает H1, дочь, topology change или геодезическую полноту.

До каждого numerical milestone заранее сохраняются профили, domain, gauge, `Δt/Δx`, нормы, критерии и expected order. Универсальный threshold `10⁻⁶` без связи с дискретизацией не применяется. Acceptance строится из: analytic truncation order, Richardson estimate на ≥3 grids, независимого residual stencil, roundoff floor и физической чувствительности observable. Абсолютная размерная норма публикуется вместе с dimensionless нормой.

Для трёх сеток `h,h/r,h/r²` order оценивается по

\[
p_{obs}=\log_r\left|\frac{Q(h)-Q(h/r)}{Q(h/r)-Q(h/r^2)}\right|.
\]

Если последовательность не asymptotic, знаки разностей меняются или величина попала в roundoff floor, заявлять `p_obs` запрещено; результат — `not in convergence regime`.

## 2. Обязательная матрица проверок

| Проверка | Независимая ожидаемая величина | Приемка/стоп |
|---|---|---|
| 5D Minkowski + `φ=const` в истинном vacuum | `R_ABCD=0`, `M_MS=0`, constraints нулевые; все координатные gauge fields сохраняют analytic solution | Residual не растёт во времени; exact-preserving variables — до floating-point accumulation. Геометрический drift выше оценённой roundoff/discretization — FAIL. |
| 5D Tangherlini vacuum | `f=1−μ/R²`, `M=3πμ/(8G₅)`, `X=f`; horizon `R_H=√μ` | Независимо реконструированная mass constant и invariant horizon должны сходиться к analytic значениям при refinement; координатная поверхность `f=0` не считается curvature singularity. |
| Weak scalar на фиксированной Minkowski geometry | Решение `ψ=R^{3/2}φ` удовлетворяет 1+1 wave equation с centrifugal potential `3/(4R²)`; проверить regular-center branch и независимое characteristic/analytic solution | Сходимость field/flux с ожидаемым order выбранного stencil; это тест matter, не gravitating backreaction. |
| Linearized weak-gravity scalar | Первое возмущение mass/metric из constraint и perturbation Einstein equations, независимо вычислено symbolic/analytic | Численный coefficient/order совпадает в области линейного regime и меняется как предсказанный expansion parameter; если сравнение не определено из-за gauge, тест не считать пройденным. |
| тест регулярных начальных данных regular Cauchy family | Constraints construction, smooth center, finite AF charge, `Q L²<1` initial untrapped round spheres | Воспроизвести initial data и независимые H/M residuals; исходный bound применяется только к точно той massless-scalar family, не к EOS/спину. |
| Constraint propagation | Независимые Hamiltonian, momentum, GH, first-order-definition residuals | Нормы сходятся при space/time refinement с предсказанным order; refinement не приводит к secular growth выше локальной error envelope. Failure одного constraint — INVALID, даже если mass ledger closes. |
| Outgoing pulse через конечный `R_out` | `ΔM_MS + ∫F_out dt + known source/boundary terms` | Ненулевой flux должен быть независимо измерен; ledger remainder уменьшается согласно truncation order, boundary reflection отделена domain-of-dependence test. Zero flux не назначать по Dirichlet. |
| Нулевые направления/горизонт | `θ±=3 k±·∇R/R`, `θ+θ−=−9X/R²` для `k±=n±s`; independent root tracking | Сходятся sign/position/trend marginal roots при refinement и gauge comparison. Сами expansions magnitude зависит от normalization. Horizon не выводится из `M` одной общей координаты. |
| Кривизна | Trace identity и независимый Riemann reconstruction (`R`, `R_ABR^AB`, Kretschmann) | Сравнивать invariants в matched events; большие координатные derivatives не являются достаточными. Divergence требует grid/time convergence и отдельного singularity analysis. |
| Cross-gauge/cross-code | GH vs independently discretized double-null на common physical initial surface and overlap events | Сначала построить coordinate/event map; сравнить `M_MS`, `X`, scalar values, invariant expansions signs. Несопоставимые initial surfaces → `NOT_COMPARABLE`. |

## 3. Размерные и относительные нормы

Для Hamiltonian constraint `H`, dimension `L⁻²`, и momentum vector residual `M_I`, также dimension `L⁻²` в выбранной convention, сообщать `L∞`, proper-volume weighted `L2`, и norm, нормированную на локальный scale

\[
S_H=|{}^{(4)}R|+|K|²+|K_{IJ}K^{IJ}|+2κ₅²|E|+2|Λ₅|+L_{char}^{-2},
\]

где `L_char` выбирается до run из initial profile width/curvature scale. Не убирать cancellation точки, печатать absolute numerator и denominator. Momentum нормировать на соответствующие отдельные connection/curvature source scales и `L_char⁻²`.

Mass ledger для `M0≠0`: `ε_M=|ΔM+F_int−S_int|/|M0|`; если `M0` сравним с numerical uncertainty, использовать фиксированный физический `M_char` из исходных данных и также absolute mass error. Flux sign convention объявить (`F_out>0` уносит массу наружу). Вычислять flux интеграл независимой формулой через stress/current и boundary normal, не переиспользовать mass RHS.

Приемка каждой величины: (1) dimensional correctness; (2) expected order in space and time; (3) independent formulation agrees within combined Richardson/error estimate; (4) boundary/gauge/domain systematics ниже заявленного observable precision. Если сигнал меньше оценки ошибки — статус `UNRESOLVED`, не pass.

## 4. Локальные и глобальные verdicts

Следующие статусы не объединять:

- `FIELD_DOMAIN`: скаляр сформировал локальный профиль/уровневую поверхность.
- `CONSTRAINTS_VALID`: независимые ограничения сходятся.
- `TRAPPED/MARGINAL/ANTI_TRAPPED`: локальные знаки двух expansions на выбранных срезах.
- `APPARENT_HORIZON_TRACKED`: marginal tube в конкретной slicing/gauge.
- `REGULAR_EXTENSION`: построено решение через ранее ограниченную chart/gauge и проверены инварианты.
- `GLOBAL_CAUSAL_SEPARATION`: требуется полный causal construction, specified observers/asymptotics и future knowledge.
- `DAUGHTER_UNIVERSE`: требует также физически интерпретируемой независимой geometry/worldline/causal relation; не synonym for trapped region.

Локальное расширение не подразумевает horizon bifurcation; bifurcation не доказывает regular continuation; regular finite-time continuation не доказывает geodesic completeness. Любая конечная численная эволюция оставляет global verdict `UNKNOWN`, пока отдельное доказательство не закрывает его.

## 5. Старые benchmarks и источники

5D Tangherlini/Minkowski и scalar linear tests — первые прямые benchmarks. Choptuik (1993) относится к 4D massless scalar critical collapse; его exponent нельзя переносить в 5D и нельзя объявлять reproduced без threshold bisection, scaling window, AMR/resolution study. Garfinkle double-null — formulation reference, но старый локальный отчёт о double-null проверке solver имел крупные coupled null-constraint/ledger defects и не считается oracle.

Исторические #40–41F/отчёт о тестах коллапса–121 выявили finite-grid constraint/mass-budget несогласованность, включая #41F incompatible discretized mass/constraint identities. Эти цифры — предупреждающий baseline, не acceptance thresholds и не ground truth. Старый solver не копировать автоматически.
