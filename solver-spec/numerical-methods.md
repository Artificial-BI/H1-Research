# Численные методы и выбранный CPU baseline

## 1. Решение для первой реализации

**Уточнение после GH-IBVP аудита:** baseline должен пониматься как full 5D Cartesian GH system restricted to SO(4)-equivariant fields, не как наивный 2D GH только для `h_ab`. В polar/warped записи full harmonic constraint содержит `−3∇^aR/R` и angular connection terms; см. [GH-IBVP/CPBC study](gh-ibvp-cpbc-proof.md). Условие `Γ[h]^a−H^a=0` в [уравнениях §3](mathematical-equations.md) — provisional orbit-gauge placeholder, не завершённый full-5D GH gauge. Gate A сейчас **BLOCKED**: reduced metric incoming CPBC map и center parity/SBP closure не доказаны.

Кандидатный baseline: SO(4)-reduced 5D Einstein–scalar в 1+1 generalized-harmonic (GH) gauge, метод линий, high-order finite differences с summation-by-parts (SBP) spatial derivative и penalty/SAT boundary treatment, explicit RK4 для начальных гладких/умеренных задач. Double-null — **отдельная независимая перекрёстная постановка** для ограниченного overlap-domain, не второй режим в общем массиве координат. До реализации требуется proof/verification strong hyperbolicity и constraint-preserving boundary system для фактически reduced GH equations (см. [уравнения](mathematical-equations.md)). Без этой проверки Gate A — STOP.

Почему не «сразу strongest code»: старые polar-areal и coupled-FV попытки имеют одновременно неудовлетворительные независимые constraints и mass budget; double-null демонстрация также partial. Их код и результаты — архивные предупреждения, не oracle. Новое ядро должно заново реализовать continuum system после независимой derivation.

## 2. Формулировки пространства-времени

| Формулировка | Что даёт | Ограничение/решение |
|---|---|---|
| Polar-areal, `R=r`, lapse radial ODE | Простая масса/constraint integration и хороший слабополевой контроль | Coordinate breakdown при `f=(∇R)²→0`; только pre-horizon reference, нельзя продолжать через marginal surface подменой `sqrt(f)`. |
| SO(4) generalized harmonic quotient | Cauchy data, динамический `R`, потенциально horizon-penetrating coordinates, волновая структура | Нужно вывести reduced equations/gauge constraints/boundary characteristics и проверить strong hyperbolicity; это выбранный основной кандидат. |
| Double-null | Явные radial null directions, полезно для trapped horizons и causal patches | Характеристические data вместо исходных тест регулярных начальных данных Cauchy data; caustics/center patch, residual null-coordinate gauge, finite domain; независимая схема и mapping initial surfaces обязательны. |
| ADM/BSSN/Z4c radial reduction | Constraint monitoring/damping; знакомые 3+1/4+1 методы | 5D extension, spherical regularization, gauge and characteristic boundaries должны быть доказаны; constraint damping не исправляет неверные equations. |

Generalized harmonic предпочтён для возможного прохождения horizon в координатах; это не гарантия, что конкретный `H^a` не выродится. При crossing failure фиксируется `GAUGE_STOP`. Excision разрешается только если все physical/numerical characteristics на excision boundary направлены внутрь удаляемой области и внешняя solution не зависит от неё; иначе не применять.

## 3. Пространственная и временная дискретизация

### Spatial

1. Начальный baseline — равномерная radial сетка, double precision; center parity задаётся явно через гладкое Cartesian extension/Taylor regularity.
2. SBP derivative должен удовлетворять дискретной summation-by-parts identity с совместимой norm matrix. SAT penalty parameters выводятся из characteristic fields/energy estimate этой системы, а не подбираются визуально.
3. Для ограниченной scalar equation допускается conservative face flux `R³`-weighted only после точного вывода из reduced divergence `R^{-3}∂_a(R³∂^aφ)`. Геометрические source terms и массо-геометрия не должны быть «поглощены» в несогласованный FV flux.
4. Все evolution RHS (включая gauge/damping) пересчитываются на каждой RK stage из одного stage-state. Constraint projection, если когда-либо понадобится, должна быть отдельным elliptic solve с измеренным projection defect; запрещён post-hoc mass correction.

### Time/CFL

Для явного RK baseline шаг удовлетворяет `Δt ≤ C_CFL Δx/max|λ_char|`, где скорости извлекаются из principal symbol фактической gauge-reduced системы. `C_CFL` задаётся stability-region/linear eigenvalue analysis выбранной spatial operator; нельзя назначать его по старому solver. Каждая benchmark серия повторяется как минимум при двух CFL и трёх последовательно измельчённых `Δx`; изменение результата только от CFL означает time error, не spatial convergence. RK4 stage order и boundary values должны быть согласованы на всех четырёх стадиях.

Если near-singular gradients/constraints делают explicit step impractical, implicit/IMEX схема — новый backend с собственными stability/consistency tests, а не незаметная замена.

### Constraint control

Эволюционируются constraints как диагностики, не заменяются тождествами из той же RHS. GH constraint damping допустим лишь при derived continuum constraint propagation system и проверке, что damping не изменяет физический solution/energy budget в пределе refinement. Отдельно считать Hamiltonian, radial momentum, GH gauge constraint и first-order definition residuals независимыми stencils/quadratures.

## 4. Boundary, center, horizon и mass budget

- **Center:** parity conditions всех tensor/scalar components; analytic limits `R→0`; никаких arbitrary floors в `1/R` terms.
- **Outer boundary:** incoming characteristic data физически заданы. Вакуумная AF граница/падающий background объявлены в run manifest. Outgoing pulse не должен быть зажат искусственным Dirichlet `φ=Π=0`. Если absorbing condition не доказана, domain of dependence урезается до времени первого отражённого сигнала.
- **Finite radius:** отдельно логировать `ΔM_MS(R_out)`, независимо вычисленный scalar/gravitational energy flux, интегральный flux через boundary и budget remainder. Равенство flux с mass change, если оно следует из выбранного continuum mass identity, дискретно оценивается как convergence observable; exact ledger closure by construction не является независимой GR validation.
- **Horizons:** marginal sphere — root `X=(∇R)²=0`/одной `θ±=0`, tracked by bracketed root on each slice. Trapped/anti-trapped classification uses both signs. Apparent horizon depends on slicing; event horizon not found by finite-time root tracking.
- **AMR/excision:** отключены в MVP. Прежде чем включать, провести refinement-interface conservation/constraint tests, moving-horizon crossing test, reflux/characteristic interface proof. Вывод за excision не называется global solution.

## 5. Независимая cross-check formulation

Для selected analytic/weak families строится double-null comparator из covariant 1+1 equations, отдельно дискретизированный (например, characteristic diamond/finite differences). Он не вызывает GH RHS, constraints, mass flux helpers или interpolation branch solver. Сравнивать разрешено gauge-invariant `M_MS`, `X`, scalar along matched geometric events and expansions after constructing a coordinate map on overlap. Raw `(t,x)` массивы и `(u,v)` не сравниваются напрямую. Если common initial hypersurface/event matching не построено, результат помечается `NOT_COMPARABLE`.

Double-null отчёт о double-null проверке в архиве не считается независимо validated strong-field reference: его own coupled null constraint/ledger residuals оставались заметными. Его equations/case служат только regressions/implementation lesson after re-derivation.

## 6. Сравнение инструментов и оборудования

Для 1D reduced CPU baseline Python 3.12 + NumPy может быть приемлем для прототипа на 7950X/64 GiB, но runtime ещё не измерен. SymPy пригоден для локальной symbolic algebra проверки, не для time evolution. Python overhead не следует оптимизировать до профилирования. C++/Numba/CUDA не обязательны; доступность текущих compiler/toolkits не подтверждалась. RTX 3090 не нужна для single-radial baseline; GPU может усложнить reproducibility и не чинит physics.

GRChombo, GRTeclyn, Einstein Toolkit и LEAN рассматриваются только как возможные references по formulations/AMR. Их текущие версии, Windows 10 buildability, D-dimensional support/license/runtime здесь **не проверялись**; ничего не скачивалось. EinsteinPy/SymPy не являются заменой nonlinear GR evolution solver. Небольшое собственное ядро оправдано лишь если оно останется узкой 1+1 reference implementation, а формулировка до этого независимо проверена; framework fork/port — отдельное решение после обзорного build test, не часть первого issue.

## 7. Численная достоверность и порядок сходимости

Для quantity `Q(h)` на трёх геометрически связанных resolution `h,h/r,h/r²` вычислить observed order

\[
p_{obs}=\log_r\left|\frac{Q(h)-Q(h/r)}{Q(h/r)-Q(h/r^2)}\right|.
\]

Использовать только если разности одного знака, не dominated by roundoff и времена/границы согласованы; иначе сообщать `p_obs undefined`, не усреднять в красивое значение. Отдельно варьировать time step, `R_out` и boundary distance. Не объявлять order по двум сеткам.

Нормы constraint остатка нормировать dimensionally: `H` имеет `L⁻²`, momentum constraint `L⁻²` для выбранных нормальных проекций; делить на сумму локальных scale terms с floor от заранее фиксированной characteristic length `L_char` из входного профиля (`L_char⁻²`), отдельно публикуя абсолютные `L∞,L2` и relative norm. Не выбирать scale после просмотра результата. Mass ledger делить на `M_initial` лишь если он надёжно ненулевой; иначе приводить абсолютную ошибку и физический reference mass. Acceptance — предсказанный refinement trend и error budget, см. [validation](validation-and-diagnostics.md), не универсальный epsilon.
