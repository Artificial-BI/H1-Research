# Физическая спецификация Solver H1 v0.2

## 1. Назначение и уровень утверждений

Solver H1 v0.2 — будущий исследовательский код для проверки ограниченной классической модели; он не должен кодировать желаемый ответ гипотезы H1. Исторические результаты маркируются `verified in original model`, `conditional`, `not derived`; passing test — численная проверка реализации заявленной системы, не подтверждение физической реальности модели.

Рабочая H1-картинка допускает отдельные local matter/metric/proper time уровней и возможный общий гравитационный сектор. Ни дочерняя область, ни общая связь между manifold, ни физический закон `dt_D/dT_P` не выведены. H2 (тёмная материя/ускорение), повторные рождения, quantum tunneling и интерпретация black-hole interior как другой universe вне первой модели.

## 2. Выбор первой модели

Модель первой реализации — **одно 5D Lorentzian manifold**, классическая GR с одним canonical real scalar field. Сигнатура `(-++++)`, `c=ℏ=1`; независимые поля — bulk metric `g_AB` и `φ`. Объёмная spin/torsion, branes, second metric, extra-dimensional mediator и высококривизные поправки отсутствуют.

\[
S=\frac1{2\kappa_5^2}\int_{M_5}d^5x\sqrt{-g}
(R-2\Lambda_5)
-\int_{M_5}d^5x\sqrt{-g}\left[\frac12g^{AB}\nabla_A\phi\nabla_B\phi+V(\phi)\right]
+S_{boundary}.
\]

Где физическая граница/cut реально введена, требуется согласованный Gibbons–Hawking–York/асимптотический член и заданная boundary condition. В базовой постановке нет thin shell/материальной границы. `κ₅²=8πG₅`, `[κ₅²]=L³`, `[Λ₅]=L⁻²`, `[φ]=L⁻³ᐟ²`, `[V]=L⁻⁵`.

Вариация даёт

\[
G_{AB}+\Lambda_5g_{AB}=\kappa_5^2T_{AB},\quad
T_{AB}=\nabla_A\phi\nabla_B\phi-g_{AB}\left[\tfrac12(\nabla\phi)^2+V\right],
\quad \Box_5\phi-V'(\phi)=0.
\]

Связность Levi–Civita; `∇^A T_AB=0` следует из scalar equation. Canonical scalar удовлетворяет `T_AB k^A k^B=(k·∇φ)²≥0` для null `k`.

### Базовый потенциал и пределы

Первый analytic/numerical milestone: `V=0`, `Λ₅=0` (massless scalar), потому что ровно этот предел поддержан конструкцией тест регулярных начальных данных и имеет Minkowski/Tangherlini benchmarks. Потенциальная версия допускается только как отдельно зафиксированный `action_id`, с гладким `V`, конечным `V,V'`, явно заданными minima и asymptotics. Для единственного AF vacuum `φ∞` требуется `Λ_eff=Λ₅+κ₅²V(φ∞)=0` и AF falloff. Ложновакуумный double well не становится базовым автоматически: нужен bounded potential, точный vacuum subtraction и пересчёт initial data/constraints.

Ограниченная первая задача не включает «вакуумный пузырь»: значение `φ` внутри некоторого радиуса не объявляется дочерней космологией. Чтобы моделировать phase transition, нужны минимум action, phase states, начальные perturbations и уравнение движения, а затем независимый causal criterion.

## 3. Что означает SO(4)-симметрия

Для сферической редукции orbit group `SO(4)` действует на трёхсферы `S³`; остаются 2D orbit space и areal radius. Это **5D spacetime с четырьмя пространственными измерениями**, не стандартная 3+1D космология и не «дочерняя 4D brane». Редукция допустима только если начальные данные, scalar и эволюционные boundary conditions инвариантны под SO(4). Она исключает nonspherical perturbations, rotation, anisotropic fragmentation и многие физические механизмы формирования дефекта.

Новая codimension-one `3+1D` brane потребовала бы embedding fields, localized matter/action, junction conditions и способа её образования. Симметричная поверхность `R=const` в сферическом bulk сама по себе не brane.

## 4. Четыре разных времени

- `t_coord` — координата выбранного численного gauge;
- `τ_P` — proper time конкретного parent observer worldline, если такой observer/geometry заданы;
- `τ_trans` — proper time переходной timelike worldline/hypersurface, если она построена;
- `t_D` — proper cosmological time только при выведенной/заданной homogeneous daughter metric.

В общем случае `dτ=αdt` для observer нормалей следует из lapse только в указанной slicing и worldline; это не превращает координатное время в независимые космические часы. Соотношение `dt_D/dT_P`, включая `K₀a_D`, не входит в физическую спецификацию без динамически определённого event map.

## 5. Модели, исключённые из первого solver

**ECSK**: требует независимой spin connection и заданного Dirac/spin matter action. вывод 5D ECSK контакта вывел контактный оператор, анализ микроскопического спинового состояния не определил коррелятор/знак макроскопического stress. До получения state-resolved `T_eff` ECSK нельзя добавить как `−C_s a^{-8}`.

**EGB/Lovelock**: другое действие с `α_GB` и иными bulk/boundary/junction equations; отдельная model branch, не `gravity toggle` в Einstein solver.

**Common bulk / brane**: одна 5D метрика может передавать causal bulk response между локализованными источниками, но это не две независимые 4D вселенные. В первой версии отсутствуют branes. Ни bimetric potential, ни RS redshift не подставляются как clock law.

**4D GR limit** не возникает из SO(4) редукции. Его можно получить лишь при дополнительной геометрии/компактификации или brane construction. Например, при product compactification с длиной `L₅` и низкоэнергетическим zero mode `κ₄²=κ₅²/L₅`; это отдельное условное dimensional-reduction утверждение, не feature этой 1+1 solver.

## 6. Граница научной интерпретации

Solver может показать regular initial data, collapse/dispersion, local trapped/marginal spheres, finite-time high-curvature growth или gauge/domain failure. Он не установит по одному run: геодезическую полноту, смену topology, существование второго asymptotic end, регулярное extension через breakdown, независимый daughter manifold или глобальную причинную изоляцию. Такие claims требуют полного future/global construction и отдельных theorem hypotheses.
