# Независимый double-null benchmark Эйнштейна–безмассового скаляра

**Область:** узкий классический numerical control для проверки метода на massless canonical scalar. Не эксперимент о рождении дочерней вселенной, не подтверждение H1 и не воспроизведение критического exponent.

## Итог

Добавлен отдельный double-null conformal-coordinate reference `double_null_reference.py` (файл не включён в публичный экспорт). Он не вызывает RHS, radial constraints, mass flux или FV operator из #40/#41F. Начальные Einstein constraints интегрируются отдельным radial RK4, а последующая эволюция использует mixed double-null equations на `(t,x)`-сетке с методом линий и RK4.

**Результат — partial reference, не валидированный сильнополевой решатель.** Minkowski vacuum и аналитический Schwarzschild exterior Misner–Sharp mass проходят ограниченные проверки. При слабом поле scalar profile сходится к независимой точной радиальной волне. Для гравитирующей семьи null constraints и mass-balance residual убывают при refinement, но на `dr=.05` они остаются большими; сравнение с исходным polar-areal решателем тоже пока расходится. Поэтому STOP до амплитудного порога, trapped-collapse поиска и Choptuik scaling.

## 1. Соглашения и уравнения

Действие — ровно исходное действие #40 в exact massless limit:

\[
S=\int d^4x\sqrt{-g}\left[\frac{R}{2\kappa}
-\frac12\nabla_a\phi\nabla^a\phi\right],
\qquad G_{ab}=\kappa T_{ab},\quad \Box\phi=0.
\]

Потенциал исходной модели `V=lambda*phi²*(phi−1)²/4+tilt*(3phi²−2phi³)` зануляется точно при `lambda=tilt=0`. Валидатор #40 теперь разрешает только эту пару как дополнительный точный контрольный предел; стандартные `lambda=1`, `tilt=.02` и их ограничения прежние. В #41F конечный объём уже принимал нулевые коэффициенты. `kappa=8 pi G=.1` в dimensionless Planck-normalized units; `mu` не откалибрована в SI.

Метрика определена как

\[
ds^2=-e^{2\sigma(u,v)}du\,dv+r^2(u,v)d\Omega^2,
\quad g_{uv}=g_{vu}=-e^{2\sigma}/2,
\quad u=t-x,\quad v=t+x.
\]

Тогда `ds²=e^(2σ)(−dt²+dx²)+r²dΩ²`; координатная скорость радиальных нулевых лучей `dx/dt=±1` точно. Не следует путать `t` с proper time: для координатного наблюдателя `x=const`, `dτ=e^σ dt`; его worldline обычно ускорена.

Из `G_uu=κ phi_u²`, `G_vv=κ phi_v²`, `R_uv=κ phi_u phi_v` и `Box phi=0` получены уравнения:

\[
r_{uu}-2\sigma_u r_u=-\frac\kappa2 r\phi_u^2,
\qquad
r_{vv}-2\sigma_v r_v=-\frac\kappa2 r\phi_v^2,
\]

\[
r_{uv}=-\frac{r_ur_v}{r}-\frac{e^{2\sigma}}{4r},
\]

\[
\sigma_{uv}=\frac{r_ur_v}{r^2}+\frac{e^{2\sigma}}{4r^2}
-\frac\kappa2\phi_u\phi_v,
\qquad
\phi_{uv}=-\frac{r_u\phi_v+r_v\phi_u}{r}.
\]

Факторы `1/4` фиксируются конвенцией `g_uv=−e^(2σ)/2`. В координатах `(t,x)` RHS evolution записан прямо в коде:

\[
r_{tt}=r_{xx}+\frac{r_x^2-r_t^2-e^{2\sigma}}r,
\]

\[
\sigma_{tt}=\sigma_{xx}+\frac{r_t^2-r_x^2+e^{2\sigma}}{r^2}
-\frac\kappa2(\phi_t^2-\phi_x^2),
\]

\[
\phi_{tt}=\phi_{xx}+\frac{2(r_x\phi_x-r_t\phi_t)}r.
\]

Для плоского вакуума `r=x`, `sigma=0`, `phi=0` ускорения тождественно нулевые. Центровые пределы не вычисляются делением на `r=0`; используются Taylor-пределы и parity. Результаты null constraints — независимые finite-difference оценки

\[
C_{uu}=r_{uu}-2\sigma_u r_u+\frac\kappa2r\phi_u^2,
\qquad
C_{vv}=r_{vv}-2\sigma_v r_v+\frac\kappa2r\phi_v^2.
\]

### Misner–Sharp масса, поток и расширения

\[
1-\frac{2m}{r}=g^{ab}r_ar_b=-4e^{-2\sigma}r_ur_v,
\qquad
m=\frac r2\left(1+\frac{4r_ur_v}{e^{2\sigma}}\right)
=\frac r2\left(1+\frac{r_t^2-r_x^2}{e^{2\sigma}}\right).
\]

С использованием двух null constraints и `r_uv` получаются уравнения массового потока

\[
m_u=-\kappa r^2 e^{-2\sigma}r_v\phi_u^2,
\qquad
m_v=-\kappa r^2 e^{-2\sigma}r_u\phi_v^2.
\]

Код отдельно вычисляет геометрическую массу из `r,r_t,r_x,sigma`, затем временно/пространственно дифференцирует её и сравнивает `m_u,m_v` с этими выражениями. Это лучше, чем просто интегрировать mass-flux equation в переменную `m`, но всё ещё использует производные того же численного решения и не заменяет сравнение с независимой reference geometry.

Будущие радиальные null expansions при координатных null normals `∂v`,`∂u`:

\[
\theta_{out}=2r_v/r=(r_t+r_x)/r,
\qquad \theta_{in}=2r_u/r=(r_t-r_x)/r.
\]

Их масштаб зависит от нормировки null vectors; знаки и произведение имеют геометрический смысл. Apparent/marginally trapped sphere — локальная характеристика среза. Event horizon требует знания глобального будущего, которого конечная сетка не предоставляет.

## 2. Начальные данные, gauge и края

Используется time-symmetric smooth family

\[
\phi(0,r)=A e^{-\frac12(r/R_0)^4},\quad \Pi=0,\quad K_{ij}=0,
\]

с регулярным центром, быстро убывающим хвостом и одним конечным внешним доменом. Это заранее внесённый скалярный профиль, а не спонтанное образование из пустого родителя.

На начальном срезе независимо интегрируется polar-areal Hamiltonian constraint

\[
m_r=\frac\kappa4r^2f\phi_r^2,\quad f=1-2m/r,
\qquad
(\ln\alpha)_r=\frac{m/r^2+(\kappa/4)r f\phi_r^2}{f}.
\]

Lapse фиксируется `alpha(r_max)=1`. В tortoise-like coordinate задаётся

\[
dx/dr=1/(\alpha\sqrt f),\qquad r_x=\alpha\sqrt f,
\qquad \sigma(0,x)=\ln\alpha(r(x)).
\]

Поэтому начальная пространственная метрика совпадает с constraint-solved polar-areal Cauchy-срезом и `t` синхронизировано с исходным временем на начальном срезе. Начальные скорости `r_t=sigma_t=phi_t=0`.

В центре imposed regularity: `r=0`, odd parity `r,r_t`, even parity `sigma,phi,sigma_t,phi_t`, `e^sigma(0)=r_x(0)`. Для удержания условий центральные значения `sigma,sigma_t` проецируются на Taylor-отношения на каждом RK stage. Поэтому нулевой center regularity residual частично enforced численно, а не независимое свидетельство точности. Используются centered second-order finite differences, их regular-center Taylor limits, method-of-lines RK4 и `dt≈CFL*dx`, где `dx` — шаг uniform optical x-coordinate. Сетка начинается от равномерного ареального `dr`, затем начальные данные интерполируются на uniform x.

Внешний `x=x_max` — конечная фиксированная coordinate boundary, значения `r,sigma,phi` заморожены; это не асимптотическая бесконечность и не доказанная поглощающая граница. В отчётных коротких run-ах волна не должна успеть вернуться к анализируемой центральной области. Продление за outer boundary, за будущую бесконечность или до глобальной причинной диаграммы не реализовано. Double-null conformal coordinates могут быть регулярны на обычном Schwarzschild horizon, но конкретная сетка ещё не протестирована на horizon crossing или `r=0` после сильного сжатия.

Две intersecting characteristic legs **не** задавались как произвольные pulses: начальные данные заданы на spacelike Cauchy surface и решено её spherical Hamiltonian constraint. Null constraints затем проверяются по эволюции. Это избегает неразрешённой неоднозначности предписания независимых данных на двух null legs, но ограничивает задачу одним конкретным time-symmetric профилем.

## 3. Предзаданный набор и измерения

### Контроли

1. `A=0`: Minkowski fixed point.
2. Положительная масса vacuum Schwarzschild: только аналитическая invariant-проверка static exterior `r>2M`, где `Omega²=f`, `r_t=0`, `r_x=f`, `sigma=.5 ln f`. В регулярном центре vacuum по Birkhoff допустима лишь `M=0`; solver не включает сингулярное начало Schwarzschild.
3. `A=.01`, `kappa=1e-8`: weak-gravity profile против exact flat radial solution с одинаковыми `A,R0`, `r<8`.
4. `A=.5`, `kappa=.1`: контроль coupled Einstein evolution при `r_max=16`, до `mu*t=1`.

Основные сетки `dr=.2,.1,.05`, `CFL=.2`; дополнительно `dr=.1` с `CFL=.1` и изменение `r_max=12,16,24`. Для исходящего pulse через близкую границу в этом iteration double-null boundary не является корректным поглощающим условием, поэтому отдельное outgoing budget через boundary к этому solver не приписывается.

Для FV был унаследован заранее установленный #41F screening gate `1e-3`; в этой новой double-null схеме он использовался как ориентир, не как preregistered/слепой тест. Численное acceptance criteria для double-null диагностики не было надёжно пререгистрировано перед экспериментальными прогонами. Ниже поэтому результаты — exploratory benchmark с сохранёнными ошибками, а не formal pass.

### Вакуум и аналитическая геометрия

- `A=0`, `r_max=16`, `dr=.2`, duration `.2`: `M_init=M_final≈−1.81e−13` (floating cancellation; физический ноль), max absolute null-constraint L2 `2.41e−14 / 2.39e−14`, max compactness `4.00e−15`, `min(1−2m/r)=1−4e−15`. Относительная норма при нулевом источнике не определена и не интерпретируется как accuracy.
- Schwarzschild exterior helper для `M=1.4`, `r=3,4,8,20` возвращает `m=1.4` с floating-point точностью. Это invariant formula check, не evolved curved solution.

### Слабое поле против точного плоского решения

Для `A=.01`, `R0=.5`, `kappa=1e−8`, `r_max=16`, `CFL=.2` и `t=.5` максимальная абсолютная ошибка `phi` на `r<8`:

| dr | max `|phi_DN−phi_flat|` |
|---:|---:|
| .20 | `2.86134e−3` |
| .10 | `1.10074e−3` |
| .05 | `3.10730e−4` |

Effective order на последовательных парах примерно `1.38`, `1.82`; поле монотонно сходится к независимой analytic flat-wave solution. Гравитационный источник здесь преднамеренно почти нулевой, поэтому relative Einstein constraints/flux в этом лимите ill-conditioned; абсолютные null-constraint L2 на трёх сетках `3.46e−14`, `1.44e−13`, `1.54e−12` для `uu` и `3.97e−14`, `1.43e−13`, `1.62e−12` для `vv`. Они растут из-за конечных разностей и геометрических cancellation, хотя остаются малыми в dimensionless units. Этот тест проверяет scalar principal part, не gravity.

### Coupled-срезы и независимые null-диагностики

Для `A=.5`, `R0=.5`, `kappa=.1`, `r_max=16`, `mu*t=1`, `CFL=.2`:

| dr | Initial Misner–Sharp M | Final M at outer coordinate edge | max `Cuu` relative L2 | max `Cvv` relative L2 | max flux `mu/mv` relative L2 | max mass ledger / M0 | min `1−2m/r` | max `2m/r` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .20 | `.00351104` | `.00351107` | `.2512` | `.5818` | `1.394 / 1.581` | `33.34%` | `.98616` | `.01384` |
| .10 | `.00352084` | `.00352085` | `.1650` | `.3882` | `.8866 / 1.074` | `18.85%` | `.98627` | `.01373` |
| .05 | `.00352092` | `.00352092` | `.07282` | `.1047` | `.6635 / .8729` | `11.18%` | `.98379` | `.01621` |

`Cuu/Cvv` уменьшаются, но остаются далеко выше `1e−3`; массовые потоки также не validated. Mass-ledger здесь — максимум по пространственным coordinate rays между `m(t,x)-m(0,x)` и временным интегралом `m_t=m_u+m_v`; это не только outer boundary loss. Его absolute maximum `1.171e−3`, `6.638e−4`, `3.936e−4` в geometry mass units. Outer edge искусственно фиксирован, поэтому почти неизменная `M_outer` сама по себе не является физической conservation validation.

Эффективные порядки max relative constraints неустойчивы: для `Cuu` около `.61` затем `1.18`; для `Cvv` около `.58` затем `1.89`; для массового ledger примерно `.82` и `.76`. На выбранных сетках это не demonstrable asymptotic convergence до точности, требуемой для физической интерпретации. Изменение `CFL=.2→.1` при `dr=.1` дало `Cuu/Cvv=.16496/.38816 → .16566/.38809`, ledger `18.8539%→18.8511%`: temporal refinement здесь почти ничего не меняет.

Минимальный `f` остаётся положительным, максимум компактности `~.016`; trapped surface не образовалась. Это лишь отсутствие локального trapped sphere на рассмотренной слабой траектории, не доказательство отсутствия event horizon в полном будущем.

## 4. Сверка с polar-areal Cauchy data

На начальном срезе профили совпадают в continuum: `A=.5`, `R0=.5`, `Pi=0`, один `kappa`, конечная asymptotically flat boundary. Double-null initial constraint интегрируется отдельным RK4, старый #40 — собственной trapezoidal constraint-рекурсией. Максимальная относительная разность Misner–Sharp `m(r)` после интерполяции на общие areal радиусы:

| dr | max initial `|m_DN−m_PA| / M_PA(R)` |
|---:|---:|
| .20 | `33.55%` |
| .10 | `8.80%` |
| .05 | `2.07%` |

Это убывает примерно вторым порядком, но fine-grid discrepancy остаётся заметной. После эволюции coordinate arrays не сравнивались как физические предсказания. Центральный observer proper time вычисляется как `tau_c=integral exp(sigma(t,0))dt`; для polar-areal — `integral alpha(t,0)dt`. `phi` при равном накопленном центральном proper time различалось на `.391`, `.117`, `.0306` на `dr=.2,.1,.05` за интервал порядка единицы. Это конвергирует, но недостаточно мало и #40 сам не валидирован. Нормировки координатных часов привязаны к конечным boundary observers, траектории которых после начального среза не идентичны; поэтому сравнение — diagnostic, не полная gauge-invariant inter-code validation.

## 5. Литература, пределы и решение

Офлайн локальная `references.md` не содержит работ Choptuik/Garfinkle; DOI ниже — предоставленные владельцем задачи библиографические данные, онлайн не проверялись:

- M. W. Choptuik, “Universality and scaling in gravitational collapse of a massless scalar field,” *Phys. Rev. Lett.* 70, 9 (1993), DOI `10.1103/PhysRevLett.70.9`.
- D. Garfinkle, “Choptuik scaling in null coordinates,” *Phys. Rev. D* 51, 5558 (1995), DOI `10.1103/PhysRevD.51.5558`.

`gamma≈0.37` — ориентир near-critical high-resolution massless problem. Здесь нет амплитудной bisection, нескольких декад scaling, AMR, apparent-horizon formation или Garfinkle cross-code benchmark; exponent **не вычислялся**.

| Gate | Итог |
|---|---|
| Exact flat vacuum and static-exterior mass identity | **PASS, restricted analytic controls only** |
| Weak scalar principal propagation | **PARTIAL; monotonic 3-grid convergence to flat solution** |
| Independent null Einstein constraints | **FAIL as accuracy gate; decreasing but large residuals** |
| Double-null mass-flux identities and mass ledger | **FAIL as accuracy gate; large residuals** |
| Equivalent initial data vs #40 | **PARTIAL; mass discrepancy decreases but remains 2.07% at dr=.05** |
| Dispersal vs black-hole threshold / strong field | **NOT TESTED** |
| Event horizon / global causal future / baby universe | **NOT COMPUTED** |
| H1 clock `dt_d/dT_p=K0 a_d` | **NOT_COMPARABLE**; no daughter FLRW `a_d` or linked parent observer |

**Решение: STOP** для amplitude scan, горизонта и физической интерпретации коллапса этой реализацией. Следующий инженерный шаг, если будет отдельно одобрен: устранить несходимость независимых null constraints/mass-flux и boundary ledger на том же фиксированном профиле; до успешной проверки независимую characteristic formulation не расширять к критичности. Даже валидированный стандартный collapse solver не подтвердит образование дочерней геометрии H1.

## 6. Тесты

Новый double-null набор: **5 passed**; вместе с предыдущими четырьмя massless collapse tests добавлено 9 тестов. Полный регрессионный запуск `\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider`: **109 passed in 105.01 s**. Он включает прежние Qt offscreen проверки. `\.venv\Scripts\python.exe -m compileall -q src tests` завершился успешно. Тестовые PASS не означают physical solver validation; ручного desktop GUI теста не было, UI не менялся.
