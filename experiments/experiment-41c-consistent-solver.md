# Experiment 41C — совместимость дискретизации и бюджет массы

## Вывод

**Связанный Einstein–scalar solver пока не восстановлен.** Непрерывные уравнения и массовый поток согласуются аналитически, однако две дискретные формы #40/#41B дают разные неустранимые в текущем коде ошибки:

1. При повторном решении радиального гамильтонова ограничения конечная Misner–Sharp масса меняется на `−2.1711%` при нулевом потоке на внешнем крае.
2. При эволюции `m` по Einstein mass-flux equation ledger закрывается, но гамильтоново ограничение остаётся с остатком `~1.45e−3` при измельчении сетки.

В этой итерации добавлен **консервативный flat-space сферический Klein–Gordon эталон** с конечными объёмами. Он подтверждает корректный scalar energy ledger и измеряет ненулевой поток через границу, но гравитация в нём выключена. Это демонстрация необходимой дискретной структуры, **не ремонт связанного solver** и не свидетельство о дочерней области.

## Непрерывная теория и проверка соглашений

Действие и потенциал #40:

`S=∫d⁴x √−g [R/(2κ) − (∇φ)²/2 − V(φ)]`, `κ=8πG`,
`V=λφ²(φ−1)²/4 + δ(3φ²−2φ³)`.

Для сигнатуры `−+++` канонический тензор энергии-импульса

`T_ab=∇_aφ∇_bφ − g_ab[(∇φ)²/2+V]`.

В полярно-ареальной gauge без shift:

`ds²=−α²dt²+A²dr²+r²dΩ²`, `f=A⁻²=1−2m/r`,
`Π=(A/α)φ_t=A n(φ)`, `Φ=φ_r`.

Тогда

`ρ=T_ab n^a n^b=f(Π²+Φ²)/2+V`,
`p_r=T^r_r=f(Π²+Φ²)/2−V`.

Поскольку `m=G M_MS/c²` — длина Misner–Sharp, сферические уравнения Эйнштейна и Klein–Gordon при `c=1` имеют вид

```
m_r = (κ/2) r² ρ
(ln α)_r = [m/r² + (κ/2) r p_r]/f
φ_t = α Π/A
Π_t = r⁻² ∂r(r² α Φ/A) − α A V_φ
```

При соглашении `K_ij=−(1/2)L_n γ_ij`:

`K^r_r=−(κ/2)rΠΦ/A`.

Из смешанного уравнения Эйнштейна (или Misner–Sharp identity)

`m_t=(κ/2)r² T_t^r=(κ/2)r² α f^(3/2)ΠΦ`.

Знак важен: `T_t^r=f φ_tφ_r=α f^(3/2)ΠΦ`. Поток наружу по энергии обозначим `F_E=−T_t^r`; тогда при `F_E>0` энергия внутри сферы убывает: `m_t=−(κ/2)r² F_E`.

Непрерывная constraint/flux identity получена напрямую. Обозначим `c=κ/2`, `S=Π²+Φ²`. Подстановка KG в производную гамильтонова constraint даёт

```
(m_rt − m_tr)/c
 = r² α sqrt(f) Π Φ [ f (ln α)_r − f_r/2 − c r f S ].
```

Из `f_r=−2m_r/r+2m/r²`, `m_r=c r²ρ` и lapse equation следует

`f(ln α)_r−f_r/2 = c r(ρ+p_r)=c r fS`,

следовательно `m_rt=m_tr`. В коде #40 параметры `Π`, знаки, потенциал, непрерывные формулы mass/lapse/KG и flux **совпадают** с этими уравнениями; явной неверной метрики или множителя в массовом потоке не найдено.

Регулярный центр требует чётного `φ,Π`, `Φ(0)=0`, `m(0)=0`; площадь сферического потока у центра равна нулю. Внешнее условие #40 `φ=Π=0` не является outgoing boundary condition и может отражать достигшую края волну. Однако в исходном интервале `μt≤2` световой сигнал от `r=2` проходит только до `r≈4.0307` при `r_max=16`; расширение края до 24/32 не меняло массовый drift существенно. Поэтому отражение не объясняет этот конкретный baseline.

## Дискретный root discrepancy

Исходный collocation solver получает `Φ` через `np.gradient`, дивергенцию `r²αΦ/A` — отдельной центральной разностью, `m` — трапецеидальной recurrence, а `m_t` после прогона — конечной временной разностью. Эти операторы не образуют дискретную пару summation-by-parts и не удовлетворяют дискретному chain/product rule, который в непрерывной теории сокращает lapse, давление и kinetic terms. Поэтому constraints, scalar update и `m_t` не сохраняют один и тот же дискретный бюджет.

В исправленной диагностике #40 mass change за шаг сравнивается с RK4-квадратурой потока на всех четырёх стадиях, а не с одним конечным stage. Intermediate max residual изменился `1.76018e−3 → 1.66388e−3`, strong `0.151851 → 0.150390`; это устраняет временную несогласованность метрики residual, но не массовый drift. Drift также практически неизменен при `dt/2`, `dt/4` и доменах `r_max=16,24,32`. Для constraint-reintegrated базы `M(0)=.0970447857396`, `M(2)=.0949378821303 μ⁻¹`, `ΔM=−.00210690360926 μ⁻¹`; измеренный boundary flux и его integral равны нулю.

Массово-эволюционная форма #41B является независимым forensic probe: она решает начальный `m_r` constraint, затем интегрирует `m_t` на RK4 стадиях и мониторит `m_r` residual. В ней внешняя масса сохраняется, но H residual не убывает к нулю: `.0018463` (`dr=.2`), `.0015543` (`.1`), `.0014751` (`.05`), `.0014559` (`.025`). Сильный вариант при `A=1.8,R=1,κ=.5` имеет `max(2m/r)=.43244` и H residual `.03094`. Нулевой mass drift в этой ветви — следствие интегрирования именно mass equation; он не считается независимым доказательством корректности. Результат указывает на bulk discrete compatibility failure, но не изолирует один конкретный операторный дефект.

## Сферический finite-volume reference: гравитация выключена

Чтобы отдельно проверить scalar flux/energy boundary, реализован reduced model для той же скалярной переменной и потенциала на фиксированном Minkowski фоне. `φ_i,Π_i` расположены в радиальных cell centers. Для faces `r_j=jh`:

`V_i=(r_{i+1}³−r_i³)/3`,
`Φ_j=(φ_j−φ_{j−1})/h` на внутренних faces,
`Φ_0=0` (регулярный нулевой area flux),
`Φ_N=−Π_{N−1}` (outgoing Sommerfeld characteristic).

Полудискретная система:

```
φ̇_i = Π_i
Π̇_i = [r_{i+1}² Φ_{i+1} − r_i² Φ_i]/V_i − V_φ(φ_i)
```

Дискретная энергия без множителя `4π`:

`E_h=Σ_i V_i[Π_i²/2+V(φ_i)] + (h/2)Σ_{j=1}^{N−1}r_j²Φ_j²`.

Суммирование по частям на конечных объёмах сокращает все внутренние face fluxes и даёт точное полудискретное тождество

`dE_h/dt=R² Π_N Φ_N=−R²Π_N²`.

Положительный outgoing power — `P_out=R²Π_N²`; budget — `E(t)−E(0)−∫R²ΠΦ dt=0`. Временная интеграция использует RK4 stage quadrature; её ledger remainder должен сходиться к нулю при CFL refinement. Это обычная scalar field на фиксированном Minkowski фоне; Einstein constraints, Misner–Sharp mass, lapse backreaction и horizons здесь **не вычисляются**. При наличии массы поля Sommerfeld — приближение к физическому бесконечному радиационному условию; flux ledger всё равно точно учитывает энергию, введённую выбранным внешним условием.

Результаты для фиксированного профиля `A=.35`, `R0=1.5`, `r_max=8`, `t_end=8`, `dt/dr=.2`:

| `dr` | `dt` | `E(0)` | `E(end)` | интеграл signed boundary power | ledger remainder | max relative ledger residual | max outward power |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .2 | .04 | .127321365 | .098552051 | −.028769077 | −2.3733e−7 | 1.991e−6 | .040208 |
| .1 | .02 | .127831607 | .097775469 | −.030056132 | −6.2262e−9 | 6.656e−8 | .038597 |
| .05 | .01 | .127960561 | .097612759 | −.030347802 | −5.7965e−11 | 2.063e−9 | .037994 |

At fixed `dr=.1`, timestep checks give max relative ledger residual `6.656e−8` for `dt=.02`, `2.005e−9` for `.01`, and `5.845e−11` for `.005`; energy loss and integrated boundary flux agree within those errors. Halving both `dr` and `dt` reduces the max relative ledger remainder by factors about 30 and 32. The regression gate `<1e−3` was added after the exploratory measurements and is therefore a future code-regression guard, not a blind/pre-registered scientific acceptance test. The data cannot be reused as evidence for the coupled solver.

## Before/after and scientific status

| Check | #41B coupled constraint form | #41B mass evolution | #41C reference result |
|---|---|---|---|
| finite-radius mass/energy ledger | fails by −2.1711% at zero boundary flux | closes by construction; Hamiltonian residual fails | scalar finite-volume ledger closes and matches nonzero outgoing power |
| Hamiltonian constraint | constraint is re-solved; discretized flux identity fails | max residual `.001456` at `dr=.025`, strong `.03094` | gravity absent; constraint N/A |
| physical strong field | invalid, max `2m/r=.4324` | invalid, same compactness and H error | no gravity/horizon claim |
| child causal region | NOT COMPUTED | NOT COMPUTED | NOT APPLICABLE |

Thus this work diagnoses a discrete compatibility class and demonstrates a correct finite-volume scalar budget, but it **does not repair simultaneous constraint and mass-budget errors in the coupled Einstein–scalar model**. A complete repair requires a compatible SBP/finite-volume formulation for radial constraints, lapse, scalar stress and mass flux with a discrete Bianchi/constraint-propagation proof; both initial and evolution constraints must converge together. The outgoing scalar boundary must be imposed by incoming characteristic data, not fixed values for both `φ` and `Π`.

The local code still uses `μt`; it does not define a homogeneous daughter `a_d` or parent observer. H1 `dt_d/dT_p=K0a_d` remains **NOT_COMPARABLE**. No `T_p` seconds, daughter formation, horizon, or H1 prediction is inferred.

## Verification

The added tests check the three-grid outgoing energy ledger, nonzero boundary power, refinement trend, fixed-profile mass/H constraint mismatch, and Qt offscreen diagnostics. Full regression command from the project root:

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
uv run --no-sync python -m pytest -q
uv run --no-sync python -m compileall -q src tests
```

Passes establish code/regression checks only; they do not validate a strong-field Einstein–scalar solution.
