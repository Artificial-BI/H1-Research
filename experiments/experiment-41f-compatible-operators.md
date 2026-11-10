# Эксперимент 41F — проверка совместимости дискретных операторов

## Вердикт

Единый compatible/SBP discretization, который одновременно сохраняет Misner–Sharp mass flux и независимо удовлетворяет Hamiltonian/momentum constraints, в этой итерации **не получен**. #41E масс-баланс и stage-tangent ошибки воспроизведены, добавлена покомпонентная/по-face диагностика и режим эволюции массы по уравнению Einstein flux. Этот режим замыкает внешний ledger до round-off, но оставляет крупные независимые constraint defects. Поэтому он не сертифицирован как исправленный гравитирующий FV решатель; точное ledger conservation не считается подтверждением GR.

## Непрерывная система и проверка конвенций

Сохранено действие #40–41E, `c=1`, `kappa=8 pi G`, `V(phi)=lambda phi^2(phi-1)^2/4+delta(3 phi^2-2 phi^3)`, metric `ds^2=-alpha^2 dt^2+f^-1 dr^2+r^2dOmega^2`, `f=1-2m/r`. Конвенции: `Pi=phi_t/(alpha sqrt(f))`, `Phi=phi_r`; следовательно `phi_t=alpha sqrt(f)Pi`. Стресс-энергия даёт `rho=f(Pi^2+Phi^2)/2+V`, `p_r=f(Pi^2+Phi^2)/2-V`. Уравнения:

`m_r=(kappa/2)r^2 rho`,

`(ln alpha)_r=[m/r^2+(kappa/2)r p_r]/f`,

`Pi_t=r^-2 partial_r(r^2 alpha sqrt(f)Phi)-alpha/sqrt(f)V_phi`,

`m_t=F_m=(kappa/2)r^2 alpha f^(3/2)Pi Phi`.

На исходящей характеристике `Pi+Phi=0`, поэтому `F_m<0`: Misner–Sharp mass внутри фиксированной areal-radius поверхности убывает. В этих формулах `4 pi` уже включено в `kappa/2`; shell FV measure `DeltaV=(r_out^3-r_in^3)/3` учитывает геометрию `r^2 dr`, но не повторяет `4 pi`. Потенциал включён в объемную плотность и работу Klein–Gordon поля, а не в дополнительный граничный поток. При конечной границе lapse нормируется `alpha(R)=1`; это координатная нормировка времени именно на `R`.

Проверка кода не выявила неверного знака потока, пропущенного `f^(3/2)`, `alpha`, `r^2`, двойного `4 pi` или пропущенного временного RK4 коэффициента. RK4 использует `(1,2,2,1)/6` для тех же четырёх stage fluxes. Основная несовместимость — пространственная: нелинейная shell-average рекурсия массы и face-divergence скалярного поля не подчиняются общей дискретной цепной/product identity. В частности, midpoint approximation `S_i=Pi_i^2+(Phi_L^2+Phi_R^2)/2`, `V(phi_i)` и metric prefactors при интегрировании ограничения не дают точного discrete Leibniz cancellation с scalar FV RHS. Ошибка растёт и внутри области, поэтому одной коррекцией внешней boundary quadrature её убрать нельзя.

## Дифференцированный constraint audit на четырёх RK стадиях

Пусть constraint mode задаёт face mass `m_(i+1)=N_i/D_i`. На каждой стадии аналитически дифференцируется тот же отображатель:

`dot m_(i+1)=(dot N_i-m_(i+1) dot D_i)/D_i`,

`dot S_i=2Pi_i dot Pi_i+Phi_L dot Phi_L+Phi_R dot Phi_R`,

`dot V_i=V_phi(phi_i) dot phi_i`.

Этот `dot m_H` на каждом face сравнивается с независимо вычисленным `F_m`; код сохраняет максимум по каждому face и временной интеграл независимого `gradient_t(m)-F_m`. Это выявляет нарушение дискретной identity, а не повторно использует outer ledger.

Для pulse `A=.35`, `R0=.5`, `kappa=.1`, `R=4`, duration `mu t=4` stagewise relative L2 mismatch в constraint mode равен `.200296` при `dr=.1,dt=.02`, `.200304` при `dr=.1,dt=.01`, `.056318` при `dr=.05,dt=.01`, `.014288` при `dr=.025,dt=.005`. Уменьшение `dt` вдвое почти не влияет; уменьшение `dr` снижает defect примерно со вторым порядком. Пространственные ошибки показывают, что схема приближается к continuum identity, но дефект не исчезает на текущих разрешениях.

Для `dr=.1` независимый time-integral локального residual до `r=3.9` около `-1.78e-5`, на `R=4` он `+3.5243e-5`; endpoint RK4 ledger mismatch `+3.4621e-5`. Stage quadrature и time-gradient дают согласующиеся integrated остатки. При `R=5`/`6` boundary power при том же duration практически нулевой, но ledger remainder остаётся `1.0043%`/`1.0243% M0`. Следовательно основной дефект — bulk operator compatibility, усиленный interface при реальном прохождении импульса.

## Предустановленный критерий и результаты

Для новых сравнений задан screening gate `1e-3` отдельно для Hamiltonian L2, momentum Linf, stage tangent L2, local mass-flux L2 и outer ledger/M0; общий numerical gate pass требует выполнения всех пяти. Исторические #41E значения были известны до #41F, поэтому эти тесты exploratory, а gate не является слепой preregistration. Это критерий отбраковки схемы, не «подтверждение физики».

### Исходящая граница: constraint-reintegrated против mass-flux evolution

Граничный тест `A=.35`, `R0=.5`, `R=4`, duration 4:

| dr | dt | constraint H L2 | constraint momentum Linf | constraint local flux L2 | constraint ledger/M0 | mass-flux H L2 | mass-flux momentum Linf | mass-flux ledger/M0 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .2 | .04 | .10994 | .63802 | .18694 | .1370% | .32335 | .50271 | 1.66e-16 |
| .1 | .02 | .03448 | .30124 | .05502 | 1.9621% | .09974 | .26600 | 3.44e-16 |
| .05 | .01 | .009536 | .10075 | .01513 | .6496% | .03191 | .09727 | 9.77e-17 |
| .025 | .005 | .002388 | .02741 | .004085 | .4023% | .01446 | .02716 | 1.56e-16 |

В constraint mode `dr=.1`: `M0=.00176454805`, `M1=.00097651779`, `DeltaM=-.00078803026`, `P_out,max=.00296230096`, stage-integrated flux `-.00082265170`, остаток `+3.46214e-5=1.9621% M0`. В mass-flux mode то же уравнение `dot m=F_m` эволюционируется прямо на faces четырьмя RK stages. Тогда интегральный бюджет замыкается до round-off, но при `dr=.025` H L2 `.01446` и momentum Linf `.02716`, оба существенно выше `1e-3`. Это закрытый консервативный ledger **с независимым constraint failure**, не совместимый успешный решатель.

Для pulse при `dr=.1`, абсолютные ошибки в solver units сохранены отдельно от нормированных: constraint mode Hamiltonian abs Linf/RMS `3.954e-5 / 2.668e-6`, momentum abs Linf/RMS `2.664e-3 / 3.572e-5`, local flux abs Linf/RMS `4.063e-4 / 4.121e-5`, integrated ledger abs `3.462e-5`. Mass-flux mode: Hamiltonian `7.425e-5 / 5.368e-6`, momentum `2.352e-3 / 3.278e-5`, local flux `3.331e-5 / 3.921e-6`, ledger `6.064e-19`. RMS Hamiltonian/momentum используют shell-volume weighting; local flux RMS использует face samples. Низкое абсолютное RMS не заменяет относительный gate для локализованного поля.

### Фиксированное поле на границе

Для сравнения реализована фиксированная координатная поверхность с Dirichlet `phi(R)=0`: нечётный ghost даёт `Phi_R=-2 phi_last/dr`, а `Pi_R=0` зануляет физический скалярный flux `F_m(R)`. Слабый профиль `A=.1`, `R0=.5`, `R=4`, `dr=.1`, `dt=.01`, duration 4 имеет точно нулевую surface power, но constraint-mode ledger drift `1.3584%`, H L2 `.03487`, momentum Linf `.3042`. Отдельный outgoing run даёт ненулевую мощность `2.47696e-4` и ledger drift `1.8488%`. Условие на границе меняет физический поток, но не убирает bulk дисбаланс.

### Слабое поле и Minkowski reference

Для `A=.1`, `R0=.5`, `R=5`, duration 2 относительная L2-разность конечного `phi` coupled solution от #41C fixed-Minkowski scalar reference равна `.004485/.002520/.002452` при `dr=.2/.1/.05`. При самой тонкой из этих сеток coupled H L2 `.009065`, momentum Linf `.09922`, local flux L2 `.02184`, stage mismatch `.06204`, ledger `.2478%`. Сходимость field curve сама по себе не проходит критерии допустимости источника.

Точный вакуум `phi=Pi=0` сохраняет `m=0, alpha=1, f=1` и получает numerical gate pass именно для этого тривиального fixed point. Новую сильнополевую эволюцию не продолжали: weak/intermediate проверки не прошли gate. Исторический короткий strong test #41E остаётся кратким полевым интегралом без горизонта, не новой проверкой #41F.

## Исходные #41D/#41E и состояние реализации

Исторический #41D constraint scheme: outgoing pulse budget residual `1.962%` при `dr=.1`; ранний #41B legacy solver: outer mass loss `-2.1711%` при наложенном нулевом потоке. #41C Minkowski FV scalar reference достигал ledger residual `1.99e-6/6.66e-8/2.06e-9` на `dr=.2/.1/.05`, но это не проверка связанной гравитации.

Добавлен опциональный `evolution="mass_flux"` и независимый `1e-3` numerical gate в `coupled_fv.py`; UI в существующей Experiment 41 странице сравнивает outgoing constraint/mass-flux варианты и Dirichlet control. Legacy #40/41D массивы и селекторы не менялись. Единого compatible/SBP backend в этой ревизии нет — он остаётся следующим численным этапом. Нельзя утверждать, что coupled solver отремонтирован.

Polar-areal chart всё ещё останавливается при `f<=1e-5`; нет horizon-penetrating continuation, apparent/event horizon audit или global causal diagram. H1 clock `dt_d/dT_p=K0 a_d` остаётся NOT_COMPARABLE: здесь используется `mu t`, не вычислена дочерняя FLRW `a_d`, не задан физически связанный parent observer. Никакой child-universe или parent-time вывод не следует.

Полный регрессионный набор: `98` прежних тестов плюс 2 новых benchmark controls, итого `100 passed` (`71.55 s`) командой `\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider`; Qt offscreen test входит в полный набор. `\.venv\Scripts\python.exe -m compileall -q src tests` успешен. Локальное окружение — проектный `.venv`, зависимости не менялись.
