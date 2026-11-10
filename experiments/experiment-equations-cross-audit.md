# Первопринципный аудит уравнений Эйнштейна–скаляра в двух калибровках

**Тип работы:** математическая сверка исходников, ручное выведение от записанного действия и один read-only initial-slice probe. Эволюционные расчёты, тесты, solver/UI изменения и новые допуски в этой работе не выполнялись. Числа прежних экспериментов используются как ранее полученные наблюдения, а не повторно подтверждённые результаты.

## Краткий вывод

1. Полярно-ареальная метрика в реальном коде — `ds²=−alpha² dt²+f⁻¹ dr²+r²dΩ²`, а **не** `−alpha² f dt²+f⁻¹dr²`. В коде `a=f⁻¹/²`, `Pi=a phi_t/alpha=phi_t/(alpha sqrt(f))`, `Phi=phi_r`. При этом соглашении KG first-order система, энергия, lapse, масса и mass flux согласуются с действием. Единственная обнаруженная ошибка в пояснительной строке — формула `K^r_r` без `sqrt(f)`; реализация и соседнее momentum constraint фактор содержат.
2. Для double-null конвенции `ds²=−Omega² du dv+r²dOmega²`, `g_uv=−Omega²/2`, уравнения `r_uv`, `sigma_uv`, `phi_uv`, null constraints, Misner–Sharp масса и `m_u/m_v` в отчёт о double-null проверке алгебраически совпадают с независимым выведением для `V=0`. Кодовая пара `u=t−x`, `v=t+x` корректна **в эволюционном conformal gauge**, который сам код задаёт для всего решения.
3. Double-null `x(r)` из начального polar-сечения не является общим преобразованием координат во всё будущее, если применить его к меняющимся `alpha(t,r),f(t,r)`. В отчёт о double-null проверке он используется только для преобразования начальных данных; дальше развивается новая double-null метрика. Это законный gauge choice, но будущие точки двух расчётов не совпадают по общему `t,x` автоматически.
4. Прежнее расхождение начальных масс **не объясняется разным `kappa`, потенциалом, начальной кинетикой или lapse-нормировкой**. Разрешённый read-only probe обнаружил разные конечносеточные initial constraints: полярный solver использует `np.gradient(phi)` и трапецеидальную неявную рекурсию, double-null — аналитический `phi_r`, RK4 radial интеграцию, затем интерполяцию `r(x),sigma(x)` и производную `r_x`. На `dr=.05` разность total mass примерно `1.67%` (max profile discrepancy `2.07%` относительно polar outer mass). Она убывает примерно как `dr²`; это измеренное численное несовпадение начальных данных, а не физическое различие gauge.
5. Double-null mass-ledger дефект `11.18%` при `A=.5, dr=.05, t=1` нельзя трактовать как физический массовый обмен. Формулы null flux правильны по непрерывной системе, однако ledger строится из конечных разностей `m_t,m_x`, глобальной нормировки и фиксированного конечного края. Он максимизируется по всей сетке, а место максимума не сохраняется в результате. Численная причина не локализована. #41F ledger и double-null ledger также являются различными диагностическими величинами и не должны сравниваться как один тест.
6. Кодовый audit не обнаружил подтверждённого непрерывного уравнения polar backend, расходящегося с massless canonical action при определениях из исходника. Проверенный локальный дефект — только пояснительная строка #40 для `K^r_r`, пропускающая `sqrt(f)`; вычисление и соседнее momentum constraint этот фактор содержат.

**Физический вердикт:** непрерывные уравнения обоих backend-ов согласуются с одним massless canonical action при заявленных определениях. Численная согласованность обратной связи не доказана: polar/FV constraints и ledger имеют прежние ошибки, а double-null constraints/balance велики. Из имеющихся результатов нельзя выбрать между начальными конечносеточными данными, центровым treatment, конечной границей, временным/пространственным усечением и оценкой ledger. Double-null формулы независимы по gauge и алгебре, но solver пока не является валидированным численным reference.

## 1. Что в проекте фактически задано

### 1.1 Действие и конвенции

В `experiment_40.py:381–389` (файл не включён в публичный экспорт) проект указывает `kappa=8 pi G`, `c=1`, действие

\[
S=\int d^4x\sqrt{-g}\left[\frac{R}{2\kappa}
-\frac12 g^{ab}\partial_a\phi\partial_b\phi-V(\phi)\right],
\]

с signature `(-,+,+,+)` (это прямо читается из `g_tt<0`, `g_rr>0`). Принят стандартный знак кривизны, при котором вариация даёт `G_ab=kappa T_ab`; сам код отдельно Riemann-sign convention не параметризует. Все численные величины в экспериментальном ядре dimensionless: `mu*r`, `mu*t`, `phi/M_Pl`, `V/(mu² M_Pl²)`; конкретная `mu` и перевод в SI не заданы. `kappa=.1` — выбранная dimensionless coupling данного вычислителя, а не измеренная физическая константа.

Полевая часть действия без неоднозначности даёт

\[
T_{ab}=\nabla_a\phi\nabla_b\phi
-g_{ab}\left[\frac12(\nabla\phi)^2+V\right],
\qquad \Box\phi-V_{,\phi}=0.
\]

Проверка знака из вариации: `delta S_phi = integral sqrt(-g) (Box phi - V') delta phi`; следовательно code source term для положительного `V'` должен входить в уравнение для ускорения со знаком минус.

Для стандартной ковариантной производной

\[
\nabla^aT_{ab}=(\Box\phi-V_{,\phi})\nabla_b\phi.
\]

Это непосредственное следствие симметрии вторых производных: члены `nabla^a phi nabla_a nabla_b phi` и `1/2 nabla_b[(nabla phi)^2]` сокращаются. На shell выполнении Klein–Gordon уравнения `nabla^a T_ab=0`; это условие Bianchi, а не отдельный выбор энергетического потока.

Потенциал, фактически заданный `potential()` / `potential_prime()` (`experiment_40.py:104–112`),

\[
V=\frac\lambda4\phi^2(\phi-1)^2+\delta(3\phi^2-2\phi^3),
\quad
V'=\frac\lambda2\phi(\phi-1)(2\phi-1)+\delta(6\phi-6\phi^2).
\]

Параметры #40 по умолчанию `lambda=1`, `delta=.02`; отчёт о тестах коллапса разрешил exact massless reduction только при `(lambda,delta)=(0,0)`. При нём `V=V'=0` для всякого `phi`; это общий physical model для отчёт о тестах коллапса и отчёт о double-null проверке. Double-null backend не содержит потенциала. Ни H2, ни дочерняя геометрия в настоящем аудите не участвуют.

### 1.2 Полевая нормировка и единицы

`experiment_40.py:385–389` документирует `Pi=(a/alpha) partial_t phi`, `Phi=partial_r phi`, `f=a^-2=1-2m/r`, `rho=f(Pi²+Phi²)/2+V`; реализация в `_metric` (`:135–170`) соответствует этому. Здесь `a` — радиальный metric coefficient `g_rr^(1/2)`, не космологический масштабный фактор. `m` имеет единицы длины в `c=1` и совпадает с `G M_phys` (если `G` явно восстановить). В double-null коде `Omega` безразмерен относительно тех же dimensionless координат, а `r` — dimensionless areal radius.

Геометрический множитель `4 pi` не пропущен: при `kappa=8 pi G` коэффициент `kappa/2=4 pi G`. Площадь сферы в mass equation уже свернулась с `G` и stress density в Misner–Sharp relation; повторно умножать эту формулу на `4 pi` нельзя. В FV recurrence используется `Delta V=(r_out³-r_in³)/3`, то есть радиальный интеграл `r²dr` без повторного множителя `4 pi` (`coupled_fv.py:140–151`). В #40 `4 pi r²` отдельно появляется только в весе объемной доли кандидата (`experiment_40.py:167`), не в уравнении `m_r`.

## 2. Независимое выведение polar-areal системы

### 2.1 Метрика и Misner–Sharp mass

Реальная source note (`experiment_40.py:385`) и код `_metric` задают

\[
ds^2=-\alpha(t,r)^2dt^2+f(t,r)^{-1}dr^2+r^2d\Omega^2,
\qquad f=1-\frac{2m(t,r)}r.
\]

Таким образом, вариант `g_tt=-alpha² f`, приведённый в задании как возможная конвенция, **не** используется. У radial null rays

\[
0=-\alpha^2dt^2+f^{-1}dr^2
\quad\Longrightarrow\quad \frac{dr}{dt}=\pm\alpha\sqrt f.
\]

Вне горизонта Misner–Sharp identity `g^{ab} partial_a r partial_b r=f` даёт `m=r(1-f)/2`.

### 2.2 Энергия, давление и Hamiltonian constraint

Нормаль к `t=const` срезу `n^a=alpha^{-1}(partial_t)^a`, outward unit radial vector `s^a=sqrt(f)(partial_r)^a`. Нормальная производная поля равна `n(phi)=phi_t/alpha`. Кодовая переменная имеет дополнительный пространственный множитель `a=f⁻¹/²`:

\[
\Pi_{code}=a\,n^a\nabla_a\phi=\frac{a\phi_t}{\alpha}
=\frac{\phi_t}{\alpha\sqrt f}.
\]

Следовательно `Pi_code` — не сама нормальная производная, а допустимая радиально перенормированная momentum variable. Из метрики следуют радиальные нулевые скорости `dr/dt=±alpha sqrt(f)`.

Для указанной в source-note метрики стандартный стресс в orthonormal frame `n,s` должен быть

\[
\rho=T_{ab}n^an^b=\frac12\left[(\phi_t/\alpha)^2+f\Phi^2\right]+V,
\quad
p_r=\frac12\left[(\phi_t/\alpha)^2+f\Phi^2\right]-V.
\]

Но code `_metric` использует

\[
\rho_{code}=\frac f2(\Pi_{code}^2+\Phi^2)+V,
\quad p_{r,code}=\frac f2(\Pi_{code}^2+\Phi^2)-V,
\quad \Pi_{code}=\frac{\phi_t}{\alpha\sqrt f}.
\]

Подстановка `Pi_code` даёт kinetic density `phi_t²/(2 alpha²)+f Phi²/2`, то есть стандартную энергию для этой метрики. `rho_code` и `p_r,code` согласованы.

Для сферического Misner–Sharp mass field equation из `G_t^t=kappa T_t^t`:

\[
m_r=\frac\kappa2r^2\rho.
\]

Установка code density даёт

\[
m_r=\frac\kappa4r^2\left[\frac{\phi_t^2}{\alpha^2}+f\Phi^2\right]
 +\frac\kappa2r^2 V.
\]

В massless case `m_r=(kappa/4)r² f(Pi_code²+Phi²)`, что кодовый комментарий иногда сокращённо пишет как `(kappa/2)r²rho` (`experiment_40.py:139–149`). #40 recurrence имеет `A=c r²(.5s+V)`, `B=c r s`, `c=kappa/2`; это trapezoidal implicit discretization соответствующего linear-in-mass constraint (`:138–148`).

**Классификация:** mass constraint coefficient/4π — PASS. `Pi` — не обычная нормальная производная поля `n^a∂_aφ`; это радиально перенормированная каноническая переменная `a φ_t/alpha`. Такая переменная допустима и именно для неё записаны кодовые first-order уравнения. Плотность и радиальный поток корректны с этим определением.

### 2.3 Lapse, KG и mass flux: проверка на согласованность

Для lapse исходя из `G_r^r` в указанной метрике стандартная формула:

\[
\partial_r\ln\alpha=\frac{m/r^2+(\kappa/2)r p_r}{f}.
\]

Она совпадает с `_metric` (`experiment_40.py:151–153`). Интеграция задаёт `alpha(r_max)=1`, т.е. конечную-boundary normalisation, а не строго ADM infinity.

Для метрики `ds²=−alpha²dt²+a²dr²+r²dΩ²`, `a=f⁻¹/²`, детерминант `sqrt(-g)=alpha*a*r² sin(theta)`. Определим кодовую переменную `Pi=a phi_t/alpha` и `Phi=phi_r`; нормальная производная равна `phi_t/alpha=Pi/a`. Из `Box(phi)-V'=0`, после умножения на `alpha*a`, непосредственно следует

\[
\phi_t=\frac{\alpha}{a}\Pi,\qquad
\Pi_t=\frac1{r^2}\partial_r\!\left(r^2\frac{\alpha}{a}\Phi\right)-\alpha a V'.
\]

Проверка факторов: временная часть KG равна `−∂t(a r² phi_t/alpha)`; пространственная — `∂r(r² alpha/a phi_r)`; деление на `r²` и потенциал дают формулу выше. Поэтому кодовые `alpha/a=alpha sqrt(f)` и `alpha*a=alpha/sqrt(f)` в `_rhs` (`experiment_40.py:205–213`) **согласованы**, включая зависимость `a(t,r)`: производная уже применяется к составной переменной `Pi=a phi_t/alpha`, и отдельный член `f_t Pi` добавлять нельзя.

Так как `rho=1/2[(phi_t/alpha)²+f Phi²]+V=f(Pi²+Phi²)/2+V`, запись в `_metric` корректна. Для `p_r` кинетические градиенты складываются, а потенциал меняет знак: `p_r=f(Pi²+Phi²)/2−V`. Никакого второго polar metric convention для объяснения расхождений вводить не требуется.

Независимая проверка через `G^r{}_t` и Misner–Sharp identity для той же метрики даёт

\[
m_t=\frac\kappa2 r^2 f\,\phi_t\Phi
=\frac\kappa2r^2\alpha f^{3/2}\Pi_{code}\Phi,
\]

Эта формула совпадает с кодовой `mass_flux_rhs` (`experiment_40.py:245–250`, `coupled_fv.py:184–190`). Для исходящего асимптотического импульса `phi≈F(t−r)/r`, `phi_t≈−phi_r`, поэтому `Pi Phi<0` и `m_t<0`: масса внутри фиксированного areal radius уменьшается. Знак и коэффициент кода совпадают с этим пределом.

Сферическая mass identity записывается как `∂_a m=(κ/2)r²(T_a{}^b−δ_a{}^b T_c{}^c)∂_b r` на двумерной орбите. Для `a=t` и `r_t=0` trace-члены исчезают, `∂_r r=1`, поэтому `m_t=(κ/2)r²T_t{}^r`; для скаляра `T_t{}^r=g^{rr}phi_t Phi=f phi_t Phi`. Это показывает происхождение знака, коэффициента и отсутствие дополнительного множителя `4π`.

**Вывод проверки KG:** полярная метрика, кодовая переменная `Pi`, уравнение Клейна–Гордона, плотность, lapse и поток массы составляют согласованный непрерывный набор. Это снимает гипотезу о доказанной ошибке непрерывного KG; оно не доказывает корректность конечносеточной реализации или эволюции constraints.

### 2.4 ADM momentum and expansions

При кодовой конвенции `K_ij=−(1/2) L_n h_ij`, нулевом shift и `h_rr=f⁻¹`,

\[
K^r{}_r=-\frac{m_t}{\alpha r f}.
\]

Подстановка потока массы даёт `K^r_r=−(kappa/2)r sqrt(f) Pi_code Phi`. Это совпадает с реализациями `_metric:154` и `_metric_from_mass:189–193`. `K^theta_theta=0`, поскольку при нулевом shift сферы координат имеют постоянный areal radius. Поэтому расширения по единичным направлениям `n±s` равны `theta_±=±2 sqrt(f)/r`, как и в `_metric:161–163`.

Поскольку единственная ненулевая компонента смешанного `K^i{}_j` — радиальная, trace `K=K^r{}_r`; полярная momentum constraint сводится к `K+(κ/2)r sqrt(f) Pi Phi=0`, что совпадает с кодовым `momentum` и пояснительной записью через `1/a` в `:158–159,389`. Это дополнительная проверка коэффициента `K^r_r`.

Текстовая заметка `experiment_40.py:387` указывает `K^r_r=−(kappa/2)r Pi Phi` **без** `sqrt(f)`; соседняя строка 388 уже содержит этот множитель в momentum constraint через `/a=sqrt(f)`. Это пропуск в документации, а не в исполняемой формуле.

## 3. Finite-volume backend: алгебра и дискретная согласованность

`coupled_fv.py:_slice` (`:126–209`) задаёт объём сферического слоя `V_i=(r_out³−r_in³)/3`, `c=kappa/2`, градиенты на гранях, `Pi` в центрах ячеек и приближённое среднее по слою

\[
\Delta m_i\simeq c V_i\left[\frac12 f_i(Pi_i^2+\langle\Phi^2\rangle_i)+V(phi_i)\right].
\]

При `f_i≈1−(m_i+m_{i+1})/r_i` решение этого соотношения относительно `m_{i+1}` даёт `aa`, `rhs`, `mass[i+1]=rhs/(1+aa)` в `:144–152`. Алгебраически это среднее приближение по сферическому слою к той же радиальной гамильтоновой связи; это не точная квадратура континуума.

Поток Эйнштейна через грань `c r_face² alpha_face f_face^(3/2) Pi_face Phi_face` (`:184–190`) имеет ожидаемые множители для polar mass flux. Дивергенция scalar-потока на гранях `(r_out² F_out−r_in²F_in)/V_i` (`:217–220`) аппроксимирует `r⁻² partial_r(r² alpha sqrt(f) Phi)`. Исходящее асимптотическое условие `Phi_R=−Pi_R` (`:137`) следует из характеристического производного `phi_t+alpha sqrt(f)phi_r=alpha sqrt(f)(Pi+Phi)`. Поэтому **знак и непрерывные множители граничного потока не первая подозреваемая причина**.

Независимое восстановление по Гауссу `hres` (`:192–208`) не повторяет точно midpoint-рекурсию по слоям, что и делает его полезной проверкой. `_constraint_mass_tangent` дифференцирует рекурсию (`:224–256`) и сравнивает её касательный вектор на стадии с независимым потоком Эйнштейна. Эти конечносеточные остатки могут быть ненулевыми; оценивать их следует по сходимости к континууму, а порог точного совпадения `1e−3` не является теоремой.

FV использует ту же пару переменных: `dphi=alpha sqrt(f)Pi`, `dpi=div−alpha/sqrt(f)V'` (`:212–220`). Поскольку это и есть KG-система, выведенная выше из действия, continuum-ошибка уравнения здесь не установлена. Пространственное уточнение должно проверять одновременно независимую Hamiltonian constraint, momentum constraint, mass flux и ожидаемый порядок дискретизации.

В ветви `mass_flux` (`run_coupled_fv:260–335`) масса `m` развивается тем же потоком на гранях, поэтому собственный mass ledger закрывается по построению. Независимо восстановленные Hamiltonian и momentum constraints всё равно обязательны: одно лишь замыкание ledger ничего не говорит о выполнении связей.

## 4. Независимый вывод double-null системы

### 4.1 Сведение уравнений с явными null-множителями

For `V=0` define

\[
ds^2=-\Omega^2 du\,dv+r^2d\Omega^2,
\quad g_{uv}=g_{vu}=-\Omega^2/2,
\quad \sigma=\ln\Omega.
\]

Then `g^{uv}=−2/Omega²`, `Box_2 r=−4 r_uv/Omega²`, `(nabla r)²=−4r_u r_v/Omega²`, and `Gamma^u_uu=2sigma_u`, `Gamma^v_vv=2sigma_v`. The warped-product base Einstein tensor is

\[
G_{ab}=-\frac2r\nabla_a\nabla_b r
+h_{ab}\left[\frac2r\Box_h r+\frac{(\nabla r)^2-1}{r^2}\right].
\]

Для безмассового поля `T_uu=phi_u²`, `T_vv=phi_v²`, `T_uv=0`. Подстановка в уравнения даёт

\[
r_{uu}-2\sigma_u r_u=-\frac\kappa2 r\phi_u^2,
\quad r_{vv}-2\sigma_v r_v=-\frac\kappa2 r\phi_v^2,
\]

\[
G_{uv}=\frac{2r_{uv}}r+\frac{2r_ur_v}{r^2}+\frac{\Omega^2}{2r^2}=0
\ \Longrightarrow\ r_{uv}=-\frac{r_ur_v}{r}-\frac{\Omega^2}{4r}.
\]

Уравнение Эйнштейна со снятым следом имеет вид `R_ab=kappa phi_a phi_b`; для warped-product разложения `R_uv=−2 sigma_uv−2 r_uv/r` получаем

\[
\sigma_{uv}=\frac{r_ur_v}{r^2}+\frac{\Omega^2}{4r^2}-\frac\kappa2\phi_u\phi_v.
\]

Волновое уравнение поля:

\[
\phi_{uv}=-\frac{r_u\phi_v+r_v\phi_u}{r}.
\]

При `u=t−x`, `v=t+x` имеем `partial_u=(partial_t−partial_x)/2`, `partial_v=(partial_t+partial_x)/2`, `partial_uv=(partial_tt−partial_xx)/4`. Отсюда

\[
r_{tt}=r_{xx}+\frac{r_x^2-r_t^2-\Omega^2}{r},
\]

\[
\sigma_{tt}=\sigma_{xx}+\frac{r_t^2-r_x^2+\Omega^2}{r^2}
-\frac\kappa2(\phi_t^2-\phi_x^2),
\quad
\phi_{tt}=\phi_{xx}+\frac{2(r_x\phi_x-r_t\phi_t)}r.
\]

Это в точности формулы из `double_null_reference.py:_rhs`, строки 245–252. **Коэффициенты и знаки: PASS по ручному выводу** для указанной null-метрики и безмассового действия.

### 4.2 Double-null масса, поток и расширения

Тот же инвариант даёт

\[
m=\frac r2\left(1+\frac{4r_ur_v}{\Omega^2}\right)
=\frac r2\left(1+\frac{r_t^2-r_x^2}{\Omega^2}\right).
\]

Дифференцируем определение и устраняем `r_uu/r_vv` через null constraints, а `r_uv` — через смешанное уравнение Эйнштейна. Члены с `r_u r_v` сокращаются, остаётся

\[
m_u=-\kappa r^2\Omega^{-2}r_v\phi_u^2,
\qquad
m_v=-\kappa r^2\Omega^{-2}r_u\phi_v^2.
\]

`misner_sharp_mass()` (`:95–103`) и предсказанные `mu_pred,mv_pred` в `_diagnostics()` (`:273–304`) используют эти коэффициенты и знаки. Численные `m_u/m_v` вычисляются конечными разностями геометрической массы, а не берутся из самих уравнений эволюции (`:374–397`). Вердикт формул: **PASS**; дискретная точность одним порогом не установлена.

Расширения по будущим null-касательным равны `theta_out=2r_v/r`, `theta_in=2r_u/r` (с точностью до положительной нормировки). Строки `:361–365` совпадают с этим. Знаки определяют локальный trapped/anti-trapped статус, но не глобальный event horizon. При `V=0` скаляр кривизны `R=kappa (nabla phi)^2=kappa Omega^-2(phi_x²−phi_t²)`; поле `ricci_scalar` использует это тождество (`:358–360`).

### 4.3 Центр, данные Коши и координатная калибровка

Для регулярного центра в conformal radial chart гладкие сферические скаляры чётны, `r` нечётна, `r(0)=0`, `Omega(0)=r_x(0)`, так что собственное радиальное расстояние и areal radius совпадают в первом порядке. Код задаёт `sigma(0)=ln(r_x(0))`, `r_t(0)=0`, чётность поля и Taylor-предел ускорений в центре (`_project_center_and_outer:213–223`, `_rhs:254–260`). Поэтому остаток центральной метрики и `phi_x(0)` навязываются projection/parity и не являются независимыми тестами точности в центре. L2-остаток constraints намеренно исключает первые две и последние две точки (`_diagnostics:286`).

Два независимых характеристических импульса на пересекающихся null-ветвях не задаются. `_initial_slice:159–210` строит одно time-symmetric пространственноподобное Cauchy-сечение: интегрирует `m_r=kappa/4 r² f phi_r²` и полярную lapse с `alpha(r_max)=1`, затем задаёт

\[
\frac{dx}{dr}=\frac1{\alpha\sqrt f},\qquad r_x=\alpha\sqrt f,
\quad \sigma(0,x)=\ln\alpha(r(x)),
\quad r_t=\sigma_t=\phi_t=0.
\]

At the initial slice this is a valid radial coordinate transformation: the polar spatial metric becomes `g_xx=alpha²`, hence `ds²=alpha²(−dt²+dx²)` there. It maps the initial areal Cauchy geometry and time-symmetric data; no overdetermined independent ingoing/outgoing pulses are prescribed.

**Различие будущих калибровок:** если полярное решение продолжается с меняющимися `alpha(t,r),f(t,r)`, фиксированная функция `x(r)` в общем случае не сохраняет `g_tt=−g_xx`. Настоящие null-функции должны решать уравнения эйконала

\[
-\alpha^{-2}(u_t)^2+f(u_r)^2=0,
\qquad -\alpha^{-2}(v_t)^2+f(v_r)^2=0.
\]

Упрощённый выбор `u=t−x(r), v=t+x(r)` требует `x_t=0`, `x_r=1/(alpha sqrt(f))`; при динамической метрике это условие обычно справедливо только на начальном сечении. **Исходник отчёт о double-null проверке не утверждает, что такой фиксированной картой преобразуется будущее решение #40.** Он отдельно эволюционирует conformal-метрику `Omega²(−dt²+dx²)` уравнениями double-null. Поэтому начальная optical-карта — допустимый выбор координат для данных, но не доказанное преобразование будущего решения. Для прямого сравнения будущих событий нужно интегрировать координатное преобразование или сопоставлять инвариантные события.

На конечном внешнем координатном крае фиксируются `r,sigma,phi` и их скорости (`_project_center_and_outer:213–223`). Это временеподобная искусственная граница, а не бесконечность; она отличается от фиксированного areal-края #40 и исходящей characteristic-грани FV. До причинного контакта можно сравнивать внутренние области; после контакта нужно отдельно проверять отражение и null constraints. На вакуумном сечении с регулярным центром `M=0`; `schwarzschild_exterior_mass()` проверяет лишь аналитическое тождество статического внешнего решения при `r>2M`, а не эволюцию Schwarzschild-геометрии с регулярным центром.

## 5. Сверка исходников с выведенными формулами

| Formula / choice | Independent result | Source trace | Status |
|---|---|---|---|
| `kappa=8 pi G`, `V` and `V'` | As printed in source; `(0,0)` is exact massless limit | `experiment_40.py:21–37,104–112`; massless harness `collapse_benchmark.py:16–78` | **PASS** |
| Metric `−alpha²dt²+f⁻¹dr²` | This gives radial null speed `±alpha sqrt(f)` | `experiment_40.py:385` | **PASS as declared metric** |
| `Pi=a phi_t/alpha`, `a=f⁻¹/²`, and energy density | `Pi` is a radially rescaled momentum; `rho=1/2[(phi_t/alpha)^2+fPhi²]+V` | `experiment_40.py:135–159,386` | **PASS** |
| Scalar time evolution | KG gives `phi_t=alpha Pi/a`, `Pi_t=r⁻²∂r(r² alpha Phi/a)−alpha*a*V'` | `experiment_40.py:205–213,388` | **PASS** |
| Lapse equation | Standard `G_r^r` polar lapse formula for declared metric | `experiment_40.py:151–153,387` | **PASS under declared metric** |
| Mass flux | `m_t=(kappa/2)r² f phi_t Phi=(kappa/2)r²alpha f^(3/2) Pi_code Phi` | `experiment_40.py:245–250,305`; `coupled_fv.py:184–190` | **PASS under declared metric/Pi** |
| ADM `K^r_r` | `−(kappa/2)r sqrt(f) Pi Phi` | implementation `experiment_40.py:154,190–193`; note `:387` omits `sqrt(f)` | **Code PASS; source-note PROVEN typo** |
| FV shell Hamiltonian recurrence | Algebraic midpoint shell-volume discretization of `m_r` | `coupled_fv.py:126–155` | **Continuum form PASS; truncation is approximate** |
| FV independent H quadrature | Separate reconstruction/quadrature from recurrence | `coupled_fv.py:192–208` | **Useful independent diagnostic; not exact** |
| Double-null mixed Einstein/KG equations | Match derived equations with stated `g_uv` | `double_null_reference.py:226–261` | **PASS algebraically for V=0** |
| Double-null null constraints/flux | Correct continuous formulas; finite-difference realization remains resolution dependent | `double_null_reference.py:273–304,374–397` | **Formula PASS; numerical accuracy UNVERIFIED** |
| Center residuals | Metric/symmetry conditions are imposed; null constraint excludes center cells | `double_null_reference.py:213–223,254–260,286` | **Not independent validation** |
| Cross-gauge initial `m(r)` | Same continuum data but disparate finite-grid quadratures | initial functions `experiment_40.py:135–149`, `double_null_reference.py:159–210` | **PROVEN finite-grid mismatch; source is discretization candidate** |

The continuous polar equations close under the stated metric and code variable definitions. The checked documentation defect is narrower: note line 387 omits `sqrt(f)` from `K^r_r`, while the implementation and momentum constraint include it. No source code is changed in this audit.

**Сверка с тестами (файлы просмотрены, тесты в этом задании не запускались):** `tests/test_experiment_40.py` проверяет локальные свойства #40; `tests/test_coupled_fv.py` включает Minkowski-контроль, наблюдаемое уточнение H и тот факт, что mass-flux-эволюция закрывает ledger при одновременно плохих Hamiltonian/momentum остатках; `tests/test_collapse_benchmark.py` проверяет massless reduction/Minkowski и слабый импульс для обоих backend-ов; `tests/test_double_null_reference.py` проверяет аналитическую постоянную Schwarzschild-массу, вакуум Minkowski, слабое сравнение с плоской волной и записанные constraint/ledger статусы. Эти тесты подтверждают собственные регрессии/контроли, но не образуют независимого доказательства гравитационной обратной связи: сильные null constraints и ledger остаются большими, и ни один тест не сравнивает обе динамические системы на общем нелинейном событии.

## 6. Read-only initial-slice probe: the 2.07% comparison

Only the initial-data helpers were evaluated: #40 `_metric()` and отчёт о double-null проверке `_initial_slice()`/`_d1()`/`misner_sharp_mass()`. No time evolution ran. Common continuum data: `A=.5`, `R0=.5`, `kappa=.1`, `r_max=16`, `Pi=0`, exact `V=0`; polar nodes and null-grid interpolation represent the same smooth profile. The finite grids do not use the same radial derivative or constraint quadrature.

| dr | Double-null radial RK4 constraint mass `m_ODE(R)` | Double-null geometric mass `m_geom(R)` | #40 polar recurrence `m_PA(R)` | `(m_geom−m_PA)/m_PA` | max radial-profile difference / `m_PA(R)` |
|---:|---:|---:|---:|---:|---:|
| .20 | `.00351155135` | `.00351104438` | `.00283562527` | `23.82%` | `33.55%` |
| .10 | `.00352094807` | `.00352083989` | `.00330451696` | `6.546%` | `8.796%` |
| .05 | `.00352094096` | `.00352092234` | `.00346313417` | `1.669%` | `2.066%` |

The difference `m_geom−m_ODE` inside the double-null initialization is only `−.0144%`, `−.00307%`, `−.000529%` of `m_ODE`; the transformed geometric identity closely preserves the mass integrated on the areal initial slice. The much larger polar-vs-RK4 difference is therefore mostly the different discrete source integration/gradient (`np.gradient(phi)` plus implicit trapezoid in #40 versus analytic `phi_r` plus radial RK4 in отчёт о double-null проверке), with additional interpolation/`r_x` finite-difference error. The finite-difference scalar-gradient RMS mismatch over `r<2` was `.11194`, `.02750`, `.00702` for `dr=.2,.1,.05`, also roughly second-order. This is a **measured finite-grid initial-data mismatch**; exact share of each quadrature component needs an intentionally designed refinement test. It is not evidence of different κ, physical normalization, lapse, or continuum initial surface.

The probe also verifies the initial optical mapping itself: its construction uses the same `alpha(r_max)=1` and `dx/dr=(alpha sqrt(f))⁻¹`. There is no sign/factor-of-`f` error in that initial transformation as coded. The 2.07% discrepancy arises when comparing two discretized mass reconstructions, not when comparing their continuum coordinate labels.

## 7. Residuals, norms, and hypotheses for the old failures

### 7.1 What can be concluded now

- #41F mass/constraint discrepancies are not explained by a demonstrated polar KG continuum mismatch: with `Pi=a phi_t/alpha`, the code RHS is the equation derived from the action. The reported failures remain numerical convergence/consistency issues or unlocalized diagnostics.
- отчёт о double-null проверке double-null massless equations, as continuous equations, are consistent. Its `Cuu/Cvv=.0728/.1047` and ledger `11.18%` at `dr=.05` indicate that this implementation has not yet propagated the derived constraints/balance accurately. They do not imply an error in GR or an H1 no-go.
- Double-null relative mass-flux norms differentiate geometric mass in both `t` and `x`; differentiation noise, center projection, boundary cells, null constraint violation, temporal truncation and normalization are all candidates. The saved output reports a maximum but not the location or component decomposition. No numerical cause is proven.
- The double-null 11.18% ledger is a max-over-coordinate-grid integral residual, not the outer ADM loss alone. The finite boundary is frozen, but at short `t=1` the pulse is initially central in a domain reaching `r=16`; boundary reflection is less likely to explain the central part, but the stored max location is unavailable, so boundary contribution is **not excluded**.
- The double-null initial `m` mismatch is not caused by the unrelated 1e−3 acceptance gate; the read-only comparison isolates a real finite-grid construction difference.

### 7.2 Why a universal `1e−3` raw gate is not meaningful

`experiment_40._summary()` compares raw maximum absolute residuals against fixed thresholds (`experiment_40.py:342–355`). Their dimensions are code-normalized, but their magnitude changes with amplitude, `kappa`, domain and grid; no truncation extrapolation is part of that gate.

FV uses several different normalizations (`coupled_fv.py:29–44,377–408`): H L2 relative to shell mass increments, momentum Linf relative to maximum expected K, local-flux L2 relative to numerical/analytic flux scale, ledger relative to initial total mass, stage tangent relative to max tangent/flux. They are dimensionless but represent different scale choices; near a vanishing source, a relative norm becomes ill-conditioned. A shared numerical value `1e−3` is a useful engineering screen only after documenting each denominator, active region, mesh order and source floor. It is not a universal physical accuracy threshold.

Double-null constraints use L2 residual divided by RMS of the three equation terms (`double_null_reference.py:286–295`), but cells 0–1 and the final two cells are excluded. At near-zero source, relative norms can be NaN or dominated by roundoff. Absolute norms and refinement orders are essential. The `11.18%` ledger denominator is total initial mass; the residual maximum is over x, so localized low-mass regions can dominate.

### 7.3 Prioritized diagnostic hypotheses (not repairs)

| Priority | Candidate | Evidence / discriminant | Status |
|---:|---|---|---|
| 1 | Different finite-grid initial constraint construction | Probe shows `np.gradient`/trapezoid vs analytic gradient/RK4/interpolation; compare against one common continuum initial-data reference. | **Measured; likely O(dr²), exact decomposition unverified** |
| 3 | Double-null near-center Taylor/projection and discrete null constraints | center is projected; constraint residual drops first two cells; compare independently evaluated constraints including first nonzero cells with analytic Taylor expansions under refinement. | **Hypothesis** |
| 4 | Double-null finite-boundary condition | `r,sigma,phi` are frozen on outer fixed-x boundary, not a radiative Einstein boundary; compare only before null domain of dependence reaches boundary and separately move it. | **Known limitation; contribution to t=1 max unlocalized** |
| 5 | Double-null mass-flux finite differences/temporal quadrature | `m_t,m_x` use `np.gradient`; predicted null fluxes use same saved fields; refine output stride/time step separately, report abs/relative norms and spatial location. | **Hypothesis** |
| 6 | FV shell-average H recurrence and scalar FV flux use different discrete product rules | Independently reconstructed H shows convergence, but exact stage equality is not required unless scheme has that identity. First make continuum convention consistent; then compare observed orders. | **Known discretization difference; not yet isolated as root cause** |
| 7 | κ, `4 pi`, outgoing sign or mass-flux coefficient | Both implementations use same κ=.1; independently derived coefficient is κ/2; outgoing wave has `Pi Phi<0` and `m_t<0`. | **Not leading suspect; formulas PASS** |
| 8 | Coordinate-array/time normalization explains mass mismatch | Misner–Sharp mass is scalar invariant; initial slices were matched by areal radius and time symmetry. Lapse normalization agrees at outer endpoint initially. | **Rejected as explanation of the initial mass difference** |
| 9 | Polar KG continuum mismatch | With `Pi=a phi_t/alpha`, code RHS equals the independently derived KG equation. | **Rejected by this audit** |

## 8. Предлагаемые минимальные различающие проверки (не выполнялись)

1. **Polar continuum closure, zero code changes:** independently derive `G_t^t`, `G_r^r`, `G_tr`, KG and ADM momentum with the source definitions; verify the code's `m_r`, lapse, `m_t`, Pi RHS and `K^r_r` as a complete tuple. This audit finds the implemented formulas match, with one explanatory-note typo. An independent symbolic tensor calculation is optional and was not run.
2. **Проверка начального источника в первом порядке по κ:** зафиксировать одинаковые аналитические `phi(r)`, `Pi=0`, `V=0`; независимо вычислить квадратурой `m_r=(κ/4)r² phi_r²+O(κ²)`. При `κ→0` величина `m/κ` должна сходиться к `1/4 integral_0^R r² phi_r² dr`; сравнить polar-рекурсию, radial RK4 double-null и геометрическую массу из `r_x`. Разность коэффициента, не исчезающая при `dr→0`, опровергнет нормировку; расхождение порядка `O(dr²)` укажет на квадратурную ошибку.
3. **Одно физическое начальное сечение Коши:** зафиксировать `A,R0,R,κ`, решить непрерывную радиальную связь с высокой точностью и спроецировать одни и те же данные `r,phi,K=0` на обе сетки. Сравнить `m(r)`, `f(r)`, `alpha(r)` общей квадратурой по areal radius. По одному менять только оператор производной или квадратуру. Так измеренное расхождение отчёт о double-null проверке в 2.07% можно разложить на ошибки градиента источника, radial solver и интерполяции `x`.
4. **Независимые null constraints:** только из сохранённых `(r,sigma,phi)` вычислить `r_uu,r_vv` stencil-производными, отличными от stencil эволюции; включить Taylor-предел центра и отдельно привести нормы в центре, bulk и внешней области. Проверить порядок минимум на трёх сетках и независимо уточнить `dt`. Плато остатка после согласования непрерывных соглашений укажет на дефект эволюции/переноса связи; убывание остатка количественно покажет truncation error.
5. **Проверка внешнего Schwarzschild-решения:** при `f=1−2M/r`, tortoise `x` с `dr/dx=f`, `Omega²=f` проверить `m=M`, `R=0`, `Cuu=Cvv=0` теми же дискретными производными при `r>2M`. Это проверит факторы null-метрики в искривлённой геометрии, не делая утверждения о регулярном центре. Нынешний helper проверяет лишь алгебраическое `m=M`, но не эволюционные уравнения.
6. **Причинная проверка границы:** определить характеристическую область до сравнения с конечной внешней границей. Сначала взять один исходящий импульс в большом домене; затем отдельно довести его до границы и измерить отражённые/входящие характеристики. Не считать замороженную metric-границу нулевой бесконечностью.

### Stop/go ordering

**Порядок диагностики:** не принимать решение о solver repair только по интегральному ledger. Непрерывные polar формулы теперь замкнуты алгебраически; следующим различающим тестом был бы общий reference для начального `m/κ`, затем независимые пространственно-временные оценки double-null constraints/flux по областям и сеткам. После этого отдельно исследовать boundary и temporal errors. До валидации обоих численных решений не переходить к near-critical search.

## 9. Часы H1 и происхождение сведений

Исследовательское соотношение `dt_d/dT_p=K0 a_d(t_d)` здесь имеет статус **NOT_COMPARABLE**: не заданы ни дочерний FLRW `a_d`, ни мировая линия родительского наблюдателя, ни связь двух многообразий. Ни `mu*t`, ни координатное время, нормированное lapse, не являются родительскими часами. Этот аудит касается только стандартных классических уравнений Эйнштейна–скаляра.

**Использованные источники:** действие, соглашения и алгоритмы прочитаны из локальных файлов проекта; прежние численные значения взяты из локальных отчётов #41F/отчёт о тестах коллапса/отчёт о double-null проверке. Тензорные формулы и переходы между координатами в этой записке выведены заново вручную из действия и тождеств warped-product, а не приняты на веру из прошлых отчётов. Сеть/API, новые пакеты, SymPy, тесты, эволюция и модификации solver не использовались. Единственный выполненный расчёт — read-only probe начальных сечений из раздела 6, использующий существующие `_metric`, `_initial_slice` и конечные разности для геометрической массы.
