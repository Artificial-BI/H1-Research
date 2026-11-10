# GH-IBVP для SO(4)-симметричной 5D Einstein–scalar системы: доказанная область и блокер

## Краткий вердикт

**Вывод не даёт GO для полного solver.** Получены три раздельных результата:

1. **CONTINUUM IVP, linearized/frozen Minkowski:** полный 5D generalized-harmonic Einstein–scalar system имеет wave principal part; первый порядок обладает вещественными характеристическими скоростями и положительным symmetrizer при `α>0`, положительной spatial metric и derivative-independent smooth gauge source. Constraint subsidiary system — волновая система. Это доказывает стандартную локальную гиперболическую часть, не глобальную регулярность/нелинейную устойчивость.
2. **SO(4) центр и frozen outer boundary:** регулярный центр корректно трактуется как smooth Cartesian parity subspace, не как физическая incoming boundary. Для линейного frozen-coefficient wave/constraint system дана maximally dissipative energy estimate и CPBC-условие на incoming subsidiary constraints. В полном Cartesian GH характерические условия должны затем быть переведены на SO(4)-симметричные metric modes.
3. **FULL nonlinear reduced IBVP и SBP-SAT Einstein solver:** не доказаны. В частности, старое orbit-space условие из `mathematical-equations.md`, если читать его как полный 5D GH, пропускает warp-factor член `−3∇^aR/R`; angular coordinate conditions нельзя просто отбросить. Точное преобразование nonlinear GH CPBC в reduced variables, boundary stability при variable coefficients/shift и совместная parity/SAT treatment на регулярном центре остаются открытым блокером.

Практическая рекомендация: сохранять GH как кандидат, но специфицировать его в полных 5D Cartesian координатах с SO(4)-эквивариантными полями; до code issue закрыть одну конкретную reduced CPBC/center energy proof. Ни одна формула ниже не объявляет численную реализацию корректной.

## 1. Действие, gauge и точная spherical geometry

Сигнатура `(-++++)`, `κ₅²=8πG₅`, действие

\[
S=\frac1{2κ_5^2}\int\sqrt{-g}(R-2Λ_5)d^5x
-\int\sqrt{-g}\left[\frac12(\nablaφ)^2+V(φ)\right]d^5x.
\]

Уравнения: `G_AB+Λ₅g_AB=κ₅²T_AB`, `T_AB=∂Aφ∂Bφ−g_AB[(∂φ)²/2+V]`, `□₅φ−V'=0`. Применяются только классические Levi–Civita геодезические. Ни ECSK, ни EGB в этот вывод не входят.

Полная GH-конвенция:

\[
C^A:=Γ^A-H^A=0,
\qquad Γ^A=g^{BC}Γ^A{}_{BC}.
\]

Рабочая координатная карта фиксирована Cartesian; gauge source выбирается с transformation law, делающим `C^A` constraint-вектором в covariant reduced equation. `H^A` — заданные smooth gauge-source functions от координат и/или metric, **без производных metric** в principal-part analysis. Если `H` зависит от metric, вариация algebraic `H(g,x)` меняет lower-order terms; если от metric derivatives — principal symbol/гиперболичность надо вывести заново. Для constraint propagation damping пока исключён; он добавляет lower-order terms и требует своего subsidiary analysis.

Точная SO(4) orbit metric:

\[
ds_5^2=h_{ab}(x)dx^a dx^b+R(x)^2dΩ_3^2,
\quad a,b\in\{t,r\},\quad X=h^{ab}R_aR_b.
\]

В этой записи 5D contracted Christoffel base-компонента равна

\[
Γ^a_{(5)}=h^{bc}Γ^a{}_{bc}[h]-\frac3R\nabla^aR.
\]

Следовательно, **полный 5D GH condition** на base имеет вид

\[
C^a=h^{bc}Γ^a{}_{bc}[h]-3R^{-1}\nabla^aR-H^a=0,
\]

а не `Γ[h]^a−H^a=0`. Для углов в обычной координатной карте `S³` contracted Christoffels unit-sphere тоже входят в `Γ^i_(5)`; соответствующий `H^i` должен быть согласован с выбранной angular/reference-coordinate connection. Нельзя объявить angular GH constraints тождественно нулевыми, просто потому что `h_ab,R,φ` не зависят от углов.

### Выбранное устранение angular-coordinate ambiguity

Использовать полные 5D Cartesian coordinates `(t,x^I)`, `r²=δ_IJx^Ix^J`, `I=1,…,4`, и SO(4)-эквивариантные metric components:

\[
g_{00}=A(t,r^2),\quad g_{0I}=B(t,r^2)x_I,\quad
g_{IJ}=C(t,r^2)δ_{IJ}+D(t,r^2)x_Ix_J,\quad φ=φ(t,r^2).
\]

Эквивалентная сферическая запись:

\[
ds²=A dt²+2rB\,dt\,dr+(C+r²D)dr²+r²C\,dΩ_3²,
\quad R=r\sqrt C.
\]

`A,B,C,D` гладкие по `r²`, `A<0`, positive spatial metric; shift cross term `rB` нечётен. GH source выбирается эквивариантно: `H^0=H^0(t,r²)`, `H^I=x^I H^1(t,r²)`. В Cartesian form все пять `C_A` сохранены: `C_0` и радиальная компонента `C_r=x^I C_I/r` независимы в SO(4) секторе; angular components не дополнительные modes, но их связь уже включена в радиальную tensor ansatz. При вычислении в spherical chart надо использовать full transformed connection/source, а не вычислять только `Γ[h]`.

## 2. Reduced field equations и principal part

Независимая covariant reduction (не gauge fixing) из предыдущей спецификации:

\[
R[h]_{ab}-3R^{-1}\nabla_a\nabla_bR
=κ_5²\left(φ_aφ_b+\frac23Vh_{ab}\right)+\frac23Λ_5h_{ab},
\]

\[
2(1-X)-R□_hR=\frac23(κ_5²V+Λ_5)R²,
\qquad □_hφ+3R^{-1}\nabla R\cdot\nabla φ-V'=0.
\]

Они точны как spherical field equations, но не являются сами по себе GH evolution system. GH reduction следует делать **из full 5D tensor equations** или доказанно эквивалентным reduced variational system; полные base/angular constraints должны совпадать с reduced equations above.

Для полного GH reduced Einstein equations удобно определить

\[
E_{AB}:=G_{AB}+Λ_5g_{AB}-κ_5²T_{AB}
-\nabla_{(A}C_{B)}+\frac12g_{AB}\nabla_CC^C=0.
\]

Принципиальная часть каждого metric компонента — `−½g^{CD}∂C∂D g_AB`; scalar equation имеет principal part `g^{CD}∂C∂Dφ`. Материальные нелинейности `∂φ∂φ` и `Γ∂φ` входят как lower-order terms в стандартной first-order reduction, поскольку `∂φ` — независимые variables.

На произвольной гладкой background principal symbol в нормальном направлении после ADM split имеет скорости

\[
λ_\pm=-β^n\pm α,
\qquad λ_0=-β^n,
\]

с кратностями волновых metric/scalar modes и advected derivative/metric modes. При fixed timelike outer boundary и `α>|β^n|` имеются входящие (`λ<0`) и исходящие (`λ>0`) световые характеристики. При Minkowski `α=1,β=0`: `λ=±1,0`.

### First-order principal block

Для каждого metric component (и отдельно scalar) `u` положить

\[
Π=-n^A∂_Au,
\qquad Φ_i=∂_iu,
\]

и использовать principal system в локальной normal frame

\[
(∂_t-β^i∂_i)u=-αΠ,\quad
(∂_t-β^i∂_i)Π+αγ^{ij}∂_iΦ_j=0,\quad
(∂_t-β^j∂_j)Φ_i+α∂_iΠ=0.
\]

Вдоль unit normal `s^i`, characteristic combinations `w_out=Π+s^iΦ_i`, `w_in=Π−s^iΦ_i` имеют скорости `−β^s+α` и `−β^s−α`, соответственно; transverse `Φ_A`/`u` advect at `−β^s`. При frozen Minkowski это обычная wave pair. Положительный symmetrizer каждого wave block — сумма `Π²+γ^{ij}Φ_iΦ_j` (для системы metric components берётся положительная Euclidean сумма по component labels). При `α>0`, `γ_ij>0`, smooth bounded coefficients он uniform на compact nonsingular domain. Definition constraints `Φ_i−∂_iu=0` дополнительно должны эволюционировать/контролироваться.

**Симметричная/сильная гиперболичность:** full Cartesian first-order GH wave system с prescribed/algebraic `H(g,x)` имеет полную характеристическую базу указанного wave block и положительную symmetrizer; это локальное continuum утверждение away from degeneracy. Скалярный principal block — такой же wave operator. В сферической reduction principal radial speeds остаются `−β^r±α√(γ^{rr})`, но коэффициенты `3/r`, tensor parity и constraint/gauge couplings — lower-order только если variables regularized как Cartesian tensor fields. Это не доказывает uniform estimates для наивных `g_ab,R` variables в `R→0`.

Из этого следует стандартная **локальная nonlinear Cauchy well-posedness** для full quasilinear wave GH+scalar system: smooth initial data в достаточном Sobolev regularity (например, `H^s×H^{s−1}`, `s>4/2+1=3`), Lorentzian nondegenerate `g`, positive spatial metric/lapse, smooth algebraic gauge source and potential, и совместимые Einstein/GH/first-order constraints дают решение на некотором времени до потери regularity/nondegeneracy. Smooth SO(4)-equivariant Cartesian subspace сохраняется уравнениями, так что регулярный центр включён как внутреннее smooth point, а не IBVP boundary. Это локальная Cauchy теорема; она не даёт глобального существования, контроля при curvature blow-up, outer finite-boundary CPBC или конкретного first-order reduced radial discretization.

## 3. GH constraint propagation

Контрактированный Bianchi identity и matter conservation `∇^AT_AB=0` (из scalar EOM) дают при `E_AB=0`:

\[
\nabla^AE_{AB}=-\frac12\left(□_5C_B+R_B{}^AC_A\right)=0,
\quad\Longrightarrow\quad □_5C_B+R_B{}^AC_A=0.
\]

Это точная homogeneous subsidiary wave equation для указанной reduced Einstein operator, если matter EOM и divergence identities выполнены. При `C_A=0` решение GH system удовлетворяет исходным Einstein equations. Если scalar EOM нарушена или boundary пропускает входящие constraint mode, заключение не работает. Constraint damping сознательно не включён; стандартные damping additions меняют subsidiary lower-order matrix, но не исправляют неверный boundary projection и требуют eigenvalue/energy proof для конкретного `κ_d`.

В 5D полный `C_A` имеет 5 компонентов; SO(4) invariant sector оставляет time scalar `C_0(t,r²)` и spatial radial vector `C_I=x_I c(t,r²)`. Их регулярная radial projection `C_r` нечётна и `C_r=O(r)` на центре. Hamiltonian/momentum constraints — normal projections исходных Einstein equations, а GH constraints — coordinate/gauge constraints; они независимы как диагностики, хотя Bianchi связывает их propagation при полной системе. First-order definition constraints — третий набор, не отождествлять их с `C_A`.

## 4. Линеаризованный Minkowski IBVP и внешняя CPBC

Линеаризуем `g_AB=η_AB+h_AB`, `H_A=0`, `φ=φ∞+δφ` в vacuum (`V=V'=Λ_eff=0`). Harmonic constraint в линейном порядке `C_A=∂^B\bar h_AB`; reduced metric equations и scalar perturbation удовлетворяют wave equations, а

\[
□_5C_A=0,
\qquad □_5δφ=0.
\]

В radial SO(4) секторе плоский wave operator на scalar component:

\[
□_5u=-∂_t²u+∂_r²u+\frac3r∂_ru.
\]

Он regular для even `u(r)=u_0+u_2r²+…`; запись `3/r` не является физическим источником/границей. Для tensor component те же radial Laplacian terms с соответствующими Cartesian parity/projection terms.

На фиксированной внешней сфере `r=R_out` frozen Minkowski normal characteristic speeds: one incoming and one outgoing wave field per evolved Cartesian tensor component (15 symmetric `h_AB`) and one scalar. Full metric decomposition: 5 incoming constraint/gauge constraints, 5 gauge data, 5 radiative gravitational data in 5D. Under strict SO(4) symmetry, only four metric functions `(A,B,C,D)` survive, two independent incoming GH constraints `(C_0,C_r)`, two gauge modes; no independent SO(4)-invariant tensor gravitational radiation mode. `δφ` carries one physical incoming scalar mode. Boundary counts refer to noncharacteristic timelike boundary; if shift causes a characteristic speed to change sign, count/project using actual `λ` not labels.

### Frozen energy estimate

For one flat wave component set `Π=−u_t`, `Φ=u_r`; then

\[
∂_t\binom{Π}{Φ}+\begin{pmatrix}0&1[локальный путь исключён]
\]

`w_out=Π+Φ` has speed `+1`, `w_in=Π−Φ` speed `−1`. Energy per unit sphere solid angle

\[
E(t)=\frac12\int_0^{R_{out}}r³(Π²+Φ²)dr,
\quad \frac{dE}{dt}=-R_{out}³ΠΦ\big|_{R_{out}}
\]

when regular center flux vanishes. Since `ΠΦ=(w_out²−w_in²)/4`, homogeneous no-incoming `w_in=0` gives nonpositive outer energy rate. For boundary data `w_in=g`, the energy influx is bounded by a constant times `R_out³g²`. The analogous estimate holds componentwise for linearized Cartesian metric waves, before spherical reduction; in radial variables vector/tensor components add angular-gradient terms to the reduced energy and those must be retained.

The subsidiary constraint wave has its own incoming characteristic. With `Π_C=−∂tC_A`, `Φ_C=∂_rC_A`, set

\[
w^{C}_{in}=Π_C−Φ_C=0\quad (r=R_{out})
\]

for homogeneous frozen CPBC (or prescribe incoming constraint data compatible with desired Einstein constraints). This kills incoming subsidiary violations and yields the analogous energy estimate for `C_A`. It **does not by itself specify metric boundary data**: a boundary operator on incoming `h_AB` must be derived whose induced `w^C_in` has this value while independent gauge/physical data close the metric system. That compatibility map is the unsolved part for the proposed SO(4)-reduced GH metric variables. “Set all metric incoming fields to zero” is not a substitute proof of CPBC and generally reflects physical waves.

For a maximally dissipative frozen wave system, prescribe `w_in=S w_out+g`, with spectral norm `||S||≤1` in the characteristic energy norm; zero data gives non-increasing energy. Constraint projector conditions must be imposed in addition and checked for compatibility/rank with the remaining incoming modes. Sommerfeld/outgoing scalar data are an approximation to reflection control, not a CPBC proof.

### Scalar radiation boundary in 5D

For high-frequency outgoing radial wave in 4 spatial dimensions, leading falloff is `δφ≈r^{-3/2}F(t−r)`. Thus a leading Sommerfeld condition is

\[
(∂_t+∂_r)(r^{3/2}δφ)=0\quad\text{at large }R_{out}.
\]

After `u=r^{3/2}δφ`, the flat radial equation has centrifugal term `3/(4r²)`; 5D spacetime has even spatial dimension and retarded waves have tails. This local leading condition is not exact transparent data. Use a derived nonlocal/exact Dirichlet-to-Neumann map or place boundary far enough and restrict to its causal domain; validate reflection separately from constraints.

## 5. Центр: parity subspace, не incoming boundary

В Cartesian-equivariant variables regularity is smoothness in `x_I`, not arbitrary data imposed at `r=0`:

\[
A,C,D,φ,Π_{AB},Π_φ\ \text{even in }r,
\quad B\ \text{radial vector coefficient even, while }g_{tr}=rB\text{ odd},
\quad ∂_ru_{even}(0)=0.
\]

For vectors `v_I=x_Iv(r²)`, the radial component is odd and `v_r=O(r)`; for tensors, radial and tangential eigenvalues coincide at center with difference `O(r²)`. These conditions enforce cancellation of apparent `1/r`/`1/r²` terms in polar equations. Center is not a timelike physical boundary: a regular radial characteristic arriving there continues by symmetry, with zero net `r³` energy flux. Homogeneous Sommerfeld or zeroing all incoming modes at `r=0` would remove/reflection-distort smooth solutions.

Using variables `R(t,r)` directly with `R(0)=0` can be regular, but coefficient `1/R` and gauge determinant must be rewritten in parity-regular combinations before discretization. `R=r e^{S}` with even smooth `S` is one useful representation; the difference between radial and angular metric components must be `O(r²)`. The full Cartesian GH condition is still enforced; choosing `S` does not remove angular terms.

## 6. Exact тест регулярных начальных данных data and GH initial gauge consistency

The local regular initial family тест регулярных начальных данных has `γ_IJ=ψ²δ_IJ` on `R⁴`, `K_IJ=0`, `Π_φ=0`, `V=Λ=0` and smooth radial scalar; its analytic Hamiltonian/momentum constraints vanish. For the GH Cauchy initialization, choose positive lapse `α`, shift `β^I`, and initial metric components

\[
g_{tt}=-α²+γ_{IJ}β^Iβ^J,\quad g_{tI}=γ_{IJ}β^J,\quad g_{IJ}=γ_{IJ}.
\]

Spatial `Φ_IAB` follow spatial derivatives. `Π_AB=−n^C∂_Cg_AB` and scalar `Π_φ=n^A∂_Aφ` require time derivatives consistent with lapse/shift and gauge. GH constraints `C_A|_Σ=0` determine the gauge-compatible combinations of `∂_t g_AB` (or initial lapse/shift derivatives); they are not ensured merely by solving Hamiltonian/momentum constraints.

For a second-order subsidiary wave, both `C_A|_Σ=0` and `n^B∇_B C_A|_Σ=0` are needed. With exact Einstein constraints, scalar EOM and the chosen reduced evolution, Bianchi identities imply the normal derivative condition in the continuum provided the initial gauge variables satisfy `C_A=0`; this implication must be checked explicitly in the reduced equations and not inferred from `K=0`. First-order definition constraints `Φ−Dg=0` must also hold initially and have compatible normal derivatives. Time symmetry of physical geometry does not force time symmetry of coordinates: lapse/shift/gauge-source derivatives may be nonzero. тест регулярных начальных данных conformal data cannot be inserted as GH arrays without this conversion.

At finite outer radius, initial data and boundary data must satisfy corner compatibility at `t=0,r=R_out` for all incoming metric/scalar/constraint fields and their evolution derivatives. No local source currently provides such a complete SO(4) GH initial-boundary construction.

## 7. SBP-SAT contract and what it proves

For a one-dimensional nodal grid with norm matrix `H=Hᵀ>0`, SBP derivative must satisfy

\[
HD+D^TH=B=e_Ne_N^T-e_0e_0^T.
\]

For frozen symmetric first-order system `u_t+A u_r=0`, `SA=(SA)^T`, discrete energy `E_h=½u^T(H⊗S)u` has only boundary flux terms plus SAT penalties. At outer boundary split `A=A^++A^-`; an upwind/maximally dissipative SAT can be written

\[
SAT_R=H^{-1}e_N(A^-\otimes I)(u_N-g),
\]

with sign convention `A^-` containing negative eigenvalues; the resulting homogeneous incoming data makes the boundary contribution non-positive. Left boundary uses the corresponding positive-speed penalty. For wave block it is equivalent to penalizing `w_in` while leaving outgoing `w_out` free. This is a frozen linear energy estimate, not nonlinear Einstein stability.

At a spherical center, ordinary endpoint SAT at `r=0` is not appropriate. Either evolve smooth Cartesian parity extensions on `[-R,R]` and restrict/identify parity, or derive a weighted radial SBP operator for measure `r³dr` with parity-compatible closure and prove its boundary/source identity. Because `3/r` terms are singular-looking, a plain Cartesian SBP `D` plus naively sampled `3/r` source does not inherit the flat Cartesian energy proof. The chosen center operator and tensor parity closure are not yet designed/proven.

For variable coefficients, SAT must use characteristic projectors of the symmetrized normal matrix and estimate coefficient variation; boundary penalties also have to preserve the **subsidiary Einstein constraints**, not only the first-order wave energy. RK time stability requires timestep inside the spatial operator’s fully discrete stability region and compatible treatment of SAT/damping. SBP identity alone does not prove a nonlinear variable-coefficient GH+SO(4) scheme stable.

## 8. Результат по уровням доказанности

| Уровень | Вердикт | Обоснование |
|---|---|---|
| (a) Linearized full Cartesian GH continuum IVP | **GO, restricted** | Wave principal symbol, complete characteristic basis, positive symmetrizer; matter scalar is a wave block. Requires smooth gauge source and nonsingular Lorentzian background. |
| Full nonlinear Cauchy IVP, without finite outer boundary | **GO locally, conditional on regular data** | Quasilinear wave local existence; Cartesian SO(4) symmetry subspace includes regular center. No global-in-time or high-curvature guarantee. |
| (b) SO(4) regular center | **CONDITIONAL** | Exact radial sector is invariant as smooth Cartesian fields; parity yields finite origin. A reduced polar implementation must prove regularized variables/operator closure. |
| (c) Frozen-coefficient outer CPBC | **CONDITIONAL / linear GO for subsidiary block** | Constraint wave incoming mode and componentwise energy estimate shown; metric incoming boundary map satisfying CPBC+gauge/physical closure for this radial reduced system not derived. |
| (d) Full nonlinear variable-coefficient SO(4) continuum IBVP | **UNSOLVED** | Need uniform symmetrizer/energy estimate with chosen H, regular center and compatible CPBC, corner data; standard 4D harmonic theorems do not prove this exact reduction automatically. |
| (e) Full nonlinear SBP-SAT semidiscrete stability | **UNSOLVED** | Frozen symmetric wave SBP-SAT estimate only; weighted center closure, variable coefficients, nonlinear constraints and time integrator remain. |
| (f) Numerically trustworthy 5D Einstein–scalar implementation | **NO-GO to implement now** | No new solver/code exists; Stage A blocker remains. Old #40–41F/отчёт о тестах коллапса–121 are not validated references. |

### Сравнение с одним fallback: double-null

Double-null 5D spherical equations offer direct characteristic directions and no lapse-collapse singularity, but start from two intersecting null legs (not тест регулярных начальных данных Cauchy data), have residual null gauge, center/caustic/domain limits and require their own CPBC/characteristic corner analysis. The old отчёт о double-null проверке-style code is partial and its own null constraints/ledger were not accurate enough for a reference. Recommendation remains GH for reusable Cauchy CPU core, **conditional on Gate A**; double-null may be an independent restricted cross-check, not automatic replacement.

## 9. Одна следующая ограниченная математическая задача

**Вывести boundary operator и energy estimate для линеаризованного SO(4)-invariant 5D GH system на half-line с Cartesian-parity center и timelike outer boundary:** начиная с exact linearization in Cartesian variables, получить explicit incoming metric characteristic combinations, выразить два incoming GH subsidiary fields through them, выбрать ровно дополняющие gauge/radiative data и доказать uniform boundary stability/rank for `α>|β^r|`. Stop criteria: mismatch in incoming count, nonpositive boundary form, or parity closure not equivalent to smooth Cartesian extension. Не переходить к nonlinear solver, пока эта задача не закрыта.

## 10. Источники и границы проверки

Ниже перечислены библиографические ориентиры, сверенные по доступным публичным страницам/abstract/full HTML в этой задаче. Их результаты относятся к общим harmonic/GH системам, в основном в стандартной 3+1D численной относительности; перенос доказательства на 5D SO(4) center/SAT требует отдельного derivation.

- H.-O. Kreiss, J. Winicour, “Problems which are well-posed in a generalized sense with applications to the Einstein equations,” *Class. Quantum Grav.* **23** (2006), S405–S420, DOI [10.1088/0264-9381/23/16/S07](https://doi.org/10.1088/0264-9381/23/16/S07); public preprint [arXiv:gr-qc/0602051](https://arxiv.org/abs/gr-qc/0602051). Abstract states a second-order harmonic Einstein CPBC well-posedness result.
- H.-O. Kreiss, O. Reula, O. Sarbach, J. Winicour, “Well-posed initial-boundary value problem for the harmonic Einstein equations using energy estimates,” *Class. Quantum Grav.* **24** (2007), 5973–5984, DOI [10.1088/0264-9381/24/23/017](https://doi.org/10.1088/0264-9381/24/23/017); public preprint [arXiv:0707.4188](https://arxiv.org/abs/0707.4188). Used only for general harmonic Einstein energy-estimate context.
- L. Lindblom, M. A. Scheel, L. E. Kidder, R. Owen, O. Rinne, “A new generalized harmonic evolution system,” *Class. Quantum Grav.* **23** (2006), S447–S462, DOI [10.1088/0264-9381/23/16/S09](https://doi.org/10.1088/0264-9381/23/16/S09). Formulation reference, not a proof of this reduced center/IBVP.
- O. Rinne, L. Lindblom, M. A. Scheel, “Testing outer boundary treatments for the Einstein equations,” *Class. Quantum Grav.* **24** (2007), 4053–4078, DOI [10.1088/0264-9381/24/16/006](https://doi.org/10.1088/0264-9381/24/16/006), public [arXiv:0704.0782](https://arxiv.org/abs/0704.0782). Full public HTML/text was inspected for the distinction between CPBC, gauge/radiation data and reflection control; numerical tests are 4D Schwarzschild perturbation tests.
- C. Gundlach, G. Calabrese, I. Hinder, J. M. Martín-García, “Constraint damping in the Z4 formulation and harmonic gauge,” *Class. Quantum Grav.* **22** (2005), 3767–3773, DOI [10.1088/0264-9381/22/17/025](https://doi.org/10.1088/0264-9381/22/17/025), [arXiv:gr-qc/0504114](https://arxiv.org/abs/gr-qc/0504114). Lower-order damping is not inserted here without new subsidiary proof.
- B. Strand, “Summation by Parts for Finite Difference Approximations for d/dx,” *J. Comput. Phys.* **110** (1994), 47–67, DOI [10.1006/jcph.1994.1005](https://doi.org/10.1006/jcph.1994.1005). Basis for the algebraic SBP identity, not proof for spherical nonlinear GH.
- O. Rinne, “Stable radiation-controlling boundary conditions for the generalized harmonic Einstein equations,” *Class. Quantum Grav.* **23** (2006), 6275–6300, DOI [10.1088/0264-9381/23/22/013](https://doi.org/10.1088/0264-9381/23/22/013), [arXiv:gr-qc/0606053](https://arxiv.org/abs/gr-qc/0606053). Frozen/high-frequency boundary stability context; no 5D radial center proof.

Локальные источники: [spec equations](mathematical-equations.md), [numerical methods](numerical-methods.md), [тест регулярных начальных данных regular data](../research/h1-v0.2/h1-v0.2-regular-initial-data-birth-test.md), [критерии коллапса и trapped surfaces collapse criteria](../research/h1-v0.2/h1-v0.2-collapse-trapping-bounce-criteria.md), [отчёт о тестах коллапса collapse benchmark](../experiments/experiment-collapse-benchmark.md), [отчёт о double-null проверке double-null reference](../experiments/experiment-double-null-reference.md), [equation cross-audit](../experiments/experiment-equations-cross-audit.md) (файл существует, но локально не подтверждено, что это именно предшествующее обсуждение), [#41F operator audit](../experiments/experiment-41f-compatible-operators.md). Старый код заморожен; ничего из него здесь не запускалось.

Веб-доступ использован только для публичных библиографических страниц/arXiv текста, без скачивания больших репозиториев/бинарников. Полный математический текст всех перечисленных работ не заявляется прочитанным; полный текст Rinne et al. 2007 просмотрен. Ни одна ссылка сама по себе не доказывает nonlinear 5D SO(4)-reduced IBVP.

**Ограниченное уточнение:** frozen Minkowski principal boundary map для reduced SO(4) сектора, её точный rank и локальная энергетическая оценка позднее выведены в [отдельной CPBC заметке](gh-linearized-cpbc-characteristics.md). Нелинейная переменно-коэффициентная задача и regular-center SBP closure по-прежнему не доказаны.
