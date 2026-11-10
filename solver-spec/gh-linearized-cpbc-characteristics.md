# Линеаризованные GH CPBC в 5D SO(4): характеристики, rank и energy gate

## Результат в одном абзаце

Для плоского 5D Einstein–scalar background и SO(4)-инвариантных возмущений удаётся получить явный **ограниченный frozen-coefficient CPBC** на внешней timelike границе. У reduced metric имеется четыре входящих wave combinations: два линейных условия зануляют входящие GH constraint modes, ещё два задают gauge. Rank равен четырём; в spherical sector нет независимых входящих гравитационных волн, но есть одна физическая incoming scalar mode. При homogeneous incoming gauge/scalar data линейный wave energy через outer boundary не возрастает. Это не доказательство nonlinear CPBC и не SBP-SAT stability с регулярным центром.

**Решение этапа:** CONTINUUM linear boundary rank/energy gates пройдены в локальном frozen Minkowski секторе; sandboxed linear toy test теперь математически осмыслен как следующий отдельно разрешаемый implementation scope. Полная nonlinear SO(4) IBVP остаётся `UNSOLVED`; solver code этим отчётом не разрешён.

## 0. Фиксированное действие и GH principal system

Берётся только базовая модель из [physics specification](physics-specification.md):

\[
S=\frac1{2\kappa_5^2}\int d^5x\sqrt{-g}(R-2\Lambda_5)
-\int d^5x\sqrt{-g}\left[\frac12(\nabla\phi)^2+V(\phi)\right],
\quad \kappa_5^2=8\pi G_5,
\]

с сигнатурой `(-++++)`, `c=\hbar=1`, `[\kappa_5^2]=L^3`, `[\phi]=L^{-3/2}`, `[V]=L^{-5}`. Если на физической границе требуется вариационная постановка, действие дополняется согласованным GHY/асимптотическим членом; численная worldtube не является физической стенкой. Вариация даёт `G_AB+Λ₅g_AB=κ₅²T_AB`, `□φ−V′=0`; в trace-reversed форме
`R_AB=κ₅²(∂_Aφ∂_Bφ+2Vg_AB/3)+2Λ₅g_AB/3`. GH constraint определён относительно reference connection в §1. Используется reduced equation
`R_AB−∇_(A C_B)=κ₅²(∂_Aφ∂_Bφ+2Vg_AB/3)+2Λ₅g_AB/3`, без constraint damping. При `Λ₅=0`, `V(φ∞)=V′(φ∞)=V′′(φ∞)=0`, `φ=φ∞` и `g=η` линеаризация даёт `□₅ h_AB=0`, `□₅ϕ=0`; trace reversal — постоянное обратимое линейное преобразование, поэтому wave principal symbol тот же. Это обоснование относится к principal part вне координатного центра, не к complete nonlinear reduced operator.

## 1. Reference-metric GH и точная 5D редукция gauge constraints

Сигнатура `(-++++)`, reference metric — 5D Minkowski

\[
d\hat s^2=-dt^2+dr^2+r^2d\Omega_3^2,
\]

и gauge-конвенция

\[
C^A:=g^{BC}\left(\Gamma^A{}_{BC}-\widehat\Gamma^A{}_{BC}\right)-H^A=0.
\]

Разность connections — tensor; `H^A` выбирается smooth SO(4)-эквивариантным reference-vector source без производных evolving metric в principal-part analysis. При `g=\hat g`, `H=0`: `C^A=0`. В Cartesian coordinates `\hat\Gamma=0`, так что это совпадает с usual GH `Γ^A−H^A`.

Для `ds²=h_abdx^adx^b+R²dΩ₃²` и `\hat h_abdx^adx^b=-dt²+dr²`, `\hat R=r` base-компоненты равны

\[
C^a=h^{bc}\left(\Gamma^a{}_{bc}[h]-\widehat\Gamma^a{}_{bc}[\hat h]\right)
-\frac3R\nabla^aR+\frac{3r}{R^2}\widehat\nabla^a r-H^a.
\]

Здесь `Γhat[h]=0` в этих base coordinates, а `\hat\nabla^a r=\hat h^{ab}∂_br`. При `R=r`, `h=\hat h` два warp terms сокращаются. Для angular indices `i` сферическая connection одинакова у `g` и reference, поэтому `C^i=-H^i`; выбрать `H^i=0` в reference-metric formulation. Это cancellation теряется, если написать только `Γ[h]^a−H^a` или использовать ordinary GH в spherical chart без angular connection.

Эквивалентная регулярная SO(4)-эквивариантная Cartesian ansatz:

\[
g_{00}=A(t,r^2),\quad g_{0I}=B(t,r^2)x_I,
\quad g_{IJ}=C(t,r^2)\delta_{IJ}+D(t,r^2)x_Ix_J,
\quad φ=φ(t,r^2),\quad I=1,…,4,
\]

\[
ds²=A dt²+2rB\,dt\,dr+(C+r²D)dr²+r²C dΩ₃²,
\quad R=r\sqrt C.
\]

Все пять Cartesian GH constraints сохранены; SO(4) оставляет `C_0` и радиальный vector component `C_r` как две независимые функции. Angular equations не отбрасываются: они представлены tangential metric component `C`/`R` и full Cartesian tensor equation.

## 2. Linearized field equations and mode count

Background: `g_AB=η_AB`, `R=r`, `φ=φ∞` at a vacuum with `V(φ∞)=V'(φ∞)=0`, `Λ_eff=0`; local outer patch `r=R_out>0`. The calculations use the massless linear scalar limit `V''(φ∞)=0`; if `V''≠0`, the scalar principal characteristics stay null but the energy and finite-frequency boundary residual acquire a lower-order mass term. Write `g=η+h`, `φ=φ∞+ϕ`, trace reversal `\bar h_AB=h_AB−η_AB h/2`. For `C_A=∂^B\bar h_AB`:

\[
□_5\bar h_{AB}=0,\qquad □_5ϕ=0,
\qquad □_5 C_A=0.
\]

The scalar stress is quadratic about a constant vacuum, so metric/scalar principal blocks decouple at linear order. `C_A` is the full five-component GH constraint; its subsidiary equation follows from the linearized Bianchi identity.

At one boundary point use a locally orthonormal frame `(t,x,y^1,y^2,y^3)` with `x` outward normal; freeze `r^{-1}` lower-order and angular-basis terms for the principal boundary symbol. SO(4)-invariant metric perturbations have four independent components

\[
h_{00}=a,\quad h_{0x}=b,\quad h_{xx}=c,
\quad h_{y^py^q}=d\,\delta_{pq}\quad(p,q=1,2,3),
\]

and trace `h=−a+c+3d`. The independent radial GH constraints are

\[
C_0=-\partial_t\frac{a+c+3d}{2}+\partial_x b,
\qquad
C_x=-\partial_t b+\partial_x\frac{a+c-3d}{2}.
\]

The other three spatial components are tangential vectors and vanish in the invariant scalar sector. The full 5D graviton has five polarizations, but none is an SO(4)-invariant spherical radiative mode; by spherical symmetry/Birkhoff-type reduction the metric modes here are constrained/gauge response to scalar matter, not five independent incoming tensor waves.

## 3. Characteristics with lapse/shift

For each wave variable `u` (metric component or scalar), choose `Π=−n^A∂_Au`, `Φ_s=s^i∂_iu`, where `n` is slice normal and `s` outward unit spatial normal. The frozen normal principal block has

\[
λ_{out}=-β^s+α,\quad w_{out}=Π+Φ_s;
\qquad
λ_{in}=-β^s-α,\quad w_{in}=Π-Φ_s,
\]

and advected/zero-speed (`−β^s`) metric/derivative fields. In coordinates with radial metric `γ_rr`, null coordinate speeds are `dr/dt=−β^r±α\sqrt{γ^{rr}}`. A fixed outer worldtube is noncharacteristic/timelike if `α>|β^s|`; then one wave family enters and one leaves. At the Minkowski background `α=1, β=0`, speeds are `+1,-1,0`. If a general background shifts either family through zero, incoming/outgoing count is recalculated from the actual signs.

Full 5D Cartesian GH contains 15 symmetric metric wave components, hence 15 incoming metric characteristics before symmetry; the full constraint vector has five components. SO(4) scalar projection leaves four metric incoming amplitudes `(q_a,q_b,q_c,q_d)` and two GH constraints. The scalar has one additional incoming physical amplitude `q_φ`; no symmetry-invariant gravitational radiation datum remains. First-order reduction definition constraints `Φ_i−∂_iu=0` are separate and need extra boundary treatment in a future first-order implementation; this derivation is for the second-order wave system/its exact first-order characteristic representation.

## 4. Explicit incoming constraint projection and exact rank

Set

\[
q_a=w^{in}_{00},\quad q_b=w^{in}_{0x},\quad
q_c=w^{in}_{xx},\quad q_d=w^{in}_{yy}=w^{in}_{zz}=w^{in}_{ww}.
\]

For an incoming plane profile `F(t+x)`, `w_in=Π−Φ=−2∂_{t+x}F`. Substitution in `C_A` gives the incoming parts of the constraints:

\[
C^{in}_0=\frac14(q_a+q_c+3q_d-2q_b),
\qquad
C^{in}_x=\frac14(-q_a-q_c+3q_d+2q_b).
\]

The homogeneous no-incoming-constraint boundary operator is

\[
L_0(q):=q_a+q_c+3q_d-2q_b=0,
\qquad
L_x(q):=-q_a-q_c+3q_d+2q_b=0.
\]

Its constraint matrix on four amplitudes is

\[
\begin{pmatrix}1&-2&1&3[локальный путь исключён]
\]

which has **rank 2**; its nullspace is exactly

\[
q_d=0,\qquad q_b=\frac{q_a+q_c}{2}.
\]

For the subsidiary wave `C_A`, use the same convention `z_A^-:=−(∂_t+∂_x)C_A`. On a pure incoming profile `F(t+x)`, `C_A^{in}=L_A(q)/4` and hence `z_A^- =−2∂_t C_A^{in}=−½∂_tL_A(q)` in the frozen constant frame. Requiring `L_A(q)=0` for every boundary time therefore kills the incoming subsidiary characteristic (with the corner compatibility stated below); it does not set the outgoing constraint mode to zero.

Thus there are exactly two free incoming metric functions. They are gauge, not radiative: a residual harmonic gauge vector `ξ_A=F_A(t+x)` produces `h_AB=2∂_(A ξ_B)` with `d=0` and `b=(a+c)/2`. A continuum boundary prescription closing all four incoming metric modes is therefore

\[
q_a=g_0(t),\quad q_c=g_r(t),\quad
q_d=0,\quad q_b=\frac{g_0(t)+g_r(t)}2.
\]

On the ordered vector `q=(q_a,q_b,q_c,q_d)`, the completed boundary operator has coefficient matrix
\[
M_{\rm in}=\begin{pmatrix}1&0&0&0[локальный путь исключён]
\qquad \operatorname{rank}M_{\rm in}=4,
\]
with right-hand side `(g_0,g_r,0,0)`. Thus the two constraint rows plus two gauge-fixing rows specify exactly the four incoming metric fields. `g_0,g_r` are the two residual incoming gauge data and must be compatible with the chosen gauge source and initial corner. For homogeneous gauge data, all four incoming metric amplitudes vanish. For nonzero gauge data they prescribe coordinate forcing; they do not represent gravitational radiation. This is a **principal frozen-coefficient** result; finite-radius `O(1/R_out)` terms and nonzero background curvature modify lower-order boundary equations and require a separate estimate.

The scalar has independent incoming physical data `q_φ=w^{in}_φ`; choose a prescribed incoming scalar pulse or homogeneous no-incoming `q_φ=0`. It is not a gravitational constraint. No extra boundary data should be imposed on the outgoing metric/scalar fields.

Since `C^{in}_A` is set to zero, its incoming subsidiary wave characteristic (a tangential/time derivative of `C^{in}_A`) also vanishes. With initially compatible corner data, this is a valid linear homogeneous CPBC for the frozen subsidiary system. Prescribing only scalar-style Sommerfeld on all four metric components would not enforce `L_0=L_x=0` and generally injects constraint modes.

Для начального угла `t=0,r=R_out` входящие данные не свободны от совместимости: initial GH data должны удовлетворять `C_A=0`, а их normal/time derivatives — соответствующему subsidiary constraint на начальной поверхности; значения `g_0,g_r` и их совместимые временные производные должны совпадать с начальными входящими характеристиками. Для тест регулярных начальных данных геометрические Cauchy data сначала задают spatial metric и `K_IJ`; GH-представление дополнительно требует выбора lapse/shift, вычисления `∂_t g_AB` по этим данным и согласования gauge source и его начального изменения так, чтобы `C_A=0` и `∂_t C_A=0`. Сам факт, что тест регулярных начальных данных решает Hamiltonian/momentum constraints, ещё не задаёт эти GH gauge data. Их явная конструкция не выведена здесь.

## 5. Frozen energy flux

Use the full Cartesian positive wave energy for linearized metric perturbations and scalar,

\[
E(t)=\frac12\int_{r<R_{out}}\left[
\sum_{A,B}(Π_{AB}^2+\partial_Ih_{AB}\partial_Ih_{AB})
+Π_φ^2+\partial_Iϕ\partial_Iϕ\right]d^4x.
\]

For the SO(4) boundary components, the incoming/outgoing metric norm has Frobenius weights

\[
\|q\|^2=q_a^2+2q_b^2+q_c^2+3q_d^2;
\]

the factor 2 is the off-diagonal pair and factor 3 is the three equal angular components. At a fixed outer sphere, the principal energy flux is

\[
\frac{dE}{dt}\Big|_{R_{out}}
=-\frac{\Omega_3R_{out}^3}{4}
\left(\|w_{out}\|^2-\|w_{in}\|^2+|w^{φ}_{out}|^2-|w^{φ}_{in}|^2\right),
\quad \Omega_3=2π^2.
\]

With the CPBC above and `g_0=g_r=q_φ=0`, `w_in=0`, hence the boundary contribution is non-positive. Nonzero gauge/scalar data contribute a bounded energy input proportional to `R_out³(\|g\|²+|q_φ|²)`. In the constraint subsystem, no incoming `C_A` means zero incoming subsidiary flux; outgoing constraint violations may leave the domain. The two conditions do not force the outgoing constraints to zero and do not overdetermine the four incoming metric modes.

This is an energy estimate for a **constant-coefficient local planar boundary symbol**, with `O(1/R)` and background-gradient terms frozen out. For fixed but finite `R_out`, those lower-order terms can be bounded by a Grönwall volume term if coefficients are smooth and `R_out>0`; that observation is not yet a uniform estimate for the exact variable-coefficient reduced system. The speed formula with shift is given in §3, but the rank-2 constraint projection above is explicitly derived only for the zero-shift Minkowski frame; for nonzero shift the incoming constraint map must be recomputed in the actual boundary-adapted frame. A global-in-time bound also needs uniformly bounded lapse, shift, spatial metric and symmetrizer.

## 6. Center regularity — separate from outer CPBC

Cartesian smoothness imposes even Taylor series

\[
A=A_0+A_2r^2+O(r^4),\quad C=C_0+C_2r^2+O(r^4),
\quad D=D_0+D_2r^2+O(r^4),
\]

\[
B=B_0+B_2r^2+O(r^4),\quad
φ=φ_0+φ_2r^2+O(r^4),\quad
H^0=H^0_0+O(r^2),\quad H^I=x^IH_1(t,r^2).
\]

Thus `g_tr=rB=O(r)` is odd; `g_rr=C+r²D` and angular coefficient `R²=r²C` agree at leading order; `R=r√C=√C_0 r+O(r³)`. For scalar first-order variables `Π_φ` is even, `Φ_r=∂rφ=2φ_2r+O(r³)` odd. Radial metric/angular eigenvalues coincide at the center; radial vector GH constraint `C_r=O(r)` is odd and vanishes there, while Cartesian `C_I=x_I c(r²)` is smooth.

The apparent `3/r` radial Laplacian terms cancel on smooth even fields because `∂r u=2u_2r+O(r³)`, giving `u_rr+3u_r/r=8u_2+O(r²)`. Similar cancellation for tensor components must be derived from the full Cartesian ansatz, not guessed component-by-component. The center is not an exterior boundary: incoming and outgoing radial characteristics continue through the origin under parity; the `r³` flux vanishes for bounded regular fields. Imposing homogeneous incoming Sommerfeld data at `r=0` is wrong.

The continuum parity subspace is smooth and invariant. A discrete weighted radial SBP closure preserving the five-dimensional energy and the coupled tensor parity has **not** been constructed here. A standard SBP endpoint penalty at `r=0` is not justified. For a discrete proof one must either (a) derive parity-compatible weighted SBP operators for `r³dr` and all tensor components including the angular source terms, or (b) use a Cartesian multi-D reference discretization; the latter is not the proposed 1D CPU core. The outer CPBC proof does not solve this center discretization problem.

## 7. Prospective SBP-SAT discrete contract (not a full radial proof)

На outer annulus (`r_min>0`) возьмём pointwise state `U=(Π,Φ_s)` для каждой frozen metric wave component и `ϕ`, с `Π=−∂_t u`, `Φ_s=∂_s u`. В локальном orthonormal frame normal principal system `U_t+A U_s=0`, `A=[[0,1],[1,0]]`, symmetrizer `I`; characteristic projectors `P_±=(I±A)/2`, `w_out=P_+U`, `w_in=P_-U` (с нормировкой `w=Π±Φ` из §3). Для SO(4) metric amplitudes использовать Frobenius inner product `W_c=diag(1,2,1,3)` на `(00,0s,ss,tt)`; scalar weight — 1. Поскольку `W_c` коммутирует с `A`, та же симметризация действует на компонентный блок.

На узлах `s_j`, `j=0,…,N`, задать diagonal-norm SBP pair `(D,H)`, `H>0`,
`HD+DᵀH=B=e_Ne_Nᵀ−e_0e_0ᵀ`. Для каждого frozen constant normal block на внешнем конце полудомена semi-discrete operator имеет вид

\[
\dot U+(D\otimes A)U
=H^{-1}e_N\otimes A^-\,(U_N-g),
\qquad A^-:=P_- A,
\]

где формула записана на один wave block; для stacked четырёхкомпонентного metric поля операторы `A,A^-` заменяются на `I_4⊗A,I_4⊗A^-`, а `g` задаётся только во входящем подпространстве. Для metric block входящие data параметризуются CPBC из §4: два ограничения `L_0=L_s=0`, два совместимых gauge data; scalar получает одно отдельно заданное `w_in^φ`. При homogeneous data энергия `E_h=½ Uᵀ(H⊗W_c⊗I_2)U` удовлетворяет
`dE_h/dt=−½ U_Nᵀ(W_c⊗A)U_N+U_Nᵀ(W_c⊗A^-)(U_N−g)` на внешнем конце; spectral splitting даёт `dE_h/dt≤−½||w_out||²_{W_c}+½||g_in||²_{W_c}`. Это локальная frozen wave/SBP estimate; boundary data are not applied to outgoing fields. Знак SAT выбран для equation `U_t+A DU=SAT` и domain ending at `s_N`.

Для radial scalar wave на annulus `r_min>0` есть более точный weighted conservative pair,
`Π_t=−r^{-3}D(r³Φ)`, `Φ_t=−DΠ`, с norm `E_h=½(ΠᵀH R³Π+ΦᵀH R³Φ)`. Обычное SBP identity даёт точное cancellation interior и только boundary flux `−[r³ΠΦ]`; это проверка scalar radial geometry, не доказательство tensor GH radial operator. Для reduced metric нужно вывести такие же compatible angular/geometric source terms and constraint variables; нельзя просто применить четыре независимых radial scalar stencils.

При переводе локального boundary estimate на полную сферу его surface term умножается на `Ω₃R_out³`; radial weighted norm выше учитывает этот множитель локально. Эта площадь не меняет incoming rank, но её нельзя опустить в физическом energy ledger.

На всём `0≤r≤R_out` требуются smooth parity из §6: `A,C,D,B,φ` even coefficient functions of `r²`, `g_tr=rB` odd, `R=r√C`, radial/tangential metric eigenvalues согласованы. Дискретное closure должно либо (i) построить weighted SBP pair для `r³dr` с отдельным аналитическим center limit и доказать cancellation tensor angular sources, либо (ii) развернуть smooth Cartesian parity extension и определить совместимую проекцию всего tensor system. Ни одна конструкция здесь не построена; не ставить penalty как на физическую входящую границу в `r=0`. Semi-discrete test обязателен: энергия/constraint boundary form и parity residual при refinement. Для explicit time step выбрать `Δt` по spectrum конкретного semi-discrete SAT operator внутри stability region интегратора; для RK4 это проверить numerically/аналитически на реальных eigenvalues (CFL не задаётся произвольной константой). Отдельно мониторить discrete `C_A`, definition constraints первого порядка и boundary data/corner residual. SBP identity для одного `D` не доказывает эти coupled свойства.

Итак, линейная frozen SAT формула выше implementable на regular annulus; продление на radial coefficients, origin, first-order reduction constraints и nonlinear GH не доказано.

## 8. Scalar transparency is not CPBC

For four spatial dimensions the asymptotic outgoing scalar behaves at high frequency as

\[
ϕ\sim r^{-3/2}F(t-r),
\qquad (∂_t+∂_r)(r^{3/2}ϕ)=0
\]

to leading Sommerfeld order. Substitution `u=r^{3/2}ϕ` yields centrifugal potential `3/(4r²)`; even spatial dimension means tails inside the light cone. Thus local Sommerfeld minimizes leading reflection only at a sufficiently remote boundary; it is neither exact nonreflecting data nor a GH CPBC. The scalar boundary flux and reflection error are measured separately.

## 9. What is and is not proved

| Claim | Status | Exact scope |
|---|---|---|
| Reference-metric full 5D GH constraint includes warp/angular terms | **PROVEN_RESTRICTED** | Formula in §1 for the stated flat spherical reference metric; Cartesian implementation is equivalent. |
| Linearized SO(4) incoming metric mode count is four; two constraints + two gauge modes | **PROVEN_RESTRICTED** | Frozen Minkowski normal symbol; explicit rank-2 matrix and 2D nullspace in §4. No spherical radiative gravitational mode. |
| Frozen linear outer CPBC rank closes all incoming metric modes and scalar data separately | **PROVEN_RESTRICTED** | Boundary matrix rank 4 plus one scalar incoming datum; corner compatibility assumed. |
| Frozen linear energy estimate with homogeneous CPBC | **PROVEN_RESTRICTED** | Constant coefficients/local planar boundary; positive Frobenius energy; no incoming metric/scalar data. |
| Semi-discrete SBP-SAT estimate for frozen Cartesian wave blocks on an annulus | **PROVEN_RESTRICTED** | Standard diagonal-norm SBP identity plus characteristic SAT gives the stated energy inequality; radial scalar conservative form also cancels exactly in the weighted volume sum. |
| Discrete metric GH operator with spherical source terms and regular center | **UNSOLVED** | No compatible tensor SBP closure, discrete Bianchi/constraint estimate, or center parity proof has been constructed. |
| Full nonlinear variable-coefficient 5D SO(4) CPBC | **UNSOLVED** | Need exact boundary operator including `R`, shift/lapse, gauge source and lower-order terms, uniform boundary stability and corner compatibility. |
| Finite-radius outgoing scalar exact transparency | **COUNTEREXAMPLE to local Sommerfeld exactness** | Tails and `3/(4r²)` potential; leading Sommerfeld is approximate. |
| SBP-SAT regular-center/full coupled numerical stability | **UNSOLVED** | Weighted radial tensor parity closure and semidiscrete energy identity absent. |
| Nonlinear 5D Einstein-scalar implementation readiness | **CONDITIONAL / BLOCKED** | Linear boundary toy is derivable; full center/operator and nonlinear CPBC proofs are missing. |

## 10. Decision gate and narrow fallback

### GH recommendation

Continue GH for the reusable Cauchy core, but only the **linearized outer-boundary toy** is ready for a future separately authorized implementation issue: Cartesian/frozen Minkowski, outer annulus `r_min>0` (so it does not pretend to solve the center), four metric modes with the explicit rank-2 CPBC plus two gauge data, one scalar incoming mode, analytic incoming/outgoing pulses and constraint injection. It must verify the predicted nonpositive energy flux and convergence of the boundary residual; this is not a production solver or approval to code now.

**Одна следующая математическая задача:** вывести variable-coefficient continuation boundary map и uniform energy estimate на SO(4)-invariant annulus `r_min≤r≤R_out`, `r_min>0`, при гладком фоновом lapse/shift и timelike outer worldtube. В задачу не включать центр. Pass gate: incoming constraint projector сохраняет rank 2, ровно два свободных gauge data дополняют его до всех четырёх incoming metric fields, а symmetrized boundary form остаётся неположительной для homogeneous data с контролируемыми lower-order terms. Иначе GH outer-boundary formulation требует пересмотра до nonlinear solver. Центровое weighted-SBP/parity доказательство остаётся отдельным более поздним implementation gate, не частью этого следующего теста.

### One fallback: double-null

Double-null evolves characteristic legs and makes radial null directions explicit; it avoids a lapse-collapse singularity in a local patch but does **not** eliminate outer boundary data if the domain is finite. It needs compatible characteristic constraints at intersecting initial legs and outer null/timelike boundary, center regularity, residual null gauge control and caustic termination rules. It does not automatically handle тест регулярных начальных данных Cauchy data or furnish global completeness. The archived отчёт о double-null проверке-style implementation had large coupled null-constraint/mass-ledger defects and is not a validated substitute. Thus GH remains preferred for Cauchy data if its nonlinear CPBC/center gates pass; double-null is a separately derived cross-check, not an automatic pivot.

## 11. Next validation gate (not performed)

For a later authorized linear toy: (1) analytic Minkowski with exact zero constraints; (2) outgoing and incoming metric pure-gauge pulses satisfying/violating the two `L` conditions; (3) boundary constraint injection showing incoming `C_A` is suppressed only by CPBC; (4) independent scalar pulse and reflected-energy comparison; (5) refinement of boundary residual/energy flux and a large-domain reference; (6) verify center parity only in a separate analytic test, not claim a discretized center proof. Error acceptance must follow formal order/Richardson estimates, not an arbitrary fixed epsilon. No tests/simulations were run for this report.

## 12. Литература и локальные источники

Локально использованы [GH-IBVP аудит GH-IBVP audit](gh-ibvp-cpbc-proof.md), [equation specification](mathematical-equations.md), [numerical methods](numerical-methods.md), [тест регулярных начальных данных](../research/h1-v0.2/h1-v0.2-regular-initial-data-birth-test.md), [критерии коллапса и trapped surfaces](../research/h1-v0.2/h1-v0.2-collapse-trapping-bounce-criteria.md), [отчёт о тестах коллапса](../experiments/experiment-collapse-benchmark.md), [отчёт о double-null проверке](../experiments/experiment-double-null-reference.md), [equation cross-audit](../experiments/experiment-equations-cross-audit.md) (его соответствие предшествующее обсуждение не установлено) и [#41F](../experiments/experiment-41f-compatible-operators.md). Все перечисленные локальные файлы существуют.

GH/CPBC ориентиры уже приведены в GH-IBVP аудит: Kreiss–Winicour (2006), Kreiss–Reula–Sarbach–Winicour (2007), Lindblom et al. (2006), Rinne et al. (2007), Gundlach et al. (2005), Strand (1994). Их область — общая harmonic/GH теория, преимущественно стандартные 3+1D Einstein IBVP и SBP алгебра; они не объявляются прямым доказательством nonlinear 5D reduced center/CPBC системы. Полный текст Rinne et al. 2007 был доступен/просмотрен в предыдущем GH-аудите; здесь библиографические details не перепроверялись по внешней сети.

**Не выполнялись:** исходный код, numerical runs, tests, solver design beyond this math gate, discretization experiment, package install или изменение старого архива. `T_P`, `t_D` и координатное время не отождествляются; H2 вне scope.
