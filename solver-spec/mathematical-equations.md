# Математическая спецификация минимальной модели

## 1. Пространство-время, действие и размерности

Используются `D=5`, `η_AB=(-++++)`, `c=ℏ=1`, `κ₅²=8πG₅`. Объёмное действие и уравнения заданы в [physics-specification.md](physics-specification.md). `x^A`, `R`, `a` имеют `L`; `g_AB` безразмерна; `[κ₅²]=L³`, `[Λ₅]=L⁻²`, `[φ]=L⁻³ᐟ²`, `[V]=[ρ]=[p]=L⁻⁵`. `K_AB=L⁻¹`, Ricci scalar `L⁻²`.

## 2. Точная SO(4)-сферическая редукция

Если вся задача сохраняет SO(4), метрика представляется как warped product

\[
ds_5^2=h_{ab}(x)dx^a dx^b+R(x)^2d\Omega_3^2,
\qquad a,b=0,1,
\]

где `h_ab` — Lorentzian metric орбитального 1+1 пространства, `R>0` — ареальный радиус, `dΩ₃²` — unit round sphere. Это не координатный выбор `R=r`; R является динамическим полем. Для scalar `φ(x^a)` точные warped-product Ricci тензоры:

\[
R^{(5)}_{ab}=R[h]_{ab}-\frac3R\nabla_a\nabla_bR,
\quad
R^{(5)}_{ij}=\left[2(1-X)-R\Box_hR\right]\gamma_{ij},
\quad X=h^{ab}\partial_aR\partial_bR,
\]

а

\[
\Box_5\phi=\Box_h\phi+\frac3R\nabla^aR\nabla_a\phi.
\]

Из trace-reversed Einstein–scalar equations, `R_AB=κ₅²[∂_Aφ∂_Bφ+2Vg_AB/3]+2Λ₅g_AB/3`, следует reduced covariant system

\[
R[h]_{ab}-\frac3R\nabla_a\nabla_bR
=\kappa_5^2\left(\nabla_a\phi\nabla_b\phi+\frac23Vh_{ab}\right)
+\frac23\Lambda_5h_{ab},
\]

\[
2(1-X)-R\Box_hR=\frac23(\kappa_5^2V+\Lambda_5)R^2,
\qquad
\Box_h\phi+\frac3R\nabla R\cdot\nabla\phi-V'(\phi)=0.
\]

Эти уравнения — точная symmetry reduction непрерывных field equations, а не готовая гиперболическая first-order схема. Bianchi identity делает часть компонент зависимой только если scalar equation и gauge constraints выполнены; конкретный независимый эволюционный набор и propagation system должны быть выведены/проверены до кодирования.

## 3. Gauge, эволюционные поля и well-posedness gate

Рекомендуемая кандидатная gauge для horizon-capable ветки — generalized harmonic (GH) на 2D orbit space с координатами `(t,x)`:

\[
C^a:=h^{bc}\Gamma^a{}_{bc}[h]-H^a(t,x,h,R,\phi)=0.
\]

`H^a` задаётся гладким gauge source; он не является физической силой. GH-reduced Ricci заменяет principal second derivatives на quasilinear wave operator и добавляет `\nabla_{(a}C_{b)}`; полное 5D reduced system с warp factor нужно преобразовать без пропуска `3/R` членов. Основные эволюционные variables: три независимых компоненты `h_ab`, `R`, `φ` и их first-order derivative variables (временные/пространственные производные); lapse/shift выводятся из выбранного 2D ADM-разложения или из `h_ab`. Constraint variables: GH `C^a`, definitions of derivative variables, Einstein constraints и scalar equation residual.

**Блокер до реализации:** локальная 1+1 GH reduction должна иметь явно выведенный principal symbol и доказанную strong hyperbolicity на разрешённом `R>0`, `det(h)<0` domain при выбранном `H^a`; boundary conditions должны быть constraint-preserving. Формулы выше не являются доказательством этого свойства. Если hyperbolicity/constraint boundary analysis не закрыта, Gate A останавливает проект; нельзя объявить систему well-posed только из well-known 5D Einstein GH результата без проверки warp-factor reduction.

Polar-areal gauge `ds²=−α²dt²+f⁻¹dr²+r²dΩ₃²` полезна как контрольный restricted backend вне `f=0`, но останавливается/вырождается у marginal sphere и не подходит как единственная сильнополевая эволюция. Double-null `ds²=−Ω²du dv+R²dΩ₃²` избегает lapse-радиальной singularity локально, но имеет характерные начальные данные, affine/residual gauge свободу, center/caustic issues и ограниченный patch; его уравнения не смешивать по массивам координат с Cauchy GH.

## 4. Cauchy data, constraints и initial-data contracts

На spacelike 4D slice `Σ`, с положительной spatial metric `γ_IJ`, future normal `n^A`, convention `K_IJ=−γ_I{}^Aγ_J{}^B∇_A n_B`, определить `Π=n^A∇_Aφ` и `D_Iφ`. При `j_I=−γ_I{}^A T_ABn^B`:

\[
E=\frac12\left(\Pi^2+D_I\phi D^I\phi\right)+V,
\quad j_I=-\Pi D_I\phi,
\]

\[
{}^{(4)}R+K^2-K_{IJ}K^{IJ}=2\kappa_5^2E+2\Lambda_5,
\quad D_J(K^J{}_I-\delta^J_IK)=\kappa_5^2j_I.
\]

Для SO(4)-симметричной slice координаты `(r,Ω₃)`:

\[
\gamma_{IJ}dx^Idx^J=L(r)^2dr^2+R(r)^2d\Omega_3^2,
\quad K^I{}_J=\mathrm{diag}(K_r,K_\Omega,K_\Omega,K_\Omega),
\quad \phi=\phi(r),\ \Pi=\Pi(r).
\]

Положить `K=K_r+3K_Ω`. Тогда

\[
{}^{(4)}R=\frac{6}{R^2}\left[1-\left(\frac{R'}L\right)^2\right]
-\frac{6R''}{L^2R}+\frac{6L'R'}{L^3R},
\]

\[
{}^{(4)}R+6K_\Omega(K_r+K_\Omega)=2\kappa_5^2E+2\Lambda_5,
\]

\[
-3K_\Omega'+3\frac{R'}R(K_r-K_\Omega)=\kappa_5^2j_r,
\quad E=\tfrac12[\Pi^2+(\phi'/L)^2]+V,
\quad j_r=-\Pi\phi'.
\]

Эти reduced constraints — два независимых radial equations. Initial-data solver обязан зафиксировать radial coordinate/conformal ansatz и свободные функции; нельзя независимо задать все `L,R,K_r,K_Ω,φ,Π`. Набор тест регулярных начальных данных `γ=ψ²(dr²+r²dΩ₃²)`, `K=Π=0`, `V=Λ=0` является одним constraint-consistent example для massless scalar, не ECSK/false-vacuum data.

Точный тест регулярных начальных данных control family следует сохранять как fixture: `q(r)=Q(r²/L₀²)e^{-r²/L₀²}`, `[Q]=L⁻²`; `I(r)=∫₀ʳs³q(s)ds`; `Δ_{R⁴}ψ=−q`, `ψ(∞)=1`, `ψ'(r)=−I(r)/r³`, `ψ(r)=1+∫ᵣ^∞I(s)s^{-3}ds`; `γ=ψ²(dr²+r²dΩ₃²)` и `φ'(r)=√[6q(r)/(κ₅²ψ(r))]`. Тогда `{}^{(4)}R=6qψ^{-3}=2κ₅²E` точно в continuum. Для этой конкретной regular time-symmetric family `Q L₀²<1` — достаточный bound, чтобы `dR/dℓ>0` и начальные круглые S³ не были trapped/anti-trapped. Это initial slice result, не future collapse test; нельзя переносить bound на иные potentials, ECSK matter или иной radial ansatz.

### Origin parity и регулярность

В регулярном Cartesian центре выбрать `R(ℓ)=ℓ+O(ℓ³)` для proper radial distance `ℓ`; `φ=φ₀+φ₂ℓ²+O(ℓ⁴)`, `Π=Π₀+Π₂ℓ²+…`; radial/tangential extrinsic-curvature eigenvalues совпадают в центре, `K_r−K_Ω=O(ℓ²)`. Odd radial derivatives of even scalars vanish. Координатные переменные с `1/R` используют аналитические parity limits, не деление на `R=0`. Любой nonzero `V(φ₀)` допустим только если constraint включает его и asymptotics согласованы.

### Внешняя граница / асимптотика

Для AF: `Λ_eff=Λ₅+κ₅²V(φ∞)=0`, scalar подходит к вакууму достаточно быстро, metric к 5D Minkowski; радиус конечной computational boundary `R_out` и causal domain должны записываться. Scalar boundary condition должна быть incoming-characteristic controlled/absorbing в выбранной gauge, а не фиксированный `φ=Π=0` на пересекающей границу волне. GH gravitational incoming characteristics также заданы constraint-preserving; если надёжная absorbing condition пока отсутствует, ограничить физический анализ domain of dependence до возврата отражения и учитывать массу/поток на boundary.

## 5. Масса, причинность, горизонты и инварианты

Определить `X=(∇R)²_h`. В 5D AF Einstein gravity без `Λ_eff` Misner–Sharp/Tangherlini mass:

\[
M_{MS}=\frac{3\Omega_3}{16\pi G_5}R^2(1-X),
\quad \Omega_3=2\pi^2,
\quad f=1-\frac{\mu}{R^2},
\quad M=\frac{3\pi}{8G_5}\mu.
\]

При `Λ_eff≠0` mass convention должен subtract background: для vacuum `f=1−μ/R²−Λ_eff R²/6`, удобно `M=3Ω₃R²(1−X−Λ_effR²/6)/(16πG₅)`. При динамическом `V(φ)` scalar potential не является константным vacuum; указывать mass aspect и matter integral отдельно, не применять эту subtraction локально без derived balance law.

Пусть `n^a` — будущая timelike normal к выбранной 2D slice, `s^a` — наружная unit spacelike вектор; `k_±=n±s`, так что `k_+·k_-=-2`. Тогда

\[
\theta_\pm=\frac3R k_\pm^a\nabla_aR,
\qquad
\theta_+\theta_-=-\frac{9X}{R^2}.
\]

Положительное rescaling каждого future null нормаля сохраняет sign, но меняет величину expansions. Trapped: оба `<0`; anti-trapped: оба `>0`; marginal: одна `=0`, при finite other. Product/`X=0` фиксирует marginal condition, но не отличает trapped от anti-trapped без обоих знаков. Apparent/marginal surfaces slicing-dependent; event horizon глобален.

Radial null characteristics — решения `ds²=h_abdx^adx^b=0` в развитой metric, не неизменные `r±t`. Causal reachability только внутри вычисленного domain; global causal disconnection требует future-complete causal chart/data, не выводится из конечной сетки.

Ricci scalar для этой модели:

\[
R_5=\frac{10}{3}\Lambda_5-\frac{2\kappa_5^2}{3}T,
\quad T=-\frac32(\nabla\phi)^2-5V,
\quad R_5=\frac{10}{3}\Lambda_5+\kappa_5^2(\nabla\phi)^2+\frac{10}{3}\kappa_5^2V.
\]

Это проверка trace equation, но не полный набор curvature invariants. Валидировать через независимый Riemann reconstruction как минимум `R_ABR^AB`, `R_ABCDR^ABCD` (с оговоркой, что Lorentzian invariant не positive definite), mass aspect и null expansions. Большой координатный компонент не равен инвариантной singularity.

## 6. Какие компоненты требуют расширения модели

- Несперические perturbations, nonspherical horizons, вращение, anisotropic collapse — multidimensional grid либо отдельная математически выведенная perturbation system; radial SO(4) solver их не содержит.
- Физическая 3+1D daughter brane — embedding `X^A(ξ)`, induced metric, surface action, codimension-one junction/formation law; не заменяется scalar `R(t,r)`.
- ECSK spin source — независимая connection и state-resolved 5D matter stress; вывод 5D ECSK контакта–133 не дают отрицательной closure.
- EGB — отдельный вариационный принцип и generalized junctions; не добавлять его к этим equations.
