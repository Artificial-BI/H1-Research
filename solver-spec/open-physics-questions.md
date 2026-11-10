# Открытые физические вопросы и блокеры

| Вопрос/требование H1 | Статус | Что известно и что требуется |
|---|---|---|
| Минимальная 5D Einstein–canonical-scalar local dynamics | **DERIVED** | Действие, field equations, SO(4) reduction и constraints записаны в [equations](mathematical-equations.md); specific numerical formulation ещё требует hyperbolicity/boundary proof. |
| тест регулярных начальных данных smooth one-ended regular data | **DERIVED in restricted model** | Massless scalar, time symmetric, `Q L²<1` sufficient for no initial trapped/anti-trapped round S³. Не доказывает future collapse. |
| Type B 4D mechanism переносится в 5D | **REJECTED as automatic transfer** | Ning et al. use pre-existing 4D false-vacuum domain/wall. 5D initial data, matter content and mass/expansion coefficients differ. |
| Старый polar-areal/FV/double-null solver — validated oracle | **REJECTED** | #40–41F/отчёт о тестах коллапса–121 show independent constraint and mass-ledger issues; отчёт о double-null проверке coupled null reference also partial. Preserve archive, do not call ground truth. |
| Full nonlinear 5D SO(4) GH Cauchy/IBVP with CPBC | **CONDITIONAL / NOT PROVEN** | The frozen Minkowski SO(4) outer-boundary principal system now has an explicit rank-2 constraint map, two residual gauge data, and a non-increasing homogeneous wave-energy flux. Variable-coefficient/nonlinear CPBC, initial-corner compatibility, and weighted regular-center closure remain open. See [restricted linear CPBC derivation](gh-linearized-cpbc-characteristics.md) and [GH-IBVP audit](gh-ibvp-cpbc-proof.md). |
| 5D ECSK contact term | **DERIVED** | вывод 5D ECSK контакта action-level 3-form current operator. |
| Sign/scaling of 5D ECSK macroscopic spin stress | **NOT_DERIVED** | анализ микроскопического спинового состояния report exists; specific `C_{ABC,DEF}`, renormalization, state-dependent metric variation and `ρ,p` not computed. |
| Spin-fluid `−C_s a⁻⁸` | **CONDITIONAL** | Requires frozen pair correlations and negative effective stress; not consequence of вывод 5D ECSK контакта or анализ микроскопического спинового состояния. |
| ECSK high-density stability/hyperbolicity | **NOT_DERIVED** | Need state-specific second variation, perturbation modes, EFT cutoff and gravitational kinetic sign. |
| 5D EGB regularizes Type B | **NOT_DERIVED** | Different theory; full branch/hyperbolicity and same-system solution absent. |
| Single SO(4)-symmetric 5D region is a 3+1D daughter brane | **REJECTED as an identification** | SO(4) orbits are S³ in a 5D bulk; brane needs embedding and localized action. |
| Formation from regular parent-only data of new topology/component | **NOT_DERIVED** | No complete action+global initial-value solution. Bubble seeds in Type B were already present initially. |
| Shared bulk gravity between distinct sectors | **CONDITIONAL** | A common bulk action can mediate causal response within its manifold; independent disconnected manifolds need a separately defined mediator/correspondence. |
| Parent asymptotic mass necessarily retained/lost at formation | **NOT_DERIVED** | Depends on bulk flux, boundary conditions, asymptotic charge and whether extra sector radiates. |
| `T_P`, `τ_trans`, `t_D` map / `dt_D/dT_P` | **NOT_DERIVED** | Must identify worldlines and events from solution; static RS redshift is not a cosmological clock law. |
| H2 dark matter/acceleration | **OUT OF SCOPE** | Separate hypothesis; not included in action or roadmap. |
| Repeated daughter births from one BH | **NOT_DERIVED** | Not an axiom and no mechanism/action specified. |

## Наиболее опасная предпосылка

Самая опасная нечисленная предпосылка — что в одной глобально регулярной и causal theory можно получить новую independently evolving geometry из parent collapse, сохранив admissible initial data/energy conditions и определив event map между уровнями. Локальный bounce или общая гравитационная law не создаёт вторую geometry сам по себе.

## Что может фальсифицировать конкретную версию H1

После фиксации единого action и observable map: отсутствие регулярного constraint-compatible solution в заявленном domain; невыполнение условий на trapped/anti-trapped transitions; отсутствие причинно доступной области/предсказанного response; конфликт рассчитанного сигнала с измерениями; либо невозможность отличить clock/metric effect от coordinate reparametrization. Пока уникальных количественных H1 predictions и соответствующего map нет, экспериментальная фальсификация именно H1 недоопределена.

## Что не должно превращаться в параметр по умолчанию

- `K₀a_D` time mapping;
- spin² sign or density scaling;
- Type B critical thresholds transplanted to 5D;
- black-hole exterior mass conservation;
- existence of brane from an areal-radius throat;
- “global daughter” status from finite-time simulation.
