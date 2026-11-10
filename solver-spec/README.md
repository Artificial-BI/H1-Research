# Solver H1 v0.2 — пакет спецификации

Этот каталог — **только переносимый пакет требований и математических заметок** для будущего самостоятельного проекта. Исходников solver здесь нет; реализацию этим заданием не начинать. Исходный проект H1 Lab и его solver остаются архивом read-only.

## Решение о базовой модели

Первый CPU solver следует ограничить классическим 5D Einstein–canonical-scalar действием, регулярными SO(4)-симметричными Cauchy-данными и проверенной 1+1 reduced geometry. Это позволит тестировать стандартную 5D GR/скалярную динамику, но **не** моделирует новую brane, смену топологии или рождение отдельной дочерней вселенной. ECSK и EGB — разные действия и не входят в первую реализацию. Общий bulk/две brane — исследовательский компонент, не обещанная функция.

Первый будущий implementation-issue (только после Gate A) должен содержать headless 1+1 Einstein–scalar ядро, constraint-solved initial data, immutable manifest, Minkowski/Tangherlini/linear-scalar/тест регулярных начальных данных controls и convergence/flux diagnostics; GUI, AMR, GPU, ECSK, EGB, branes и daughter-birth claims исключить. Условия его выдачи: аналитическая редукция/эволюционная формулировка и strong-hyperbolicity/constraint-propagation gates закрыты; интерфейсы и acceptance profile зафиксированы. Этот пакет сам по себе не разрешает писать solver.

## Порядок чтения

1. [Физическая спецификация](physics-specification.md)
2. [Математические уравнения](mathematical-equations.md)
3. [Численные методы](numerical-methods.md)
4. [Архитектура ПО](software-architecture.md)
5. [Валидация и диагностики](validation-and-diagnostics.md)
6. [Дорожная карта](development-roadmap.md)
7. [Открытые физические вопросы](open-physics-questions.md)
8. [GH-IBVP/CPBC: исходный аудит](gh-ibvp-cpbc-proof.md)
9. [Линеаризованные GH CPBC в 5D SO(4): характеристический вывод и пределы доказательства](gh-linearized-cpbc-characteristics.md)

## Трассировка требования → уравнение → проверка → источник

| Требование | Уравнение/раздел | Обязательная проверка | Локальная опора |
|---|---|---|---|
| Базовое действие и AF вакуум | § equations 1–3 | flat/vacuum exact solution; potential vacuum check | `docs/h1-v0.2-regular-initial-data-birth-test.md`; тест регулярных начальных данных |
| Регулярные SO(4) Cauchy data | § equations 5–7 | constraints, parity, origin limits, тест регулярных начальных данных reproduction | тест регулярных начальных данных, критерии коллапса и trapped surfaces |
| 5D mass normalization | § equations 8 | Tangherlini mass and invariant `∇R·∇R` | диагностика Type B в 5D; `docs/experiment-double-null-reference.md` |
| Trapping diagnosis | § equations 9 | independent null expansions and horizon convergence | диагностика Type B в 5D, критерии коллапса и trapped surfaces |
| Constraint propagation | numerical methods + validation | independent H/M residuals under refinement | #40–41F; отчёт о тестах коллапса–121 |
| Outgoing energy budget | validation | mass change + measured boundary flux + remainder | #41B–F; отчёт о тестах коллапса–121 |
| Daughter geometry and clocks | physics + open questions | global causal construction; worldline/event map | #31, #34–38, анализ общего гравитационного сектора–130 |
| ECSK spin closure | open questions | renormalized state correlator/stress; no current result | вывод 5D ECSK контакта report; анализ микроскопического спинового состояния report |

| GH CPBC and regular center | [линейный аудит CPBC](gh-linearized-cpbc-characteristics.md) | frozen 5D SO(4) boundary rank and energy estimate; variable-coefficient and regular-center proofs still required | restricted derivation; no nonlinear or SBP-SAT proof |

## Происхождение требований и пробелы

Локально найдены/сверены четыре документа математический реестр H1 v0.2, документы сравнение моделей-кандидатов–132 и доступный отчёт анализ микроскопического спинового состояния `docs/h1-v0.2-ecsk-5d-spin-state-closure.md`: state-dependent correlator не вычислен, знак spin stress не установлен. Архивные `experiment-31.md`–`experiment-41.md`, `experiment-41b-solver-repair.md`–`experiment-41f-compatible-operators.md`, `experiment-collapse-benchmark.md`, `experiment-double-null-reference.md` и `experiment-equations-cross-audit.md` доступны как предупреждения/история, не эталон solver. Отдельного локального файла, однозначно помеченного как предшествующее обсуждение, не найдено; equation cross-audit не объявляется его эквивалентом. Внешние версии GRChombo, GRTeclyn, Einstein Toolkit, LEAN, CUDA toolkit и C++ компиляторов **не проверялись**; ничего не скачивалось и не устанавливалось.

**Не выполнено:** кодирование, компиляция, симуляция, установка пакетов, тестовый прогон, GUI или глобальная causal reconstruction.
