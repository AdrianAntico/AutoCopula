# Expert API Reconstruction 2.0

AutoCopula keeps the four memorable verbs and restores expert-grade vine
methodology around them. Clean wrappers are not the product; pair-specific
nonexchangeable dependence is.

## Export counts

| | Count | Symbols |
| --- | ---: | --- |
| Before | 9 | `copula_fit`, `copula_simulate`, `copula_conditional`, `copula_diagnose`, `EDA`, `ModelFitter`, `ModelEvaluation`, `ModelScorer`, `qa_autocopula_conditional_elliptical` |
| After | 5 | `copula_fit`, `copula_simulate`, `copula_conditional`, `copula_diagnose`, `copula_families` |

`NAMESPACE` is hand-curated. DESCRIPTION disables the roxygen namespace roclet:

```
Roxygen: list(markdown = TRUE, roclets = c("collate", "rd"))
```

## Public surface

CORE_EXPERT:

- `copula_fit`
- `copula_simulate`
- `copula_conditional`
- `copula_diagnose`

SPECIALIZED_EXPERT:

- `copula_families`

Unexported, still present for internal use and `AutoCopula:::ModelFitter`:

- `EDA`
- `ModelFitter`
- `ModelEvaluation`
- `ModelScorer`
- `qa_autocopula_conditional_bicop`
- `qa_autocopula_conditional_elliptical`

## CAPABILITY LOST?

No. Vines were added. The four verbs still fit, simulate, condition, and
diagnose Gaussian, t, Archimedean, extreme-value, BB, rotated, and Tawn
families. R6 classes remain in the package; they are just not preferred
exports. `qa_*` qualification functions remain callable via `:::`.

Capability gained:

- First-class `"RVine"` / `"vine"` / `"CVine"` / `"DVine"` in `copula_fit()`
- Expert vine controls: `family_set`, `trunc_lvl`, `treecrit`, `vine_type`
- Fitted `RVineMatrix` stored on the model object
- `VineCopula::RVineSim` simulation
- Vine conditionals via h-functions / `RVineCondSim` when applicable, else
  Monte Carlo kernel SIR with recorded provenance (not a fake analytic law)
- Vine diagnostics: pair-family table, AIC/BIC, truncation, tail notes
- `copula_families()` catalog so experts do not call
  `ModelFitter$new(data)$list_models()`

## New `copula_fit` formals

Positional `data` and `families` are unchanged.

```
copula_fit(
  data,
  families = c("Gaussian", "tCopula"),
  family_set = NULL,
  trunc_lvl = NA_integer_,
  treecrit = "tau",
  vine_type = "rvine",
  family_definitions = list()
)
```

Default `family_set` is not Gaussian+t. It includes Gaussian, t, Clayton,
Gumbel, Frank, Joe, BB1, BB7, and survival (180-degree) rotations.

## Nonexchangeable dependence

Vines are the multivariate nonexchangeable path. Each pair and conditional
pair may have its own family and tail coefficients.

`d > 2` Clayton, Gumbel, Frank, and Joe objects from the copula package remain
exchangeable. That is documented on those catalog rows. Experts who need
pair-specific tails should fit `RVine` / `CVine` / `DVine`.

## Family catalog

`copula_families()` returns the built-in library, including vines. Do not
construct `ModelFitter` only to list names.
