![Version:1.0.0](https://img.shields.io/static/v1?label=Version&message=1.0.0&color=blue&?style=plastic)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/last-week/AutoCopula)](https://cran.r-project.org/package=AutoCopula)

<img src="https://raw.githubusercontent.com/AdrianAntico/AutoCopula/master/inst/Logo.PNG" align="center" width="800" />

# **AutoCopula**
Automated Copula Modeling and Exploratory Data Analysis in R


[Download the AutoCopula Reference Manual](docs/AutoCopula-manual.pdf)


---

## **Overview**

AutoCopula is an R package designed for automating copula-based modeling, exploratory data analysis (EDA), and multivariate simulations. Whether you're working in finance, insurance, or any field requiring advanced dependency modeling, AutoCopula simplifies complex tasks through its intuitive API and advanced visualization tools.

---

## **Key Features**

- **Copula Modeling**:
  - Supports 26 copula models, (e.g. Gaussian, tCopula, Clayton, Gumbel, Frank, Joe, etc.).
  - Includes both bivariate and multivariate copulas.
  - Handles rotated and asymmetric copula models.
  - Automated model evaluation with metrics like AIC, BIC, and log-likelihood.

- **Exploratory Data Analysis (EDA)**:
  - Automated correlation analysis.
  - Interactive visualizations using `echarts4r`.
  - Empirical copula plots, scatterplots, and 3D density plots.

- **Scoring and Simulation**:
  - Generate conditional and hybrid predictions based on fitted copula models.
  - Perform large-scale simulations for stress testing and scenario analysis.
  - Interactive and back-transformed predictions for real-world interpretability.

### Conditional elliptical laws

Gaussian and Student-t conditional simulation is performed in latent
elliptical space and transformed back to copula uniforms before marginal
inversion. Student-t conditioning uses the fitted degrees of freedom, the
number of conditioned variables, and the conditional scale adjustment. This
preserves fitted tail behavior and keeps values in `(0, 1)` before marginal
back-transformation.

---

## **Installation**

### From GitHub

To install the development version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")
install.packages("R6")
install.packages("data.table")
install.packages("dplyr")
install.packages("echarts4r")
install.packages("copula")
install.packages("VineCopula")
install.packages("future")
install.packages("future.apply")

# Install AutoCopula
devtools::install_github("AdrianAntico/AutoCopula")
```

## Code Usage

### Preferred expert API

The preferred workflow is `copula_fit()` -> `copula_simulate()` or
`copula_conditional()` -> `copula_diagnose()`. The original R6 classes remain
available for compatibility and low-level mutable workflows.

```r
library(AutoCopula)
library(data.table)

set.seed(42)
z <- matrix(rnorm(2000), ncol = 2L)
data <- data.table(
  demand = z[, 1L],
  cost = 0.7 * z[, 1L] + sqrt(1 - 0.7^2) * z[, 2L]
)

fit <- copula_fit(
  data,
  families = c("Gaussian", "tCopula", "Clayton", "Gumbel", "Frank")
)

draws <- copula_simulate(
  fit,
  family = "tCopula",
  n = 5000L,
  seed = 42L
)

conditional <- copula_conditional(
  fit,
  family = "tCopula",
  known_ranges = list(demand = c(-1, 0, 1)),
  n = 1000L,
  seed = 42L
)

diagnostics <- copula_diagnose(fit)
```

Family identity is never erased. Gaussian, Student-t, Archimedean, extreme
value, BB, rotated, and Tawn families remain selectable by name. Advanced users
can extend or override a family through `family_definitions`, supplying an
engine-native `fit_function(data)` with its consequential controls captured in
the closure. Conditional simulation uses each fitted family's qualified law;
unsupported combinations fail rather than silently substituting independence.

### EDA

```r
library(AutoCopula)
library(MASS)

means <- c(0, 2, -1)
cov_matrix <- matrix(c(
  1.0, 0.9, 0.8,  # Correlations for Var1
  0.9, 1.0, 0.9,  # Correlations for Var2
  0.8, 0.9, 1.0   # Correlations for Var3
), nrow = 3, byrow = TRUE)

# Generate multivariate normal data
correlated_data <- MASS::mvrnorm(n = 1000, mu = means, Sigma = cov_matrix)

# Convert to data.table
data <- data.table(
  Var1 = correlated_data[, 1],
  Var2 = correlated_data[, 2],
  Var3 = correlated_data[, 3]
)

# Initialize EDA class
eda <- EDA$new(data)

# Generate correlation matrix
eda$correlate()

# Visualize distributions
eda$visualize_distributions()

# Generate scatterplots
eda$visualize_scatterplots()

# Generate a 3D density plot
eda$generate_3d_density_plot("Var1", "Var2", "Var3")

```

### Model Fitting

```r
# Initialize ModelFitter
fitter <- ModelFitter$new(data)

# List available copula models and their info
fitter$list_models()

# List just the model names
fitter$list_models()$Model

# Fit selected copula models
fitter$fit_models(c("Gaussian", "Clayton", "Frank"))
```

### Model Evaluation

```r
# Initialize ModelEvaluation
evaluator <- ModelEvaluation$new(fit_results = fitter$fit_results, data = data)

# Generate evaluation metrics
metrics <- evaluator$generate_metrics()

# Generate overlay plots
plots <- evaluator$generate_overlay_plot(model_name = "Gaussian")
```

### Model Scoring

```r
# Initialize ModelScorer
scorer <- ModelScorer$new(fit_results = fitter$fit_results, data = data)

# Single instance prediction
scorer$single_instance_prediction(model_name = "Gaussian")

# Batch predictions
scorer$batch_prediction(model_name = "Gaussian", n = 100)

# Large scale simulation
large_scale_pred <- scorer$large_scale_simulation(
  model_name = "Gaussian",
  batches = 1000,
  batch_size = 1000,
  parallel = FALSE,
  threads = 16
)

# Conditional range predictions
scorer$conditional_range_prediction(
  model_name = "Gaussian",
  known_ranges = list(
    Var1 = c(0.1, 0.2),
    Var2 = c(0.5, 1.5, 2.5)
  ),
  n = 10000,
  parallel = FALSE,
  threads = 2
)

# Get both conditional and unconditional predictions
hybrid_pred <- scorer$hybrid_range_simulation(
  model_name = "Gaussian",
  known_ranges = list(
    Var1 = seq(0.1, 1.0, by = 0.1),
    Var2 = c(0.5, 1.5, 2.5)
  ),
  n = 100,
  parallel = FALSE,
  threads = 2
)
```

## Dependencies
AutoCopula relies on the following R packages:

- R6
- data.table
- dplyr
- echarts4r
- copula
- VineCopula
- future
- future.apply

The Shiny App relies on the following R packages:
- shiny
- bs4Dash
- readxl
- DT

## Contributing
We welcome contributions! If you'd like to contribute, please:

1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request.

For bugs or feature requests, please open an issue on https://github.com/AdrianAntico/AutoCopula/issues.


## License
This project is licensed under the AGPL-3.0 License with additional conditions.
