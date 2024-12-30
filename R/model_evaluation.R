#' ModelEvaluation
#'
#' An R6 class to evaluate copula models.
#' Includes tools to generate tables of evaluation metrics and visualizations
#' to assess model performance and dependence structure.
#'
#' @export
ModelEvaluation <- R6::R6Class(
  "ModelEvaluation",
  public = list(

    #' @field fit_results A list of fitted copula model objects.
    fit_results = NULL,

    #' @field evaluation_metrics A data.table containing model evaluation metrics.
    evaluation_metrics = NULL,

    #' @field plots A list of visualizations for evaluating copula models.
    plots = list(),

    #' @field data The original dataset used for fitting copulas.
    data = NULL,

    #' Initialize the ModelEvaluation class
    #'
    #' @param fit_results A list of fitted copula model objects.
    #' @param data The dataset used for fitting models, in uniform margins.
    #' @return A new instance of the ModelEvaluation class.
    initialize = function(fit_results, data) {
      if (!is.list(fit_results)) stop("fit_results must be a list of copula model objects.")
      if (!data.table::is.data.table(data)) stop("data must be a data.table.")
      self$fit_results <- fit_results
      self$data <- data
    },

    #' @description Compute and summarize evaluation metrics for fitted copulas.
    #' Metrics include log-likelihood, AIC, BIC, and dependence measures.
    #'
    #' @return A data.table of evaluation metrics.
    generate_metrics = function() {
      if (is.null(self$fit_results)) stop("No fitted models to evaluate.")

      # Compute static metrics for each model
      metrics <- lapply(names(self$fit_results), function(model_name) {
        fit <- self$fit_results[[model_name]]
        if (is.null(fit)) return(NULL)

        tryCatch({
          copula_model <- fit@copula  # Extract the actual copula object

          log_likelihood <- logLik(fit)
          aic <- AIC(fit)
          bic <- BIC(fit)
          tau <- copula::tau(copula_model)

          # Only compute rho if the method is available
          rho <- tryCatch({
            copula::rho(copula_model)
          }, error = function(e) NA)

          tail_dependence <- if (inherits(copula_model, "archmCopula")) {
            list(lambda_lower = copula::lambda(copula_model)[1],
                 lambda_upper = copula::lambda(copula_model)[2])
          } else {
            list(lambda_lower = NA, lambda_upper = NA)
          }

          list(
            `Model Name` = model_name,
            LogLik = as.numeric(log_likelihood),
            AIC = aic,
            BIC = bic,
            Kendall_Tau = tau,
            Spearman_Rho = rho,
            Lambda_Lower = tail_dependence$lambda_lower,
            Lambda_Upper = tail_dependence$lambda_upper
          )
        }, error = function(e) {
          message("Error processing model: ", e$message)
          NULL
        })
      })

      # Combine metrics into a data.table
      self$evaluation_metrics <- data.table::rbindlist(metrics, fill = TRUE)
      return(self$evaluation_metrics)
    },

    #' @description Generate overlay plot for observed vs simulated pseudo-observations.
    #' @param model_name The name of the model for which to generate the plot.
    #' @param var the variable you want to compare to simulated values
    #' @param theme name of theme for `echarts4r`
    #' @return An echarts4r plot object.
    generate_overlay_plot = function(model_name, var, theme = "westeros") {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        # Generate pseudo-observations from the data
        pseudo_obs <- copula::pobs(as.matrix(self$data))
        colnames(pseudo_obs) <- colnames(self$data)
        pseudo_obs_dt <- data.table::as.data.table(pseudo_obs)

        # Simulate pseudo-observations from the fitted copula
        simulated_obs <- copula::rCopula(nrow(self$data), fit@copula)
        colnames(simulated_obs) <- colnames(self$data)
        simulated_obs_dt <- data.table::as.data.table(simulated_obs)

        # Stack the observed and simulated data
        stacked_data <- rbindlist(lapply(colnames(pseudo_obs_dt), function(col) {
          data.table(
            Observed = sort(pseudo_obs_dt[[col]]),
            Simulated = sort(simulated_obs_dt[[col]]),
            Variable = col
          )
        }))

        # Create overlay plot
        echarts4r::e_charts(data = stacked_data |> dplyr::group_by(Variable), x = Simulated) |>
          echarts4r::e_scatter(serie = Observed) |>
          echarts4r::e_title(text = paste("Observed vs Simulated:", model_name)) |>
          echarts4r::e_theme(theme) |>
          echarts4r::e_legend(orient = "horizontal", top = 35, type = "scroll")
      }, error = function(e) {
        message("Error generating overlay plot for model '", model_name, "': ", e$message)
        NULL
      })
    }
  )
)
