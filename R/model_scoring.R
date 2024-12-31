#' ModelScorer
#'
#' An R6 class to score copula models on new data, perform predictions,
#' generate simulations, and visualize the results.
#' @export
ModelScorer <- R6::R6Class(
  "ModelScorer",
  public = list(
    #' @field fit_results A list of fitted copula model objects.
    fit_results = NULL,

    #' @field scored_data A list of scored results (e.g., predictions, simulations).
    scored_data = list(),

    #' @field data Scoring data
    data = list(),

    #' @field score_plots A list of plots visualizing scored data.
    score_plots = list(),

    #' @param fit_results A list of fitted copula model objects.
    #' @param data data for scoring a single instance of copula models
    #' @return A new instance of the ModelScorer class.
    initialize = function(fit_results, data) {
      if (!is.list(fit_results)) stop("fit_results must be a list of copula model objects.")
      if (!data.table::is.data.table(data)) stop("data must be a data.table.")
      self$fit_results <- fit_results
      self$data <- data  # Store the data for column names and other uses
    },

    #' @description Generate a single instance prediction using the entire model.
    #' @param model_name Name of the model to use.
    #' @return A single joint prediction as a named list.
    single_instance_prediction = function(model_name) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        pred <- copula::rCopula(1, fit@copula)
        colnames(pred) <- colnames(self$data)
        as.vector(pred)
      }, error = function(e) {
        message("Error in single instance prediction for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Generate batch predictions using the entire model.
    #' @param model_name Name of the model to use.
    #' @param n Number of instances to generate.
    #' @return A data.table of batch predictions.
    batch_prediction = function(model_name, n = 100) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        pred <- copula::rCopula(n, fit@copula)
        colnames(pred) <- colnames(self$data)
        result <- data.table::as.data.table(pred)
        self$scored_data[[model_name]] <- result
        return(result)
      }, error = function(e) {
        message("Error in batch prediction for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Perform large-scale simulations using fitted copula models.
    #' The simulation can be processed sequentially or in parallel, with user-defined batch sizes and thread counts.
    #'
    #' @param model_name A string specifying the name of the fitted model to use.
    #' @param batches An integer specifying the number of batches to process. Default is 10.
    #' @param batch_size An integer specifying the number of samples per batch. Default is 1000.
    #' @param parallel A logical value indicating whether to use parallel processing. Default is FALSE.
    #' @param threads An integer specifying the number of threads to use for parallel processing.
    #' If NULL, defaults to the number of available cores minus one.
    #'
    #' @return A `data.table` containing all simulated values, with a `BatchID` column to indicate the batch each row belongs to.
    #'
    #' @details This method generates large-scale simulations by splitting the total number of samples
    #' into smaller batches. If `parallel` is set to TRUE, it uses the `future.apply` package to
    #' process batches in parallel. On sequential processing, batches are processed one at a time.
    #' The function ensures that all parallel sessions are closed upon completion or in the case of an error.
    #'
    #' @export
    large_scale_simulation = function(model_name, batches = 10, batch_size = 1000, parallel = FALSE, threads = NULL) {
      fit <- self$fit_results[[model_name]]

      # Function to process a single batch
      process_batch <- function(batch_id) {
        message(sprintf("Processing batch %d", batch_id))
        batch <- copula::rCopula(batch_size, fit@copula)
        colnames(batch) <- colnames(self$data)
        batch_dt <- data.table::as.data.table(batch)
        batch_dt[, BatchID := batch_id]
        return(batch_dt)
      }

      tryCatch({
        if (parallel) {
          # Set the number of threads
          num_threads <- if (is.null(threads)) future::availableCores() - 1 else min(threads, future::availableCores())
          message(sprintf("Using %d threads for parallel processing.", num_threads))

          # Configure future backend
          future::plan(future::multisession, workers = num_threads)

          # Parallel processing with future_lapply
          results <- future.apply::future_lapply(seq_len(batches), process_batch)
        } else {
          # Sequential processing
          results <- lapply(seq_len(batches), process_batch)
        }

        # Combine all batches into one data.table
        final_result <- data.table::rbindlist(results)
        self$scored_data[[paste0(model_name, "_large_simulation")]] <- final_result
        return(final_result)
      }, error = function(e) {
        message("Error in large-scale simulation for model '", model_name, "': ", e$message)
        NULL
      }, finally = {
        # Ensure that all workers are shut down
        if (parallel) {
          message("Closing all parallel workers.")
          future::plan("sequential")
        }
      })
    },

    #' @description Generate conditional predictions for a subset of variables.
    #' @param model_name Name of the model to use.
    #' @param known_values A named list of known variable values.
    #' @return A data.table of conditional predictions.
    conditional_prediction = function(model_name, known_values) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        # Perform conditional sampling (to be implemented)
        # Placeholder logic: Replace with real conditional copula sampling
        conditioned_data <- copula::conditionalCopula(fit@copula, known_values)
        data.table::as.data.table(conditioned_data)
      }, error = function(e) {
        message("Error in conditional prediction for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Generate hybrid simulations (conditional + unconditional).
    #' @param model_name Name of the model to use.
    #' @param known_values A named list of known variable values.
    #' @param n Number of instances to simulate.
    #' @return A data.table of hybrid simulations.
    hybrid_simulation = function(model_name, known_values, n = 100) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        conditional <- copula::conditionalCopula(fit@copula, known_values)
        unconditional <- copula::rCopula(n, fit@copula)
        list(Conditional = conditional, Unconditional = unconditional)
      }, error = function(e) {
        message("Error in hybrid simulation for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Perform importance sampling for the copula model.
    #' @param model_name Name of the model to use.
    #' @param n Number of samples to draw.
    #' @return Weighted samples.
    importance_sampling = function(model_name, n = 1000) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        samples <- copula::rCopula(n, fit@copula)
        weights <- exp(-abs(samples))  # Example weighting; adjust as needed
        list(Samples = samples, Weights = weights)
      }, error = function(e) {
        message("Error in importance sampling for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Perform stress testing by simulating extreme events.
    #' @param model_name Name of the model to use.
    #' @param n Number of extreme event scenarios to simulate.
    #' @return A data.table of extreme event simulations.
    stress_testing = function(model_name, n = 100) {
      fit <- self$fit_results[[model_name]]
      tryCatch({
        extreme_events <- copula::rCopula(n, fit@copula)
        extreme_events <- pmax(extreme_events, 0.95)  # Simulate upper tail
        data.table::as.data.table(extreme_events)
      }, error = function(e) {
        message("Error in stress testing for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Visualize single instance prediction.
    #' @param predictions A named list of single instance predictions.
    #' @return An echarts4r bar chart visualization.
    visualize_single_instance_prediction = function(predictions) {
      pred_dt <- data.table::as.data.table(predictions)
      pred_dt[, Variable := names(predictions)]
      pred_dt |>
        echarts4r::e_charts(Variable) |>
        echarts4r::e_bar(value = V1, name = "Prediction") |>
        echarts4r::e_title("Single Instance Prediction") |>
        echarts4r::e_theme("dark")
    },

    #' @description Visualize batch predictions.
    #' @param batch_predictions A data.table of batch predictions.
    #' @return An echarts4r scatter plot visualization.
    visualize_batch_predictions = function(batch_predictions) {
      batch_predictions |>
        data.table::as.data.table() |>
        echarts4r::e_charts(Var1) |>
        echarts4r::e_scatter(Var2, name = "Var2") |>
        echarts4r::e_scatter(Var3, name = "Var3") |>
        echarts4r::e_title("Batch Predictions") |>
        echarts4r::e_theme("dark")
    },

    #' @description Visualize large-scale simulation.
    #' @param large_scale_predictions A data.table of large-scale predictions.
    #' @return An echarts4r density plot visualization.
    visualize_large_scale_simulation = function(large_scale_predictions) {
      large_scale_predictions |>
        data.table::as.data.table() |>
        echarts4r::e_charts(Var1) |>
        echarts4r::e_density(name = "Var1 Density") |>
        echarts4r::e_title("Large-Scale Simulation Density for Var1") |>
        echarts4r::e_theme("dark")
    },

    #' @description Visualize conditional predictions.
    #' @param conditional_predictions A data.table of conditional predictions.
    #' @return An echarts4r scatter plot visualization.
    visualize_conditional_predictions = function(conditional_predictions) {
      conditional_predictions |>
        data.table::as.data.table() |>
        echarts4r::e_charts(Var1) |>
        echarts4r::e_scatter(Var2, name = "Conditioned Var2") |>
        echarts4r::e_scatter(Var3, name = "Conditioned Var3") |>
        echarts4r::e_title("Conditional Predictions") |>
        echarts4r::e_theme("dark")
    },

    #' @description Visualize hybrid simulations.
    #' @param hybrid_predictions A list with conditional and unconditional components.
    #' @return A list of echarts4r density plots.
    visualize_hybrid_simulation = function(hybrid_predictions) {
      conditional <- hybrid_predictions$Conditional |> data.table::as.data.table()
      unconditional <- hybrid_predictions$Unconditional |> data.table::as.data.table()

      plots <- list(
        Conditional = conditional |>
          echarts4r::e_charts(Var1) |>
          echarts4r::e_density(name = "Conditional Density") |>
          echarts4r::e_title("Hybrid Simulation - Conditional Component") |>
          echarts4r::e_theme("dark"),

        Unconditional = unconditional |>
          echarts4r::e_charts(Var1) |>
          echarts4r::e_density(name = "Unconditional Density") |>
          echarts4r::e_title("Hybrid Simulation - Unconditional Component") |>
          echarts4r::e_theme("dark")
      )
      return(plots)
    },

    #' @description Visualize importance sampling.
    #' @param importance_samples A list with samples and weights.
    #' @return An echarts4r bubble plot visualization.
    visualize_importance_sampling = function(importance_samples) {
      samples <- data.table::as.data.table(importance_samples$Samples)
      weights <- importance_samples$Weights
      samples[, Weight := weights]
      samples |>
        echarts4r::e_charts(Var1) |>
        echarts4r::e_bubble(Var2, Var3, Weight, name = "Weighted Samples") |>
        echarts4r::e_title("Importance Sampling Visualization") |>
        echarts4r::e_theme("dark")
    },

    #' @description Visualize stress testing results.
    #' @param stress_test_results A data.table of stress test results.
    #' @return An echarts4r density plot visualization for extreme events.
    visualize_stress_testing = function(stress_test_results) {
      stress_test_results |>
        data.table::as.data.table() |>
        echarts4r::e_charts(Var1) |>
        echarts4r::e_density(name = "Extreme Event Density") |>
        echarts4r::e_title("Stress Testing - Extreme Events") |>
        echarts4r::e_theme("dark")
    }
  )
)
