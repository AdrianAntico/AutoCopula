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
        # Check if the model is a Vine Copula model
        if (inherits(fit, "BiCop")) {
          # Generate a single prediction using VineCopula's BiCopSim
          pred <- VineCopula::BiCopSim(1, family = fit$family, par = fit$par, par2 = fit$par2)
        } else if (!is.null(fit@copula)) {
          # Generate a single prediction using the copula package
          pred <- copula::rCopula(1, fit@copula)
        } else {
          stop("Unsupported copula type for model '", model_name, "'.")
        }

        # Assign column names if data is available
        if (!is.null(self$data)) {
          colnames(pred) <- colnames(self$data)
        }
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
        if (inherits(fit, "BiCop")) {
          # Generate batch predictions using VineCopula's BiCopSim
          pred <- VineCopula::BiCopSim(n, family = fit$family, par = fit$par, par2 = fit$par2)
        } else if (!is.null(fit@copula)) {
          # Generate batch predictions using the copula package
          pred <- copula::rCopula(n, fit@copula)
        } else {
          stop("Unsupported copula type for model '", model_name, "'.")
        }

        # Assign column names if data is available
        if (!is.null(self$data)) {
          colnames(pred) <- colnames(self$data)
        }

        # Convert to data.table and return
        result <- data.table::as.data.table(pred)
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

      tryCatch({
        if (parallel) {
          # Load the future.apply package
          require(future.apply)

          # Set the number of threads
          num_threads <- if (is.null(threads)) future::availableCores() - 1 else min(threads, future::availableCores())
          message(sprintf("Using %d threads for parallel processing.", num_threads))

          # Configure future backend
          future::plan(future::multisession, workers = num_threads)

          # Parallel processing with future_lapply, ensuring proper random number seeding
          results <- future.apply::future_lapply(seq_len(batches), private$process_batch, future.seed = TRUE)
        } else {
          # Sequential processing
          results <- lapply(seq_len(batches), private$process_batch)
        }

        # Combine all batches into one data.table
        final_result <- data.table::rbindlist(results)
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
    #'
    #' @details
    #' The `conditional_prediction()` function generates conditional samples from a fitted copula model
    #' based on specified known values for a subset of variables. This method is particularly useful in
    #' scenarios where some variables are observed (known), and predictions or simulations are required
    #' for the remaining variables, conditioned on the observed values.
    #'
    #' The function supports a wide range of copula families, including both multivariate and bivariate
    #' copulas, as well as rotated and asymmetric copula models. The appropriate private method is
    #' dynamically selected based on the copula type to handle the specifics of each family.
    #'
    #' - **Input Parameters**:
    #'   - `model_name`: The name of the fitted copula model from which predictions will be generated.
    #'   - `known_values`: A named list where the names correspond to variable names, and the values
    #'     represent the observed values for those variables.
    #'   - `n`: The number of conditional samples to generate.
    #'
    #' - **Process**:
    #'   1. Validate that the `model_name` corresponds to a fitted copula model and that the `known_values`
    #'      variables match the names of the dataset used for model fitting.
    #'   2. Convert the known values into pseudo-observations using their empirical cumulative
    #'      distribution functions (ECDFs).
    #'   3. Call the appropriate private conditional sampling method for the copula family, which computes
    #'      the conditional distribution of the remaining variables given the known values.
    #'   4. Combine the conditional samples with the known values to produce the final result.
    #'
    #' - **Output**:
    #'   - A `data.table` containing the conditional samples for the remaining variables, along with the
    #'     known values for the conditioned variables.
    #'
    #' - **Error Handling**:
    #'   If the copula type is not supported or if the input parameters are invalid, the function will
    #'   return `NULL` with an appropriate warning or error message.
    #'
    #' - **Applications**:
    #'   - Use cases include scenarios in finance, insurance, or any domain requiring predictions under
    #'     partially observed conditions.
    #'
    #' **Note**: The method currently supports conditioning on a subset of variables for both multivariate
    #' and bivariate copulas. Some copula families (e.g., rotated copulas) require special handling, which
    #' is implemented through the private methods.
    #'
    #' @param model_name Name of the model to use.
    #' @param known_values A named list of known variable values.
    #' @param n Number of samples to get
    #'
    #' @return A data.table of conditional predictions.
    #'
    #' @export
    #' @description Generate conditional predictions for a subset of variables.
    #'
    #' @details
    #' The `conditional_prediction()` function generates conditional samples from a fitted copula model
    #' based on specified known values or ranges for a subset of variables. This method is particularly
    #' useful in scenarios where some variables are observed (known) with exact values or are constrained
    #' within specific ranges, and predictions or simulations are required for the remaining variables.
    #'
    #' - **Input Parameters**:
    #'   - `model_name`: The name of the fitted copula model from which predictions will be generated.
    #'   - `known_values`: A named list where:
    #'     - **Single Values**: The value represents an exact observation of the variable (e.g., `X = 0.5`).
    #'       These are converted into pseudo-observations using the empirical cumulative distribution function (ECDF).
    #'     - **Ranges**: The value is a vector of length 2 (e.g., `X = c(0.4, 0.6)`), representing the range within
    #'       which the variable is constrained. Samples are drawn uniformly within the range during each iteration.
    #'   - `n`: The number of conditional samples to generate.
    #'
    #' - **Process**:
    #'   1. Validate that `model_name` corresponds to a fitted copula model and that `known_values` variables
    #'      match the dataset used for model fitting.
    #'   2. Process `known_values`:
    #'      - For single values, convert them into pseudo-observations using their ECDF.
    #'      - For ranges, sample uniformly within the range for each iteration.
    #'   3. Dynamically select the appropriate private conditional sampling method for the copula family.
    #'   4. Compute the conditional distribution of the remaining variables given the known values or sampled ranges.
    #'   5. Combine the conditional samples with the known values/ranges to produce the final result.
    #'
    #' - **Output**:
    #'   - A `data.table` containing the conditional samples for the remaining variables, along with the
    #'     known values or sampled values for the conditioned variables.
    #'
    #' - **Error Handling**:
    #'   - If `model_name` is not a supported copula type or `known_values` contains invalid inputs, the
    #'     function returns `NULL` with an appropriate warning or error message.
    #'
    #' - **Applications**:
    #'   - Use cases include finance, insurance, and other domains requiring predictions under partially
    #'     observed or constrained conditions.
    #'
    #' @param model_name Name of the model to use.
    #' @param known_values A named list of known variable values or ranges.
    #'   - For single values: Provide a scalar value (e.g., `X = 0.5`).
    #'   - For ranges: Provide a vector of length 2 (e.g., `X = c(0.4, 0.6)`).
    #' @param n Number of conditional samples to generate.
    #'
    #' @return A `data.table` of conditional predictions.
    #'
    #' @export
    conditional_prediction = function(model_name, known_values, n = 1) {
      fit <- self$fit_results[[model_name]]

      # Determine copula type
      if (model_name %in% c(
        "Gaussian",
        "tCopula",
        "Clayton",
        "Gumbel",
        "Frank",
        "Joe"
      )) {
        copula_model <- fit@copula
      } else {
        copula_model <- "vine"
      }

      tryCatch({
        # Ensure known_values is valid
        if (!all(names(known_values) %in% colnames(self$data))) {
          warning("Some known_values do not match column names in the dataset. Returning NULL.")
          return(NULL)
        }

        # Delegate to the appropriate private method based on copula type
        if (inherits(copula_model, "normalCopula")) {
          return(private$conditional_gaussian(fit, known_values, n))
        } else if (inherits(copula_model, "tCopula")) {
          return(private$conditional_t(fit, known_values, n))
        } else if (inherits(copula_model, "claytonCopula")) {
          return(private$conditional_clayton(fit, known_values, n))
        } else if (inherits(copula_model, "gumbelCopula")) {
          return(private$conditional_gumbel(fit, known_values, n))
        } else if (inherits(copula_model, "frankCopula")) {
          return(private$conditional_frank(fit, known_values, n))
        } else if (inherits(copula_model, "joeCopula")) {
          return(private$conditional_joe(fit, known_values, n))
        } else if (model_name %in% c("BB1", "BB6", "BB7", "BB8")) {
          conditional_function <- private[[paste0("conditional_", tolower(model_name))]]
          return(conditional_function(fit, known_values, n))
        } else if (grepl("Rotated", model_name)) {
          # Rotated copulas (e.g., Rotated Clayton (180))
          model_key <- gsub("[()\\s]", "_", tolower(model_name))
          conditional_function <- private[[paste0("conditional_", model_key)]]
          return(conditional_function(fit, known_values, n))
        } else {
          message("Conditional sampling is not implemented for this copula type.")
          return(NULL)
        }
      }, error = function(e) {
        message("Error in conditional prediction for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Generate conditional predictions for a range of known values.
    #'
    #' @details
    #' This function generates conditional predictions over a range of values for a single known variable.
    #' It iteratively calls `conditional_prediction()` for each value in the range, with an option to
    #' parallelize the computations for efficiency.
    #'
    #' - **Parallel Processing**:
    #'   - If `parallel = TRUE`, the function uses the `future.apply` package to process the range
    #'     values concurrently.
    #'   - The number of threads used can be controlled via the `threads` parameter.
    #'
    #' @param model_name Name of the model to use.
    #' @param known_variable The name of the variable to condition on.
    #' @param value_range A numeric vector specifying the range of values for the known variable.
    #' @param n Number of samples to generate for each value in the range.
    #' @param parallel Logical, whether to process the range values in parallel. Default is `FALSE`.
    #' @param threads Integer, the number of threads to use if `parallel = TRUE`. Defaults to all available cores minus one.
    #'
    #' @return A `data.table` containing conditional predictions for the specified range.
    #'
    #' @export
    conditional_range_prediction = function(model_name, known_ranges, n = 100, parallel = FALSE, threads = NULL) {
      # Validate inputs
      if (!is.list(known_ranges) || length(known_ranges) == 0) {
        stop("known_ranges must be a non-empty named list of variable ranges.")
      }
      if (!all(sapply(known_ranges, is.numeric))) {
        stop("All ranges in known_ranges must be numeric vectors.")
      }

      # Generate all combinations of known variable values
      combinations <- do.call(expand.grid, known_ranges)

      # Parallel or sequential processing
      results_list <- tryCatch({
        if (parallel) {
          # Load the future.apply package
          options(future.rng.onMisuse = "ignore")
          RNGkind("L'Ecuyer-CMRG")
          set.seed()

          # Set the number of threads
          num_threads <- if (is.null(threads)) future::availableCores() - 1 else min(threads, future::availableCores())
          message(sprintf("Using %d threads for parallel processing.", num_threads))
          future::plan(future::multisession, workers = num_threads)

          # Process combinations in parallel
          future.apply::future_lapply(seq_len(nrow(combinations)), function(i) {
            known_values <- as.list(combinations[i, , drop = FALSE])
            private$process_value_conditional(known_values = known_values, batch_id = i)
          }, future.seed = TRUE)
        } else {
          # Process combinations sequentially
          lapply(seq_len(nrow(combinations)), function(i) {
            known_values <- as.list(combinations[i, , drop = FALSE])
            private$process_value_conditional(known_values = known_values, batch_id = i)
          })
        }
      }, error = function(e) {
        message("Error during conditional range prediction: ", e$message)
        return(NULL)
      }, finally = {
        # Ensure all workers are shut down if parallel
        if (parallel) {
          future::plan("sequential")
        }
      })

      # Combine results into a single data.table
      if (length(results_list) > 0) {
        final_results <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
        return(final_results)
      } else {
        warning("No predictions were generated. Please check your inputs.")
        return(NULL)
      }
    },

    #' @description Generate hybrid simulations (conditional + unconditional).
    #' @param model_name Name of the model to use.
    #' @param known_values A named list of known variable values.
    #' @param n Number of instances to simulate.
    #' @return A list containing two `data.table` objects: one for conditional predictions
    #' and another for unconditional simulations.
    #' @export
    hybrid_simulation = function(model_name, known_values, n = 100) {
      tryCatch({

        # Perform conditional predictions
        conditional <- self$conditional_prediction(
          model_name = model_name,
          known_values = known_values,
          n = n
        )

        # Perform unconditional predictions
        unconditional <- self$batch_prediction(
          model_name = model_name,
          n = n
        )

        # Return both results
        list(Conditional = conditional, Unconditional = unconditional)
      }, error = function(e) {
        message("Error in hybrid simulation for model '", model_name, "': ", e$message)
        NULL
      })
    },

    #' @description Generate hybrid simulations for a range of known values.
    #'
    #' @details
    #' The `hybrid_range_simulation()` function extends `hybrid_simulation()` by enabling simulations
    #' over a range of values for one known variable. It generates both conditional and unconditional
    #' samples for each value in the specified range and combines the results.
    #'
    #' - **Input Parameters**:
    #'   - `model_name`: The name of the fitted copula model to use for simulations.
    #'   - `known_variable`: The name of the variable to condition on.
    #'   - `value_range`: A numeric vector specifying the range of values for the known variable.
    #'   - `n`: Number of instances to simulate for each known value.
    #'   - `parallel`: Logical indicating whether to run the simulation in parallel. Default is `FALSE`.
    #'   - `threads`: Number of parallel threads to use if `parallel = TRUE`. Defaults to available cores minus one.
    #'
    #' - **Output**:
    #'   - A `data.table` containing the hybrid simulation results for all values in the specified range.
    #'
    #' @param model_name Name of the model to use.
    #' @param known_variable Name of the variable to condition on.
    #' @param value_range A numeric vector specifying the range of known values.
    #' @param n Number of instances to simulate for each known value.
    #' @param parallel Logical, whether to use parallel processing.
    #' @param threads Number of threads to use for parallel processing. Defaults to available cores minus one.
    #'
    #' @return A data.table of hybrid simulations for the specified range of values.
    #'
    #' @export
    hybrid_range_simulation = function(model_name, known_ranges, n = 100, parallel = FALSE, threads = NULL) {
      # Validate inputs
      if (!is.list(known_ranges) || length(known_ranges) == 0) {
        stop("known_ranges must be a non-empty named list of variable ranges.")
      }
      if (!all(sapply(known_ranges, is.numeric))) {
        stop("All ranges in known_ranges must be numeric vectors.")
      }

      # Generate all combinations of known variable values
      combinations <- do.call(expand.grid, known_ranges)

      # Parallel or sequential processing
      results_list <- tryCatch({
        if (parallel) {
          # Load the future.apply package
          options(future.rng.onMisuse = "ignore")
          RNGkind("L'Ecuyer-CMRG")
          set.seed()

          # Set the number of threads
          num_threads <- if (is.null(threads)) future::availableCores() - 1 else min(threads, future::availableCores())
          message(sprintf("Using %d threads for parallel processing.", num_threads))
          future::plan(future::multisession, workers = num_threads)

          # Process combinations in parallel
          future.apply::future_lapply(seq_len(nrow(combinations)), function(i) {
            known_values <- as.list(combinations[i, , drop = FALSE])
            private$process_value_hybrid(value = known_values, batch_id = i)
          }, future.seed = TRUE)
        } else {
          # Process combinations sequentially
          lapply(seq_len(nrow(combinations)), function(i) { # i = 1
            known_values <- as.list(combinations[i, , drop = FALSE])
            private$process_value_hybrid(value = known_values, batch_id = i)
          })
        }
      }, error = function(e) {
        message("Error during hybrid range simulation: ", e$message)
        return(NULL)
      }, finally = {
        # Ensure all workers are shut down if parallel
        if (parallel) {
          future::plan("sequential")
        }
      })

      # Combine results into a single data.table
      if (length(results_list) > 0) {
        final_results <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)
        return(final_results)
      } else {
        warning("No simulations were generated. Please check your inputs.")
        return(NULL)
      }
    },

    #' @description Perform importance sampling for the copula model.
    #' @param model_name Name of the model to use.
    #' @param n Number of samples to draw.
    #' @return Weighted samples.
    #' @export
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
    #' @export
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
    }
  ),

  private = list(

    # Private method for processing a single batch in batch_predictions()
    process_batch = function(batch_id) {
      message(sprintf("Processing batch %d", batch_id))

      # Generate batch data based on the copula type
      if (inherits(fit, "BiCop")) {
        batch <- VineCopula::BiCopSim(batch_size, family = fit$family, par = fit$par, par2 = fit$par2)
      } else if (!is.null(fit@copula)) {
        batch <- copula::rCopula(batch_size, fit@copula)
      } else {
        stop("Unsupported copula type for model '", model_name, "'.")
      }

      # Assign column names and add BatchID
      if (!is.null(self$data)) {
        colnames(batch) <- colnames(self$data)
      }
      batch_dt <- data.table::as.data.table(batch)
      batch_dt[, BatchID := batch_id]
      return(batch_dt)
    },

    # Private method for processing a single batch in conditional_range_predictions()
    process_value_conditional = function(known_values, batch_id) {
      # set.seed(batch_id)
      message(sprintf("Processing conditional prediction for batch_id = %d", batch_id))
      predictions <- self$conditional_prediction(model_name = model_name, known_values = known_values, n = n)

      if (!is.null(predictions)) {
        # Add the conditioned value to the predictions
        return(predictions)
      } else {
        return(NULL)
      }
    },

    # Private method for processing a single batch in hybrid_range_predictions()
    process_value_hybrid = function(value, batch_id) {
      # Set seed for reproducibility
      set.seed(batch_id)

      # Log the process
      message(sprintf("Processing hybrid simulation for batch_id = %d, %f", batch_id, value))

      # Run the hybrid simulation
      sim_output <- self$hybrid_simulation(model_name = model_name, known_values = value, n = n)

      # Validate simulation output
      if (!all(c("Conditional", "Unconditional") %in% names(sim_output))) {
        message("Simulation output must contain 'Conditional' and 'Unconditional' data.tables.")
        return(NULL)
      }

      # Process Conditional data
      Conditional <- sim_output$Conditional
      data.table::setnames(Conditional, names(value), paste0(names(value), "_fixed"))
      vars_cond <- names(Conditional)[!names(Conditional) %in% paste0(names(value), "_fixed")]
      data.table::setnames(Conditional, vars_cond, paste0(vars_cond, "_cond"))

      # Process Unconditional data
      Unconditional <- sim_output$Unconditional
      Unconditional <- Unconditional[, .SD, .SDcols = setdiff(names(Unconditional), names(value))]
      data.table::setnames(Unconditional, paste0(names(Unconditional), "_uncond"))

      # Combine results
      final <- cbind(Conditional, Unconditional)

      # Add batch_id to the result
      final[, batch_id := batch_id]

      return(final)
    },

    # Gaussian Conditional Sampling
    conditional_gaussian = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      dim <- dim(copula_model)

      # Reconstruct the correlation matrix
      corr <- copula::getSigma(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Extract and partition the correlation matrix
      corr_11 <- corr[known_indices, known_indices, drop = FALSE]
      corr_12 <- corr[known_indices, remaining_indices, drop = FALSE]
      corr_22 <- corr[remaining_indices, remaining_indices, drop = FALSE]
      corr_21 <- t(corr_12)

      # Compute the conditional mean and covariance
      inv_corr_11 <- solve(corr_11)
      conditional_mean <- corr_21 %*% inv_corr_11 %*% known_u
      conditional_cov <- corr_22 - corr_21 %*% inv_corr_11 %*% corr_12

      # Generate conditional samples
      conditional_samples <- mvtnorm::rmvnorm(n, mean = as.vector(conditional_mean), sigma = conditional_cov)

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # t Conditional Sampling
    conditional_t = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      dim <- dim(copula_model)
      df <- copula_model@parameters[[2]]  # Degrees of freedom

      # Reconstruct the correlation matrix
      corr <- copula::getSigma(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Extract and partition the correlation matrix
      corr_11 <- corr[known_indices, known_indices, drop = FALSE]
      corr_12 <- corr[known_indices, remaining_indices, drop = FALSE]
      corr_22 <- corr[remaining_indices, remaining_indices, drop = FALSE]
      corr_21 <- t(corr_12)

      # Compute the conditional mean and covariance
      inv_corr_11 <- solve(corr_11)
      conditional_mean <- corr_21 %*% inv_corr_11 %*% known_u
      conditional_cov <- corr_22 - corr_21 %*% inv_corr_11 %*% corr_12

      # Generate conditional samples using mvtnorm::rmvt
      conditional_samples <- mvtnorm::rmvt(
        n = n,
        delta = as.vector(conditional_mean),
        sigma = conditional_cov,
        df = df
      )

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Clayton Conditional Sampling
    conditional_clayton = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Clayton parameter
      dim <- dim(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Compute the conditional distribution
      V <- sum(known_u^(-theta)) - length(known_u) + 1  # Intermediate value
      conditional_samples <- replicate(n, {
        W <- stats::runif(length(remaining_indices))
        (V + W)^(-1 / theta)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- t(conditional_samples)
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Gumbel Conditional Sampling
    conditional_gumbel = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Gumbel parameter
      dim <- dim(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Compute the conditional distribution
      V <- sum((-log(known_u))^(1 / theta))  # Intermediate value
      conditional_samples <- replicate(n, {
        W <- stats::runif(length(remaining_indices))
        exp(-(V - log(W))^theta)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- t(conditional_samples)
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Frank Conditional Sampling
    conditional_frank = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Frank parameter
      dim <- dim(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Compute intermediate value (joint generator evaluation)
      D_inv <- function(x) -log((exp(-theta * x) - 1) / (exp(-theta) - 1))
      V <- -log(1 + (exp(-theta * sum(D_inv(known_u))) - 1) / (exp(-theta) - 1))

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- stats::runif(length(remaining_indices))
        -log(1 + (exp(-theta * V) * (W - 1)) / (W * (exp(-theta) - 1))) / theta
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- t(conditional_samples)
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Joe Conditional Sampling
    conditional_joe = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Joe parameter
      dim <- dim(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Compute intermediate value (generator evaluation for known values)
      V <- sum((1 - known_u)^(-theta)) - length(known_u) + 1  # Joint dependence value

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- stats::runif(length(remaining_indices))
        1 - (V + (1 - W)^(-theta))^(-1 / theta)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- t(conditional_samples)
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Galambos Conditional Sampling
    conditional_galambos = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Galambos parameter
      dim <- dim(copula_model)

      # Match known variable indices
      known_vars <- names(known_values)
      known_indices <- match(known_vars, colnames(self$data))
      remaining_indices <- setdiff(seq_len(dim), known_indices)

      # Convert known_values to pseudo-observations
      known_u <- sapply(known_vars, function(col) {
        value <- known_values[[col]]
        if (!col %in% colnames(self$data)) {
          message(paste0("Column '", col, "' not found in dataset."))
          return(NULL)
        }
        ecdf(self$data[[col]])(value)  # Convert to pseudo-observations
      })

      # Galambos copula generator function
      generator <- function(t) (-log(t))^(1 / theta)

      # Compute conditional dependence for known values
      conditional_cdf <- function(remaining_var, known_u) {
        partial_sum <- sum(generator(known_u))
        (-partial_sum + generator(remaining_var))^(-theta)
      }

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(length(remaining_indices))  # Generate random uniform values
        conditional_cdf(W, known_u)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = dim)
      result[, known_indices] <- matrix(rep(known_u, each = n), nrow = n)
      result[, remaining_indices] <- t(conditional_samples)
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Husler Reiss Conditional Sampling
    conditional_huslerreiss = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      alpha <- copula_model@parameters[[1]]  # Alpha parameter
      if (length(known_values) != 1) stop("Hüsler-Reiss copula is bivariate; only one variable can be conditioned.")

      # Extract known variable and its index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Simulate conditional samples
      conditional_samples <- replicate(n, {
        # Generate a random value for the remaining variable
        W <- runif(1)  # Random uniform value

        # Conditional CDF for the remaining variable
        l2 <- 1 / alpha
        if (remaining_index == 1) {
          # Solve for u1 given u2
          u1 <- W
          u2 <- known_u
          l1 <- log(u1 * u2)
          l3 <- log(u2)
          l32 <- l3 / l1
          l33 <- 1 - l32
          l4 <- log(l32 / (l33))
          l5 <- l2 + 0.5 * alpha * l4
          l6 <- l2 - 0.5 * alpha * l4
        } else {
          # Solve for u2 given u1
          u2 <- W
          u1 <- known_u
          l1 <- log(u1 * u2)
          l3 <- log(u1)
          l32 <- l3 / l1
          l33 <- 1 - l32
          l4 <- log(l32 / (l33))
          l5 <- l2 + 0.5 * alpha * l4
          l6 <- l2 - 0.5 * alpha * l4
        }
        exp(l1 * (l32 * pnorm(l5) + (l33) * pnorm(l6)))
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- as.data.table(result)
      setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # tEV Conditional Sampling
    conditional_tev = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      rho <- copula_model@parameters[[1]]  # Scalar correlation parameter
      df <- copula_model@parameters[[2]]  # Degrees of freedom

      # Match known variable indices
      known_vars <- names(known_values)
      if (length(known_vars) != 1) stop("t-EV copula currently supports conditioning on a single variable.")
      known_index <- match(known_vars, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_vars]])(known_values[[1]])
      known_t <- qt(known_u, df = df)  # Transform to t-distribution quantile

      # Compute conditional mean and variance for the remaining variable
      cond_mean <- rho * known_t
      cond_var <- 1 - rho^2

      # Generate samples from the conditional t-distribution
      conditional_samples <- rt(n, df = df) * sqrt(cond_var) + cond_mean
      conditional_u <- pt(conditional_samples, df = df)  # Transform back to pseudo-observations

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_u
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Plackett Conditional Sampling
    conditional_plackett = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # Plackett association parameter
      if (length(known_values) != 1) {
        message("Plackett copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Function to compute the conditional CDF of the remaining variable
      conditional_cdf <- function(w, u, theta) {
        num <- theta * (w * (1 - u) + u * (1 - w))
        denom <- 1 + (theta - 1) * (w + u - 2 * w * u)
        num / denom
      }

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling
        conditional_cdf(W, known_u, theta)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # FGM Conditional Sampling
    conditional_fgm = function(fit, known_values, n = 1) {
      copula_model <- fit@copula
      theta <- copula_model@parameters[[1]]  # FGM dependence parameter
      if (length(known_values) != 1) {
        message("FGM copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Function to compute the conditional CDF of the remaining variable
      conditional_cdf <- function(w, u, theta) {
        w + theta * u * (1 - u) * (1 - 2 * w)
      }

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling
        conditional_cdf(W, known_u, theta)
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # BB1 Conditional Sampling
    conditional_bb1 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter theta
      delta <- fit$par2  # Tail dependence parameter delta

      if (length(known_values) != 1) {
        message("BB1 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 7, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 7, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # BB6 Conditional Sampling
    conditional_bb2 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter theta
      delta <- fit$par2  # Tail dependence parameter delta

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the BB2 copula.")
        return(NULL)
      }
      if (delta <= 0) {
        message("Delta must be greater than 0 for the BB2 copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("BB2 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 8, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 8, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # BB7 Conditional Sampling
    conditional_bb3 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter theta
      delta <- fit$par2  # Tail dependence parameter delta

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the BB3 copula.")
        return(NULL)
      }

      if (delta <= 1) {
        message("Delta must be greater than 1 for the BB3 copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("BB3 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 9, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 9, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # BB8 Conditional Sampling
    conditional_bb8 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter theta
      delta <- fit$par2  # Tail dependence parameter delta

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the BB8 copula.")
        return(NULL)
      }
      if (delta <= 0) {
        message("Delta must be greater than 0 for the BB8 copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("BB8 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 10, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 10, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated Clayton (180) Conditional Sampling
    conditional_rotated_clayton_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated Clayton (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated Clayton (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated Clayton (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          ((1 - W) * (1 - known_u^(-theta)) + known_u^(-theta))^(-1 / theta)
        } else {
          # u2 is known, solve for u1
          (known_u^(-theta) - (1 - W) * (known_u^(-theta) - 1))^(-1 / theta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated Gumbel (180) Conditional Sampling
    conditional_rotated_gumbel_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated Gumbel (180)

      # Validate parameters
      if (theta <= 1) {
        message("Theta must be greater than 1 for the Rotated Gumbel (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated Gumbel (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 14, par = theta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 14, par = theta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated Joe (180) Conditional Sampling
    conditional_rotated_joe_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated Joe (180)

      # Validate parameters
      if (theta <= 1) {
        message("Theta must be greater than 1 for the Rotated Joe (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated Joe (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 16, par = theta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 16, par = theta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated BB1 (180) Conditional Sampling
    conditional_rotated_bb1_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated BB1 (180)
      delta <- fit$par2  # Additional tail parameter for Rotated BB1 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated BB1 (180) copula.")
        return(NULL)
      }
      if (delta <= 0) {
        message("Delta must be greater than 0 for the Rotated BB1 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated BB1 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 17, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 17, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated BB6 (180) Conditional Sampling
    conditional_rotated_bb6_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated BB6 (180)
      delta <- fit$par2  # Additional tail parameter for Rotated BB6 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated BB6 (180) copula.")
        return(NULL)
      }
      if (delta <= 0) {
        message("Delta must be greater than 0 for the Rotated BB6 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated BB6 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 18, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 18, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated BB7 (180) Conditional Sampling
    conditional_rotated_bb7_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated BB7 (180)
      delta <- fit$par2  # Additional tail parameter for Rotated BB7 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated BB7 (180) copula.")
        return(NULL)
      }
      if (delta <= 1) {
        message("Delta must be greater than 1 for the Rotated BB7 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated BB7 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 19, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 19, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated BB8 (180) Conditional Sampling
    conditional_rotated_bb8_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated BB8 (180)
      delta <- fit$par2  # Additional tail parameter for Rotated BB8 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated BB8 (180) copula.")
        return(NULL)
      }
      if (delta <= 0) {
        message("Delta must be greater than 0 for the Rotated BB8 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated BB8 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 20, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 20, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Tawn Type 1 Conditional Sampling
    conditional_tawn_type1 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Tawn Type 1
      delta <- fit$par2  # Asymmetry parameter for Tawn Type 1

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Tawn Type 1 copula.")
        return(NULL)
      }
      if (delta < 0 || delta > 1) {
        message("Delta must be in [0, 1] for the Tawn Type 1 copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Tawn Type 1 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 104, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 104, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated Tawn Type 1 (180) Conditional Sampling
    conditional_rotated_tawn_type1_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated Tawn Type 1 (180)
      delta <- fit$par2  # Asymmetry parameter for Rotated Tawn Type 1 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated Tawn Type 1 (180) copula.")
        return(NULL)
      }
      if (delta < 0 || delta > 1) {
        message("Delta must be in [0, 1] for the Rotated Tawn Type 1 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated Tawn Type 1 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 114, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 114, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Tawn Type 2 Conditional Sampling
    conditional_tawn_type2 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Tawn Type 2
      delta <- fit$par2  # Asymmetry parameter for Tawn Type 2

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Tawn Type 2 copula.")
        return(NULL)
      }
      if (delta < 0 || delta > 1) {
        message("Delta must be in [0, 1] for the Tawn Type 2 copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Tawn Type 2 copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 204, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 204, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    },

    # Rotated Tawn Type 2 (180) Conditional Sampling
    conditional_rotated_tawn_type2_180 = function(fit, known_values, n = 1) {
      theta <- fit$par  # Dependence parameter for Rotated Tawn Type 2 (180)
      delta <- fit$par2  # Asymmetry parameter for Rotated Tawn Type 2 (180)

      # Validate parameters
      if (theta <= 0) {
        message("Theta must be greater than 0 for the Rotated Tawn Type 2 (180) copula.")
        return(NULL)
      }
      if (delta < 0 || delta > 1) {
        message("Delta must be in [0, 1] for the Rotated Tawn Type 2 (180) copula.")
        return(NULL)
      }

      if (length(known_values) != 1) {
        message("Rotated Tawn Type 2 (180) copula supports conditioning on a single variable.")
        return(NULL)
      }

      # Match known variable index
      known_var <- names(known_values)[1]
      known_index <- match(known_var, colnames(self$data))
      remaining_index <- setdiff(1:2, known_index)

      # Convert known_values to pseudo-observations
      known_u <- ecdf(self$data[[known_var]])(known_values[[1]])

      # Generate conditional samples
      conditional_samples <- replicate(n, {
        W <- runif(1)  # Random uniform value for conditional sampling

        if (known_index == 1) {
          # u1 is known, solve for u2
          VineCopula::BiCopHfunc2(u1 = known_u, u2 = W, family = 214, par = theta, par2 = delta)
        } else {
          # u2 is known, solve for u1
          VineCopula::BiCopHfunc(u1 = W, u2 = known_u, family = 214, par = theta, par2 = delta)
        }
      })

      # Combine known and conditional values
      result <- matrix(NA, nrow = n, ncol = 2)
      result[, known_index] <- rep(known_u, n)
      result[, remaining_index] <- conditional_samples
      result_dt <- data.table::as.data.table(result)
      data.table::setnames(result_dt, colnames(self$data))
      return(result_dt)
    }
  )
)
