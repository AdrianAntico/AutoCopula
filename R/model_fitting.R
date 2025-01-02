#' ModelFitter
#'
#' An R6 class for fitting and evaluating copula models.
#' Provides a library of pre-defined copula models and tools for custom models.
#'
#' @export
ModelFitter <- R6::R6Class(
  "ModelFitter",
  public = list(
    #' @field data A data.table containing the dataset for modeling.
    data = NULL,

    #' @field copula_library A pre-defined library of copula models.
    copula_library = list(

      # Multivariate Copulas from copula
      Gaussian = list(
        description = "Gaussian copula with a correlation matrix.",
        fit_function = function(data) copula::fitCopula(copula::normalCopula(dim = ncol(data)), data, method = "ml")
      ),
      tCopula = list(
        description = "t-copula with degrees of freedom and a correlation matrix.",
        fit_function = function(data) copula::fitCopula(copula::tCopula(dim = ncol(data)), data, method = "ml")
      ),
      Clayton = list(
        description = "Clayton copula for lower tail dependence.",
        fit_function = function(data) copula::fitCopula(copula::claytonCopula(dim = ncol(data)), data, method = "ml")
      ),
      Gumbel = list(
        description = "Gumbel copula for upper tail dependence.",
        fit_function = function(data) copula::fitCopula(copula::gumbelCopula(dim = ncol(data)), data, method = "ml")
      ),
      Frank = list(
        description = "Frank copula with no tail dependence.",
        fit_function = function(data) copula::fitCopula(copula::frankCopula(dim = ncol(data)), data, method = "ml")
      ),
      Joe = list(
        description = "Joe copula with upper tail dependence.",
        fit_function = function(data) copula::fitCopula(copula::joeCopula(dim = ncol(data)), data, method = "ml")
      ),

      # Bivariate Copulas from copula
      Galambos = list(
        description = "Galambos copula for modeling extreme values.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Galambos copula requires exactly 2 dimensions.")
          copula::fitCopula(copula::galambosCopula(), data, method = "ml")
        }
      ),
      HuslerReiss = list(
        description = "Husler-Reiss copula for extreme value modeling.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Husler-Reiss copula requires exactly 2 dimensions.")
          copula::fitCopula(copula::huslerReissCopula(), data, method = "ml")
        }
      ),
      tEV = list(
        description = "t-EV copula for extreme value dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("t-EV copula requires exactly 2 dimensions.")
          copula::fitCopula(copula::tevCopula(), data, method = "ml")
        }
      ),
      Plackett = list(
        description = "Plackett copula for symmetric dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Plackett copula requires exactly 2 dimensions.")
          copula::fitCopula(copula::plackettCopula(), data, method = "ml")
        }
      ),
      FGM = list(
        description = "FGM (Farlie-Gumbel-Morgenstern) copula for weak dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("FGM copula requires exactly 2 dimensions.")
          copula::fitCopula(copula::fgmCopula(), data, method = "ml")
        }
      ),

      # Additional Copulas from VineCopula
      BB1 = list(
        description = "BB1 copula for flexible upper and lower tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB1 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 7)
        }
      ),
      BB6 = list(
        description = "BB6 copula for flexible tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB6 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 8)
        }
      ),
      BB7 = list(
        description = "BB7 copula for asymmetrical dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB7 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 9)
        }
      ),
      BB8 = list(
        description = "BB8 copula for modeling both tail dependencies.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB8 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 10)
        }
      ),

      # Rotated Copulas (180-degree)
      "Rotated Clayton (180)" = list(
        description = "180-degree rotated Clayton copula (survival Clayton) for upper tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated Clayton (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 13)
        }
      ),

      "Rotated Gumbel (180)" = list(
        description = "180-degree rotated Gumbel copula (survival Gumbel) for lower tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated Gumbel (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 14)
        }
      ),

      "Rotated Joe (180)" = list(
        description = "180-degree rotated Joe copula (survival Joe) for upper tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated Joe (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 16)
        }
      ),

      "Rotated BB1 (180)" = list(
        description = "180-degree rotated BB1 copula for upper and lower tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated BB1 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 17)
        }
      ),

      "Rotated BB6 (180)" = list(
        description = "180-degree rotated BB6 copula for tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated BB6 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 18)
        }
      ),

      "Rotated BB7 (180)" = list(
        description = "180-degree rotated BB7 copula for asymmetrical dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated BB7 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 19)
        }
      ),

      "Rotated BB8 (180)" = list(
        description = "180-degree rotated BB8 copula for modeling tail dependencies.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated BB8 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 20)
        }
      ),

      # Tawn Copulas
      "Tawn Type 1" = list(
        description = "Tawn Type 1 copula for asymmetrical lower tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Tawn Type 1 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 104)
        }
      ),

      "Rotated Tawn Type 1 (180)" = list(
        description = "180-degree rotated Tawn Type 1 copula for asymmetrical upper tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated Tawn Type 1 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 114)
        }
      ),

      "Tawn Type 2" = list(
        description = "Tawn Type 2 copula for asymmetrical upper tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Tawn Type 2 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 204)
        }
      ),

      "Rotated Tawn Type 2 (180)" = list(
        description = "180-degree rotated Tawn Type 2 copula for asymmetrical lower tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Rotated Tawn Type 2 (180) copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 214)
        }
      )
    ),

    #' @field fit_results A list to store results of copula fits.
    fit_results = list(),

    #' Initialize the ModelFitter class
    #'
    #' @param data A data.table containing the dataset for modeling.
    #' @return A new instance of the ModelFitter class.
    initialize = function(data) {
      if (!data.table::is.data.table(data)) stop("Input data must be a data.table.")
      self$data <- data
    },

    #' @description Fits multiple copula models in a loop.
    #'
    #' @param model_names A character vector of model names to fit.
    #' @return A data.table summarizing the results for each model.
    #' @export
    fit_models = function(model_names) {
      if (length(model_names) == 0) stop("You must specify at least one model name.")

      lapply(model_names, function(model_name) {
        tryCatch(
          {
            # Fit the model using the private method
            fit <- private$fit_model(model_name)

            # Store the result in the fit_results list
            self$fit_results[[model_name]] <- fit
          },
          error = function(e) {
            message(sprintf("Error fitting model '%s': %s", model_name, e$message))
          }
        )
      })

      # No return value; results are stored in self$fit_results
      invisible(NULL)
    },

    #' @description Lists available copula models in the library.
    #'
    #' @return A data.table summarizing the copula models.
    #' @export
    list_models = function() {
      data.table::data.table(
        Model = names(self$copula_library),
        Description = sapply(self$copula_library, function(x) x$description)
      )
    }
  ),

  private = list(
    fit_model = function(model_name) {
      # Validate that the model exists in the copula library
      if (!model_name %in% names(self$copula_library)) {
        stop("Model not found in library. Use list_models() to see available models.")
      }

      # Transform data to pseudo-observations
      data_transformed <- copula::pobs(as.matrix(self$data))

      # Retrieve model fitting function from the copula library
      model_info <- self$copula_library[[model_name]]

      # Fit the model using the specified fitting function
      fit <- tryCatch({
        model_info$fit_function(data_transformed)
      }, error = function(e) {
        message("Error fitting model '", model_name, "': ", e$message)
        NULL
      })

      # Return the fitted model (or NULL if fitting failed)
      if (!is.null(fit)) {
        message("Successfully fitted model: ", model_name)
      }
      return(fit)
    }
  )
)
