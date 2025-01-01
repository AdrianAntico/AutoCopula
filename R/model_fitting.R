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
      BB2 = list(
        description = "BB2 copula for flexible tail dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB2 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 8)
        }
      ),
      BB3 = list(
        description = "BB3 copula for asymmetrical dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB3 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 9)
        }
      ),
      BB4 = list(
        description = "BB4 copula for modeling both tail dependencies.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("BB4 copula requires exactly 2 dimensions.")
          VineCopula::BiCopEst(data[, 1], data[, 2], family = 10)
        }
      ),

      # Additional copulas from the evd package
      Logistic = list(
        description = "Logistic copula for symmetric extreme value dependence.",
        fit_function = function(data) {
          if (ncol(data) != 2) stop("Logistic copula requires exactly 2 dimensions.")
          evd::fbvevd(data, model = "log")
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
