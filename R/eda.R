#' @title EDA (Exploratory Data Analysis) Class for AutoCopula
#'
#' @description Provides tools for exploratory data analysis tailored to copula modeling.
#'
#' @export
EDA <- R6::R6Class(
  "EDA",
  public = list(
    #' @field data A `data.table` containing the dataset for analysis.
    data = NULL,
    #' @field plots A list of `echarts4r` plots generated during the analysis.
    plots = list(),
    #' @field correlation_matrix A `data.table` of correlations
    correlation_matrix = NULL, # Initialize as NULL

    #' Initialize the EDA class
    #'
    #' @param data A `data.table` containing the dataset for analysis.
    initialize = function(data) {
      if (!"data.table" %in% class(data)) {
        stop("Input data must be a data.table object.")
      }
      self$data <- data
    },

    #' @description Calculates mean, median, sd, and the count of missing values for each column.
    #'
    #' @return A `data.table` containing the summary statistics.
    #' @export
    summarize = function() {
      # Process numeric columns
      numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
      numeric_summary <- if (length(numeric_cols) > 0) {
        lapply(numeric_cols, function(col_name) {
          col <- self$data[[col_name]]
          data.table::data.table(
            Variable = col_name,
            Mean = mean(col, na.rm = TRUE),
            Median = median(col, na.rm = TRUE),
            StDev = sd(col, na.rm = TRUE),
            NA_Count = sum(is.na(col))
          )
        })
      } else {
        list()
      }

      # Process categorical columns
      categorical_cols <- names(self$data)[sapply(self$data, function(col) is.factor(col) || is.character(col))]
      categorical_summary <- if (length(categorical_cols) > 0) {
        lapply(categorical_cols, function(col_name) {
          col <- self$data[[col_name]]
          data.table::data.table(
            Variable = col_name,
            Mean = NA_real_,
            Median = NA_real_,
            StDev = NA_real_,
            NA_Count = sum(is.na(col))
          )
        })
      } else {
        list()
      }

      # Combine numeric and categorical summaries
      self$summary_stats <- rbindlist(c(numeric_summary, categorical_summary), fill = TRUE)

      return(self$summary_stats)
    },

    #' @description Calculates Pearson, Spearman, and Kendall's Tau correlations between all numeric columns in the dataset.
    #'
    #' @param input_cols Names of numeric variables to correlate. If NULL then all numeric columns in the data.table will be utilized
    #' @return A `data.table` with the pairwise Pearson, Spearman, and Kendall's Tau correlation values for all numeric columns.
    #' @export
    correlate = function(input_cols = NULL) {

      # Identify numeric columns excluding the target column
      if(is.null(input_cols)) {
        numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
        if (length(numeric_cols) < 2) {
          stop("Not enough columns are numeric")
        }
      } else {
        numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
        numeric_cols <- numeric_cols[numeric_cols %in% input_cols]
        if (length(numeric_cols) < 2) {
          stop("Not enough columns are numeric")
        }
      }

      # Initialize an empty list to collect results
      correlation_results <- list()

      # Compute correlations for all numeric pairs
      for (i in seq_along(numeric_cols)) {
        for (j in seq_along(numeric_cols)) {
          if (i < j) {
            col_x <- numeric_cols[i]
            col_y <- numeric_cols[j]
            pearson_corr <- tryCatch(
              stats::cor(self$data[[col_x]], self$data[[col_y]], use = "complete.obs", method = "pearson"),
              error = function(e) NA
            )
            spearman_corr <- tryCatch(
              stats::cor(self$data[[col_x]], self$data[[col_y]], use = "complete.obs", method = "spearman"),
              error = function(e) NA
            )
            kendall_corr <- tryCatch(
              stats::cor(self$data[[col_x]], self$data[[col_y]], use = "complete.obs", method = "kendall"),
              error = function(e) NA
            )
            # Append results
            correlation_results <- append(correlation_results, list(
              data.table::data.table(
                Variable_X = col_x,
                Variable_Y = col_y,
                Pearson = pearson_corr,
                Spearman = spearman_corr,
                Kendall = kendall_corr,
                Difference_Pearson_Spearman = pearson_corr - spearman_corr
              )
            ))
          }
        }
      }

      # Combine results into a single data.table
      if (length(correlation_results) > 0) {
        self$correlation_matrix <- data.table::rbindlist(correlation_results, fill = TRUE)
      } else {
        self$correlation_matrix <- data.table::data.table(
          Variable_X = character(),
          Variable_Y = character(),
          Pearson = numeric(),
          Spearman = numeric(),
          Kendall = numeric(),
          Difference_Pearson_Spearman = numeric()
        )
      }

      return(self$correlation_matrix)
    },

    #' @description Generates histograms for numeric columns and optionally overlays density lines.
    #'
    #' @param input_cols Names of numeric variables to plot
    #' @param title_prefix Character. Prefix for the plot title.
    #' @param bins Integer. Number of bins for the histogram. Defaults to Sturges' formula.
    #' @param add_density Logical. Whether to add a density line. Defaults to `TRUE`.
    #' @param tooltip_trigger "axis"
    #' @param theme Character. Theme for the plot
    #' @param density_opacity numeric. default 0.4
    #' @return A list of `echarts4r` histogram plots.
    #' @export
    visualize_distributions = function(
    input_cols = NULL,
    title_prefix = "Distribution of",
    bins = 20,
    add_density = TRUE,
    tooltip_trigger = "axis",
    theme = "westeros",
    density_opacity = 0.4) {

      # Clear self$plots to avoid mixing states
      self$plots <- list()

      # Identify numeric columns excluding the target column
      if(is.null(input_cols)) {
        numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
        if (length(numeric_cols) == 0) {
          stop("No columns are numeric")
        }
      } else {
        numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
        numeric_cols <- numeric_cols[numeric_cols %in% input_cols]
        if (length(numeric_cols) == 0) {
          stop("No input_cols are numeric")
        }
      }

      # Validate numeric columns
      if (length(numeric_cols) == 0) {
        return(list())
      }

      for (col in numeric_cols) {
        # Prepare the dataset for plotting
        plot_data <- self$data[, .(Value = get(col))]

        # Validate the column data
        if (nrow(plot_data) == 0 || !is.numeric(plot_data$Value)) {
          stop(paste("Column", col, "is not numeric or contains no data."))
        }

        # Create histogram with optional density overlay
        plot <- plot_data |>
          echarts4r::e_charts() |> # Initialize the plot
          echarts4r::e_histogram(Value, name = "Histogram", breaks = bins) |>
          echarts4r::e_title(text = paste(title_prefix, col)) |>
          echarts4r::e_tooltip(trigger = tooltip_trigger) |>
          echarts4r::e_theme(theme) |>
          echarts4r::e_x_axis(name = col) |>
          echarts4r::e_y_axis(name = "Hist") |>
          echarts4r::e_legend(show = TRUE, type = "scroll", orient = "horizontal", right = 50, top = 30) |>
          echarts4r::e_datazoom(x_index = c(0,1)) |>
          echarts4r::e_toolbox_feature(feature = c("saveAsImage","dataZoom"))

        if (add_density) {
          plot <- plot |>
            echarts4r::e_density(
              Value,
              areaStyle = list(opacity = density_opacity),
              smooth = TRUE,
              name = "Density",
              y_index = 1) |>
            echarts4r::e_theme(theme)
        }

        # Save the plot
        self$plots[[col]] <- plot
      }
      return(self$plots)
    },

    #' @description Generates scatterplots for the percentile ranks of all numeric column pairs.
    #'
    #' @param title_prefix Character. Prefix for the plot title.
    #' @param theme Character. Theme for the plot. Defaults to "westeros".
    #' @return A list of `echarts4r` scatter plots for the percentile ranks.
    #' @export
    visualize_scatterplots = function(
    title_prefix = "Empirical Copula View of",
    theme = "westeros") {
      # Identify numeric columns
      numeric_cols <- names(self$data)[sapply(self$data, is.numeric)]
      if (length(numeric_cols) < 2) {
        stop("Not enough numeric columns to create scatterplots.")
      }

      # Reset plots
      self$plots <- list()

      # Compute percentile ranks
      rank_data <- self$data[, lapply(.SD, function(x) {
        data.table::frank(x, na.last = "keep") / sum(!is.na(x))
      }), .SDcols = numeric_cols]

      # Generate scatterplots for all pairs
      for (i in seq_along(numeric_cols)) {
        for (j in seq_along(numeric_cols)) {
          if (i < j) {
            col_x <- numeric_cols[i]
            col_y <- numeric_cols[j]

            # Prepare data for plotting
            plot_data <- data.table::data.table(
              X = rank_data[[col_x]],
              Y = rank_data[[col_y]]
            )

            # Generate scatterplot
            plot <- plot_data |>
              echarts4r::e_charts(X) |>
              echarts4r::e_scatter(Y, name = "Empirical Copula View") |>
              echarts4r::e_title(text = paste(title_prefix, col_x, "and", col_y)) |>
              echarts4r::e_x_axis(name = paste("Percentile Rank of", col_x), min = 0, max = 1) |>
              echarts4r::e_y_axis(name = paste("Percentile Rank of", col_y), min = 0, max = 1) |>
              echarts4r::e_tooltip(trigger = "axis") |>
              echarts4r::e_theme(name = theme) |>
              echarts4r::e_legend(show = TRUE, type = "scroll", orient = "horizontal", right = 50, top = 30) |>
              echarts4r::e_datazoom(x_index = c(0, 1)) |>
              echarts4r::e_toolbox_feature(feature = c("saveAsImage", "dataZoom"))

            # Save the plot
            plot_key <- paste(col_x, col_y, sep = "_vs_")
            self$plots[[plot_key]] <- plot
          }
        }
      }

      return(self$plots)
    },

    #' @description Generates a 3D scatter plot for three numeric variables.
    #'
    #' @param col_x The name of the first numeric column.
    #' @param col_y The name of the second numeric column.
    #' @param col_z The name of the third numeric column.
    #' @param rank_values Logical. Whether to transform variables to their percentile ranks. Defaults to TRUE.
    #' @return An `echarts4r` 3D scatter plot.
    #' @export
    generate_3d_scatter_plot = function(
    col_x,
    col_y,
    col_z,
    rank_values = TRUE,
    theme = "westeros") {
      if (!(col_x %in% names(self$data) && col_y %in% names(self$data) && col_z %in% names(self$data))) {
        stop("Columns not found in the dataset.")
      }
      if (!is.numeric(self$data[[col_x]]) || !is.numeric(self$data[[col_y]]) || !is.numeric(self$data[[col_z]])) {
        stop("All specified columns must be numeric.")
      }

      # Optionally rank variables
      if (rank_values) {
        plot_data <- self$data[, .(
          X = frank(get(col_x), na.last = "keep") / sum(!is.na(get(col_x))),
          Y = frank(get(col_y), na.last = "keep") / sum(!is.na(get(col_y))),
          Z = frank(get(col_z), na.last = "keep") / sum(!is.na(get(col_z)))
        )]
      } else {
        plot_data <- self$data[, .(
          X = get(col_x),
          Y = get(col_y),
          Z = get(col_z)
        )]
      }

      # Generate 3D scatter plot
      plot <- plot_data |>
        echarts4r::e_charts(X) |>
        echarts4r::e_scatter_3d(Y, Z, symbol_size = 5) |>
        echarts4r::e_title(text = paste("3D Scatter Plot of", col_x, col_y, "and", col_z)) |>
        echarts4r::e_x_axis_3d(name = col_x) |>
        echarts4r::e_y_axis_3d(name = col_y) |>
        echarts4r::e_z_axis_3d(name = col_z) |>
        echarts4r::e_theme(theme)

      self$plots[[paste(col_x, col_y, col_z, sep = "_vs_")]] <- plot
      return(plot)
    }
  )
)
