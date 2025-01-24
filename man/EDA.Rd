% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eda.R
\name{EDA}
\alias{EDA}
\title{EDA (Exploratory Data Analysis) Class for AutoCopula}
\description{
Provides tools for exploratory data analysis tailored to copula modeling.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{data}}{A \code{data.table} containing the dataset for analysis.}

\item{\code{plots}}{A list of \code{echarts4r} plots generated during the analysis.}

\item{\code{correlation_matrix}}{A \code{data.table} of correlations
Initialize the EDA class}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-EDA-new}{\code{EDA$new()}}
\item \href{#method-EDA-summarize}{\code{EDA$summarize()}}
\item \href{#method-EDA-correlate}{\code{EDA$correlate()}}
\item \href{#method-EDA-visualize_distributions}{\code{EDA$visualize_distributions()}}
\item \href{#method-EDA-visualize_scatterplots}{\code{EDA$visualize_scatterplots()}}
\item \href{#method-EDA-generate_3d_scatter_plot}{\code{EDA$generate_3d_scatter_plot()}}
\item \href{#method-EDA-clone}{\code{EDA$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-new"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$new(data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A \code{data.table} containing the dataset for analysis.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-summarize"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-summarize}{}}}
\subsection{Method \code{summarize()}}{
Calculates mean, median, sd, and the count of missing values for each column.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$summarize()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A \code{data.table} containing the summary statistics.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-correlate"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-correlate}{}}}
\subsection{Method \code{correlate()}}{
Calculates Pearson, Spearman, and Kendall's Tau correlations between all numeric columns in the dataset.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$correlate(input_cols = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{input_cols}}{Names of numeric variables to correlate. If NULL then all numeric columns in the data.table will be utilized}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{data.table} with the pairwise Pearson, Spearman, and Kendall's Tau correlation values for all numeric columns.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-visualize_distributions"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-visualize_distributions}{}}}
\subsection{Method \code{visualize_distributions()}}{
Generates histograms for numeric columns and optionally overlays density lines.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$visualize_distributions(
  input_cols = NULL,
  title_prefix = "Distribution of",
  bins = 20,
  add_density = TRUE,
  tooltip_trigger = "axis",
  theme = "westeros",
  density_opacity = 0.4
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{input_cols}}{Names of numeric variables to plot}

\item{\code{title_prefix}}{Character. Prefix for the plot title.}

\item{\code{bins}}{Integer. Number of bins for the histogram. Defaults to Sturges' formula.}

\item{\code{add_density}}{Logical. Whether to add a density line. Defaults to \code{TRUE}.}

\item{\code{tooltip_trigger}}{"axis"}

\item{\code{theme}}{Character. Theme for the plot}

\item{\code{density_opacity}}{numeric. default 0.4}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of \code{echarts4r} histogram plots.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-visualize_scatterplots"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-visualize_scatterplots}{}}}
\subsection{Method \code{visualize_scatterplots()}}{
Generates scatterplots for the percentile ranks of all numeric column pairs.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$visualize_scatterplots(
  title_prefix = "Empirical Copula View of",
  theme = "westeros"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{title_prefix}}{Character. Prefix for the plot title.}

\item{\code{theme}}{Character. Theme for the plot. Defaults to "westeros".}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of \code{echarts4r} scatter plots for the percentile ranks.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-generate_3d_scatter_plot"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-generate_3d_scatter_plot}{}}}
\subsection{Method \code{generate_3d_scatter_plot()}}{
Generates a 3D scatter plot for three numeric variables.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$generate_3d_scatter_plot(
  col_x,
  col_y,
  col_z,
  rank_values = TRUE,
  theme = "westeros"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{col_x}}{The name of the first numeric column.}

\item{\code{col_y}}{The name of the second numeric column.}

\item{\code{col_z}}{The name of the third numeric column.}

\item{\code{rank_values}}{Logical. Whether to transform variables to their percentile ranks. Defaults to TRUE.}

\item{\code{theme}}{Name of theme for \code{echarts4r} plot}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An \code{echarts4r} 3D scatter plot.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-EDA-clone"></a>}}
\if{latex}{\out{\hypertarget{method-EDA-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EDA$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
