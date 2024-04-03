#' Stacked Bar Plots Based on Relative Abundance Data
#'
#' @description This function generates stacked bar plots for visualizing
#' the relative abundance data of different operational taxonomic units (OTUs)
#' in various samples.
#'
#' @usage stacked_bar_plots(data, num_samples_per_plot)
#'
#' @param data A data frame containing the relative abundance data for the OTUs. The first column should contain the OTU IDs, and the subsequent columns should represent samples.
#' @param num_samples_per_plot The number of samples to be displayed in each stacked bar plot.
#'
#' @returns
#'  A list of interactive stacked bar plots, one for each group of samples, showing the relative abundance of OTUs in the samples.
#'
#' @export
#' @examples
#' #To run input data
#' core_1 <- CoreMicrobiome(
#'  otu_table = demo_otu,
#'  tax_table = demo_tax,
#'  metadata_table = demo_md,
#'  filter_type = "occupancy_fun_filter", #Or "abundance_fun_filter", Or "combined_filter"
#'  percent = 0.5,
#'  method = "css",  # Or "srs", "rrarefy", "tmm", "tmmwsp", "rle", "upperquartile", "none"
#'  beta_diversity_method = "jaccard",
#'  top_percentage = 10  # Adjust the percentage as needed for core/non-core OTUs
#')
#' #To run the stacked bar plots function
#' stacked_plots <- stacked_bar_plots(core_1$core_otus_relative_abundance, 10)
#' #To view the stacked bar plot
#' stacked_plots[[1]]


stacked_bar_plots <- function(data, num_samples_per_plot) {
  variable <- NULL
  Relative_Abundance <- NULL
  `OTU ID` <- NULL
  sample_columns <- colnames(data)

  # Calculate the number of plots needed based on the number of samples and the number of samples per plot
  num_plots <- ceiling(length(sample_columns) / num_samples_per_plot)

  # Create a list to store the plots
  plots_list <- list()

  # Loop through each group of samples and create the plots
  for (i in 1:num_plots) {
    start_idx <- (i - 1) * num_samples_per_plot + 1
    end_idx <- min(i * num_samples_per_plot, length(sample_columns))

    sample_subset <- sample_columns[start_idx:end_idx]

    plot_data <- reshape2::melt(data[, c(1, start_idx:end_idx)], varnames = c("OTU_ID", "Sample"), value.name = "Relative_Abundance")

    plot <- plotly::ggplotly(ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = Relative_Abundance, fill = `OTU ID`)) +
                               ggplot2::geom_bar(stat = "identity") +
                               ggplot2::labs(title = paste("Stacked Bar Plot for Samples", start_idx, "to", end_idx),
                                             x = "Samples",
                                             y = "Relative Abundance",
                                             fill = "OTU ID") +
                               ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))) # Optional: Rotate x-axis labels for better readability

    # Add the plot to the list
    plots_list[[i]] <- plot
  }

  # Return the list of plots
  return(plots_list)
}
