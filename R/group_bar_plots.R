#' Grouped Bar Plots Based on Sample Size
#'
#' @description The grouped_bar_plots function is designed for generating grouped
#' bar plots to visualize data. It takes a OTU table before filtering and OTU table
#' after filtering as input containing data for multiple samples and creates a
#' series of grouped bar plots, each representing a specific group of samples.
#'
#' @usage group_bar_plots(otu_table_bef_filtering, otu_table_aft_filtering,
#' num_samples_per_plot)
#'
#' @param otu_table_bef_filtering A data frame of OTUs before filtering where the first row is the OTU ID and column names refer to sites/sample names
#' @param otu_table_aft_filtering A data frame of OTUs after filtering where the first row is the OTU ID and column names refer to sites/sample names
#' @param num_samples_per_plot The number of samples to be displayed in each grouped bar plot.
#'
#' @returns
#'  A list of interactive grouped bar plots, showing the change in sample size before and after filtering OTU table
#'
#' @export
#'
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
#' #To run grouped bar plot function
#' plot_group_bar <- group_bar_plots(core_1$final_otu_table_bef_filter,
#' core_1$final_otu_aft_filter, 10)
#' #To view the grouped bar plot
#' plot_group_bar[[1]]

group_bar_plots <- function(otu_table_bef_filtering, otu_table_aft_filtering, num_samples_per_plot) {
  column_name <- NULL
  column_sum_value <- NULL
  sample_columns_p1 <- colnames(otu_table_bef_filtering[,-1])
  sample_columns_p2 <- colnames(otu_table_aft_filtering[,-1])

  # Calculate the number of plots needed based on the number of samples and the number of samples per plot
  num_plots <- max(ceiling(length(sample_columns_p1) / num_samples_per_plot),
                   ceiling(length(sample_columns_p2) / num_samples_per_plot))

  # Create a list to store the plots
  plots_list <- list()

  # Loop through each group of samples and create the plots
  for (i in 1:num_plots) {
    start_idx_p1 <- (i - 1) * num_samples_per_plot + 1
    end_idx_p1 <- min(i * num_samples_per_plot, length(sample_columns_p1))
    sample_subset_p1 <- sample_columns_p1[start_idx_p1:end_idx_p1]

    start_idx_p2 <- (i - 1) * num_samples_per_plot + 1
    end_idx_p2 <- min(i * num_samples_per_plot, length(sample_columns_p2))
    sample_subset_p2 <- sample_columns_p2[start_idx_p2:end_idx_p2]

    # Calculate column sums (before and after filtering) for both data sets
    column_sums_p1 <- colSums(otu_table_bef_filtering[, sample_subset_p1])
    column_sums_p2 <- colSums(otu_table_aft_filtering[, sample_subset_p2])

    # Convert the column sums into data frames
    column_sum_df_p1 <- data.frame(column_name = names(column_sums_p1),
                                   column_sum_value = column_sums_p1)

    column_sum_df_p2 <- data.frame(column_name = names(column_sums_p2),
                                   column_sum_value = column_sums_p2)

    # Create the plot using ggplot2 for both data sets
    p1 <- ggplot2::ggplot(column_sum_df_p1, ggplot2::aes(x = column_name, y = column_sum_value, fill = column_name)) +
      ggplot2::geom_bar(stat = 'identity', position = 'dodge') +
      ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0.5),
                     axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = "Before filtering", y = "Sample size")

    p2 <- ggplot2::ggplot(column_sum_df_p2, ggplot2::aes(x = column_name, y = column_sum_value, fill = column_name)) +
      ggplot2::geom_bar(stat = 'identity', position = 'dodge') +
      ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0.5),
                     axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = "After filtering", y = "Sample size", title = "Grouped Bar Plot")

    # Convert ggplot to plotly for both data sets
    p1 <- plotly::ggplotly(p1)
    p2 <- plotly::ggplotly(p2)

    # Arrange the plots side by side with a common legend using subplot
    combined_plot_side_by_side <- plotly::subplot(p1, p2, shareY = TRUE, shareX = TRUE, nrows = 1)

    # Add the plot to the list
    plots_list[[i]] <- combined_plot_side_by_side
  }

  # Return the list of plots
  return(plots_list)
}
