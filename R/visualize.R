#' Visualizing the effect of minimum count on the core size
#'
#' @description The visualize function generates interactive line plots that
#' allow users to explore the impact of different min_count values on the number of
#' core OTUs. Users can interact with the plots to examine the relationship between
#' filtering criteria and core OTU identification visually.
#'
#' @usage visualize(filtered_otu, min_count_val, max_count_val, count_val_interval,
#' prop, min_total_count, method, top_percentage)
#'
#' @param filtered_otu A dataframe of OTUs obtained before filtering which is retrieved from CoreMicrobiome function where the first row is the OTU ID and column names refer to sites/sample names
#' @param min_count_val A numeric value of Minimum count for each OTU to be present in each to be included after the filtering
#' @param max_count_val A numeric value of Maximum count for each OTU to be present in each to be included after the filtering
#' @param count_val_interval Count value interval for each OTU to be present in each to be included after the filtering
#' @param prop Minimum proportion of samples in which an OTU must be present
#' @param min_total_count Minimum total count for each OTU to be included after the filtering
#' @param method Different normalization methods, includes "rrarefy", "srs", "css", "tmm", or "none"
#' @param top_percentage Percentage used for obtaining the Core OTUs
#'
#'
#' @returns
#'  This function gives a line plot which shows change in number of core OTUs with minimum count
#'
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
#' #To view the line plot
#' visualize(filtered_otu = core_1[["final_otu_table_bef_filter"]],
#'          min_count_val = 5,
#'          max_count_val = 25,
#'          count_val_interval = 5,
#'          prop = 0.1,
#'          min_total_count = 10,
#'          method = "srs",
#'          top_percentage =10)

visualize <- function(filtered_otu, min_count_val, max_count_val, count_val_interval,
                       prop, min_total_count, method, top_percentage) {
  type <- NULL

  # min_count_values: Vector of different min_count values you want to test
  min_count_values <- as.vector(seq(min_count_val, max_count_val, by = count_val_interval))  # Adjust the range as needed

  #Getting the core_otu values for the given count values
  core_otu <- function (x) {
    core_otus_values <- numeric(length(min_count_values))

    for (i in seq_along(min_count_values)) {
      abundance_fun1 <- function(x, min_count, prop, min_total_count){
        # Function 'abundance_fun1' logic here
        id1 <- which(x >= min_count)
        pr1 <- length(id1) / length(x)
        sel1 <- ifelse(pr1 >= prop, 1, 0)
        sel2 <- ifelse(sum(x) >= min_total_count, 1, 0)
        sel <- sel1 & sel2
        return(sel)
      }

      sel.ind <- apply(filtered_otu[, -1], 1, abundance_fun1, min_count = min_count_values[i],
                       prop = prop, min_total_count = min_total_count)
      selected_otus <- filtered_otu[sel.ind, ]

      # Check if any OTUs remain after filtering
      if (nrow(selected_otus) == 0) {
        message("All OTUs were filtered out during the '", type, "' step. Please adjust the filtering parameters.")
        return(NULL)  # Return NULL to indicate no OTUs remain after filtering
      }
      # Step 2: normalize_data
      switch(method,
             "rrarefy" = {
               rarefaction_depth <- min(colSums(selected_otus[,-1]))
               rarefied_table <- vegan::rrarefy(selected_otus[,-1], sample = rarefaction_depth)
               rarefied_table_1 <- cbind(selected_otus[, 1, drop = FALSE], rarefied_table)
               colnames(rarefied_table_1) <- colnames(selected_otus)
               selected_otus_1 <- rarefied_table_1
             },
             "srs" = {
               Cmin <- min(colSums(selected_otus[,-1]))
               SRS_output <- SRS::SRS(selected_otus[,-1], Cmin)
               SRS_output_1 <- cbind(selected_otus[, 1, drop = FALSE], SRS_output)
               colnames(SRS_output_1) <- colnames(selected_otus)
               selected_otus_1 <- SRS_output_1
             },
             "css" = {
               cumulative_sums <- colSums(selected_otus[, -1])
               scaling_factors <- cumulative_sums / sum(cumulative_sums)
               css_output <- selected_otus[, -1] * scaling_factors
               css_output_1 <- cbind(selected_otus[, 1, drop = FALSE], css_output)
               colnames(css_output_1) <- colnames(selected_otus)
               selected_otus_1 <- css_output_1
             },
             "tmm" = {
               d <- edgeR::DGEList(counts = selected_otus[,-1])
               d$lib.size <- colSums(d$counts)
               d_tmm <- edgeR::calcNormFactors(d, method = "TMM", refColumn = NULL, geoMeans = d$lib.size)
               tmm_output <- edgeR::cpm(d_tmm, normalized.lib.sizes = TRUE)
               tmm_output_1 <- cbind(selected_otus[, 1, drop = FALSE], tmm_output)
               colnames(tmm_output_1) <- colnames(selected_otus)
               selected_otus_1 <- tmm_output_1
             },
             "none" = {
               # No normalization
               selected_otus_1 <- selected_otus
             },
             "default" = {
               stop("Invalid normalization method. Please choose one of the available methods: rrarefy, srs, css, tmm, none.")
             }
      )

      # Step 3: get_core_noncore_rows
      row_sums <- rowSums(selected_otus_1[, -1])
      otu_summary <- data.frame(OTU_ID = rownames(selected_otus_1), RowSum = row_sums)
      sorted_otu_summary <- otu_summary[order(otu_summary$RowSum, decreasing = TRUE), ]
      num_rows <- ceiling(nrow(sorted_otu_summary) * (top_percentage / 100))
      top_core_rows <- selected_otus_1[sorted_otu_summary$OTU_ID[1:num_rows], ]

      core_otus <- nrow(top_core_rows)

      core_otus_values[i] <- core_otus
    }

    return(core_otus_values)
  }

  core_otu_values <- core_otu(filtered_otu)

  # Create data for plot
  val <- data.frame(min_count_values, core_otu_values)

  # Define a custom theme
  my_theme <- function() {
    ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 14, face = "italic"),
        axis.text = ggplot2::element_text(size = 12),
      )
  }

  #Creating the line plot
  plot <- ggplot2::ggplot(data = val, ggplot2::aes(x = min_count_values, y = core_otu_values, label = round(core_otu_values, 2))) +
    ggplot2::geom_line(color = "green", linewidth = 1.5) +
    ggplot2::geom_point() +
    ggrepel::geom_label_repel() +
    ggplot2::labs(x = "min_count", y = "Number of Core OTUs") +
    ggplot2::ggtitle("Change in Number of Core OTUs with min_count") +
    my_theme()

  return(plot)
}
