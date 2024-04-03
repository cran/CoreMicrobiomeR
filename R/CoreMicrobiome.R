CoreMicrobiome <- function(otu_table, tax_table, metadata_table, filter_type, ...,
                           method, beta_diversity_method, top_percentage) {
  type <- NULL
  # Step 1: sorted_otu_md_func
  colnames(otu_table)[1] <- "OTU ID"
  colnames(tax_table)[1] <- "OTU ID"
  id.common <- fastmatch::fmatch(tax_table[, 1][[1]], otu_table[, 1][[1]])
  common_otu_ids_table1 <- otu_table[id.common, ]
  common_otu_ids_table <- common_otu_ids_table1[, -1]
  common_otu_ids_table <- as.data.frame(common_otu_ids_table)
  arab_md <- as.data.frame(metadata_table[, 1])
  colnames(arab_md) <- "Sample"
  id.com.md <- fastmatch::fmatch(arab_md[, 1], colnames(common_otu_ids_table))
  id.col <- stats::na.omit(id.com.md)
  final_otu_table <- common_otu_ids_table[, id.col]
  final_otu_table_with_row <- cbind(common_otu_ids_table1[, 1], final_otu_table)
  colnames(final_otu_table_with_row)[1] <- "OTU ID"
  # Add an "id" column to the tax_table using the OTU IDs
  #tax_table$id <- rownames(tax_table)
  id.md <- which(!is.na(id.com.md))
  filtered_arab_md_old <- metadata_table[id.md, ]

  # Step 2: filtered_non_zero_otu
  row_sums <- rowSums(final_otu_table_with_row[,-1])
  non_zero <- row_sums != 0
  filtered_otu_tab <- final_otu_table_with_row[non_zero, ]

  # Step 3: combined_filter_function
  selected_otus <- NULL
  if (filter_type == "abundance_fun_filter") {
    # Check if the required arguments are provided
    if (!("min_count" %in% names(list(...))) ||
        !("prop" %in% names(list(...))) ||
        !("min_total_count" %in% names(list(...)))) {
      stop("For 'abundance_fun_filter', please provide values for 'min_count', 'prop', and 'min_total_count' parameters.")
    }

    abundance_fun1 <- function(x, min_count, prop, min_total_count) {
      # Function 'abundance_fun1' logic here
      id1 <- which(x >= min_count)
      pr1 <- length(id1) / length(x)
      sel1 <- ifelse(pr1 >= prop, 1, 0)
      sel2 <- ifelse(sum(x) >= min_total_count, 1, 0)
      sel <- sel1 & sel2
      return(sel)
    }
    sel.ind <- apply(filtered_otu_tab[, -1], 1, abundance_fun1,
                     min_count = list(...)$min_count,
                     prop = list(...)$prop,
                     min_total_count = list(...)$min_total_count)
    selected_otus <- filtered_otu_tab[sel.ind, ]
  } else if (filter_type == "occupancy_fun_filter") {
    # Check if the required arguments are provided
    if (!("percent" %in% names(list(...)))){
      stop("For 'occupancy_fun_filter', please provide a value for the 'percent' parameter.")
    }

    nc <- ncol(filtered_otu_tab) - 1
    nonzero_counts <- rowSums(filtered_otu_tab[, -1] != 0)
    filtered_ids <- rownames(filtered_otu_tab)[nonzero_counts >= (list(...)$percent * nc)]
    selected_otus <- filtered_otu_tab[filtered_ids, ]
  } else if (filter_type == "combined_filter") {
    # Check if the required arguments are provided
    if (!("percent" %in% names(list(...))) ||
        !("min_count" %in% names(list(...))) ||
        !("prop" %in% names(list(...))) ||
        !("min_total_count" %in% names(list(...)))) {
      stop("For 'combined_filter', please provide a value for the 'percent', 'min_count', 'prop', and 'min_total_count' parameter.")
    }

    # Apply occupancy filter
    nc <- ncol(filtered_otu_tab) - 1
    nonzero_counts <- rowSums(filtered_otu_tab[, -1] != 0)
    filtered_ids <- rownames(filtered_otu_tab)[nonzero_counts >= (list(...)$percent * nc)]
    selected_otus_occupancy <- filtered_otu_tab[filtered_ids, ]

    # Apply abundance filter on the results of occupancy filter
    abundance_fun1 <- function(x, min_count, prop, min_total_count) {
      # Function 'abundance_fun1' logic here
      id1 <- which(x >= min_count)
      pr1 <- length(id1) / length(x)
      sel1 <- ifelse(pr1 >= prop, 1, 0)
      sel2 <- ifelse(sum(x) >= min_total_count, 1, 0)
      sel <- sel1 & sel2
      return(sel)
    }
    sel.ind_abundance <- apply(selected_otus_occupancy[, -1], 1, abundance_fun1,
                               min_count = list(...)$min_count,
                               prop = list(...)$prop,
                               min_total_count = list(...)$min_total_count)
    selected_otus_abundance <- selected_otus_occupancy[sel.ind_abundance, ]

    selected_otus <- selected_otus_abundance
  } else {
    stop("Invalid 'filter_type' argument. Please provide 'abundance_fun_filter' or 'occupancy_fun_filter' or 'combined_filter'.")
  }

  # Check if any OTUs remain after filtering
  if (nrow(selected_otus) == 0) {
    message("All OTUs were filtered out during the '", type, "' step. Please adjust the filtering parameters.")
    return(NULL)  # Return NULL to indicate no OTUs remain after filtering
  }

  # Step 4: normalize_data
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
         "tmmwsp" = {
           d <- edgeR::DGEList(counts = selected_otus[,-1])
           d$lib.size <- colSums(d$counts)
           d_tmmwsp <- edgeR::calcNormFactors(d, method = "TMMwsp", refColumn = NULL, geoMeans = d$lib.size)
           tmmwsp_output <- edgeR::cpm(d_tmmwsp, normalized.lib.sizes = TRUE)
           tmmwsp_output_1 <- cbind(selected_otus[, 1, drop = FALSE], tmmwsp_output)
           colnames(tmmwsp_output_1) <- colnames(selected_otus)
           selected_otus_1 <- tmmwsp_output_1
         },
         "rle" = {
           d <- edgeR::DGEList(counts = selected_otus[,-1])
           d$lib.size <- colSums(d$counts)
           d_rle <- edgeR::calcNormFactors(d, method = "RLE", refColumn = NULL, geoMeans = d$lib.size)
           rle_output <- edgeR::cpm(d_rle, normalized.lib.sizes = TRUE)
           rle_output_1 <- cbind(selected_otus[, 1, drop = FALSE], rle_output)
           colnames(rle_output_1) <- colnames(selected_otus)
           selected_otus_1 <- rle_output_1
         },
         "upperquartile" = {
           d <- edgeR::DGEList(counts = selected_otus[,-1])
           d$lib.size <- colSums(d$counts)
           d_upperquartile <- edgeR::calcNormFactors(d, method = "upperquartile", refColumn = NULL, geoMeans = d$lib.size)
           upperquartile_output <- edgeR::cpm(d_upperquartile, normalized.lib.sizes = TRUE)
           upperquartile_output_1 <- cbind(selected_otus[, 1, drop = FALSE], upperquartile_output)
           colnames(upperquartile_output_1) <- colnames(selected_otus)
           selected_otus_1 <- upperquartile_output_1
         },
         "none" = {
           # No normalization
           selected_otus_1 <- selected_otus
         },
         "default" = {
           d <- edgeR::DGEList(counts = selected_otus[,-1])
           d$lib.size <- colSums(d$counts)
           d_tmm <- edgeR::calcNormFactors(d, method = "TMM", refColumn = NULL, geoMeans = d$lib.size)
           tmm_output <- edgeR::cpm(d_tmm, normalized.lib.sizes = TRUE)
           tmm_output_1 <- cbind(selected_otus[, 1, drop = FALSE], tmm_output)
           colnames(tmm_output_1) <- colnames(selected_otus)
           selected_otus_1 <- tmm_output_1
         }
  )

  # Step 5: calculate_diversity
  datafromstepstep4 <- selected_otus_1
  alpha_diversity <- NULL
  bd <- NULL
  tryCatch({
    # Calculate the alpha diversity measures
    richness <- vegan::specnumber(t(datafromstepstep4[, -1]))
    evenness <- vegan::diversity(t(datafromstepstep4[, -1])) / log(richness)
    shannon <- vegan::diversity(t(datafromstepstep4[, -1]), index = "shannon")
    simpson <- vegan::diversity(t(datafromstepstep4[, -1]), index = "simpson")

    # Create a matrix to store the alpha diversity results
    alpha_diversity <- rbind(richness, evenness, shannon, simpson)

    # Choose beta diversity method (default: "bray")
    beta_diversity_method <- "bray"  # Default method
    if ("beta_diversity_method" %in% names(list(...))) {
      beta_diversity_method <- list(...)$beta_diversity_method
    }

    # Calculate chosen beta diversity matrix
    if (beta_diversity_method == "bray") {
      bd <- vegan::vegdist(t(datafromstepstep4[, -1]), method = "bray", diag = TRUE)
    } else if (beta_diversity_method == "jaccard") {
      bd <- vegan::vegdist(t(datafromstepstep4[, -1]), method = "jaccard", diag = TRUE)
    } else if (beta_diversity_method == "mountford") {
      bd <- vegan::vegdist(t(datafromstepstep4[, -1]), method = "mountford", diag = TRUE)
    } else {
      stop("Invalid beta diversity method. Please choose 'bray', 'jaccard', or 'mountford'.")
    }
  }, error = function(e) {
    # Handle error (if any) and print a message
    message("Error encountered during diversity analysis:", e$message)
    alpha_diversity <- NULL
    bd <- NULL
  }, warning = function(w) {
    # Handle warning (if any) and print a message
    message("Warning encountered during diversity analysis:", w$message)
    alpha_diversity <- NULL
    bd <- NULL
  })

  # Step 6: get_core_noncore_rows
  row_sums <- rowSums(selected_otus_1[, -1])
  otu_summary <- data.frame(OTU_ID = rownames(selected_otus_1), RowSum = row_sums)
  sorted_otu_summary <- otu_summary[order(otu_summary$RowSum, decreasing = TRUE), ]
  num_rows <- ceiling(nrow(sorted_otu_summary) * (top_percentage / 100))
  top_core_rows <- selected_otus_1[sorted_otu_summary$OTU_ID[1:num_rows], ]
  non_core_rows <- selected_otus_1[sorted_otu_summary$OTU_ID[(num_rows + 1):nrow(sorted_otu_summary)], ]

  # Step 7: taxonomy of core otus

  # Find common rows based on 'OTU ID' column
  common_otu <- fastmatch::fmatch(top_core_rows[, 1], final_otu_table_with_row[, 1])
  common_otu_tax <- tax_table[common_otu, ]

  #step 8 : core otus original count data

  match_core <- fastmatch::fmatch(top_core_rows[, 1], selected_otus[, 1])
  core_otus_demo_count <- selected_otus[match_core, ]

  #relative abundance function

  # Step 1: Calculate the total counts for each sample (sum across all OTUs)
  sample_total_counts <- colSums(core_otus_demo_count[,-1])

  # Step 2: Calculate the relative abundance for each OTU in each sample
  otu_relative_abundance <- core_otus_demo_count[,-1] / sample_total_counts

  otu_relative_abundance_1 <- cbind(core_otus_demo_count[, 1, drop = FALSE], otu_relative_abundance)
  colnames(otu_relative_abundance_1) <- colnames(core_otus_demo_count)

  # Return the final results as a list
  list(
    final_otu_table_bef_filter = final_otu_table_with_row,
    filtered_md_table = filtered_arab_md_old,
    final_otu_aft_filter = selected_otus,
    normalized_table = selected_otus_1,
    alpha_diversity = alpha_diversity,
    beta_diversity = bd,
    core_otus = top_core_rows,
    non_core_otus = non_core_rows,
    core_otus_tax = common_otu_tax,
    core_otus_count_data = core_otus_demo_count,
    core_otus_relative_abundance = otu_relative_abundance_1
  )
}
