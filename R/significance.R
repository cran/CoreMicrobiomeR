#' Testing the Significance of the Identified Core Microbiome
#'
#' @description This function performs a two-sample variance test to assess the
#' statistical significance of differences in abundance between core OTUs and
#' non-core OTUs. It takes two data frames as input, representing the abundance
#' of core OTUs and non-core OTUs, and returns the results of the variance test.
#' It tells whether the identified core represents the particular environment or
#' habitat.
#'
#' @usage significance(core_ids, non_core_ids)
#'
#' @param core_ids A Dataframe of core OTUs where the first row is the OTU ID and column names refer to sites/sample names
#' @param non_core_ids A Dataframe of non_core OTUs where the first row is the OTU ID and column names refer to sites/sample names
#'
#' @returns
#'  This function gives the list  which consist of following results.
#'
#' `statistic`  Calculated F test statistic
#'
#' `parameter` The numerator degrees of freedom (num df), and the denominator degrees of freedom (denom df)
#'
#' `p-value` Probability value
#'
#' `alternative` The alternative hypothesis for this test is that the true ratio of variances is not equal to 1. This suggests that the variances of the two data sets are different
#'
#' `conf.int` 95 percent confidence interval limit for the ratio of variances
#'
#' `estimate` Ratio of variances between core_data and non_core_data calculated
#'
#' `method` The test performed is an F test, which compares the variances of the two data sets
#'
#' `data.name` The data used for the test are core_ids and non_core_ids
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
#' #To run significance test
#' f_test <- significance(core_1[["core_otus"]] , core_1[["non_core_otus"]] )
#'
#' #To view the significance test result
#' f_test

significance <- function(core_ids, non_core_ids) {
  # Extract the first column of each data frame (excluding the first column)
  core_data <- as.vector(data.matrix(core_ids[,-1]))
  non_core_data <- as.vector(data.matrix(non_core_ids[,-1]))

  # Perform the two-sample variance test
  vt_result <- stats::var.test(core_data, non_core_data, alternative = "two.sided")

  return(vt_result)
}
