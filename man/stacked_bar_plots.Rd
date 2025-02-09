% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stacked_bar_plots.R
\name{stacked_bar_plots}
\alias{stacked_bar_plots}
\title{Stacked Bar Plots Based on Relative Abundance Data}
\usage{
stacked_bar_plots(data, num_samples_per_plot)
}
\arguments{
\item{data}{A data frame containing the relative abundance data for the OTUs. The first column should contain the OTU IDs, and the subsequent columns should represent samples.}

\item{num_samples_per_plot}{The number of samples to be displayed in each stacked bar plot.}
}
\value{
A list of interactive stacked bar plots, one for each group of samples, showing the relative abundance of OTUs in the samples.
}
\description{
This function generates stacked bar plots for visualizing
the relative abundance data of different operational taxonomic units (OTUs)
in various samples.
}
\examples{
#To run input data
core_1 <- CoreMicrobiome(
 otu_table = demo_otu,
 tax_table = demo_tax,
 metadata_table = demo_md,
 filter_type = "occupancy_fun_filter", #Or "abundance_fun_filter", Or "combined_filter"
 percent = 0.5,
 method = "css",  # Or "srs", "rrarefy", "tmm", "tmmwsp", "rle", "upperquartile", "none"
 beta_diversity_method = "jaccard",
 top_percentage = 10  # Adjust the percentage as needed for core/non-core OTUs
)
#To run the stacked bar plots function
stacked_plots <- stacked_bar_plots(core_1$core_otus_relative_abundance, 10)
#To view the stacked bar plot
stacked_plots[[1]]
}
