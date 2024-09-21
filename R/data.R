#' OTU Data from HMP2 Gut Dataset
#'
#' This dataset contains operational taxonomic unit (OTU) abundance data from the Human Microbiome Project Phase 2 (HMP2),
#' focusing on gut microbiome samples. The dataset features 1122 samples (rows) across 1953 OTUs (columns).
#'
#' @format A data frame with 1122 rows and 1953 columns.
#' @source Human Microbiome Project Phase 2
#' @usage data(otu_HMP2)
"otu_HMP2"

#' Sample Metadata for HMP2 Gut Dataset
#'
#' This dataset provides the experimental and metadata information corresponding to the samples in the `otu_HMP2` dataset.
#' Each row matches a sample from `otu_HMP2`, and columns include various metadata attributes such as Subject ID, dates,
#' insulin resistance status, and health status indicators.
#'
#' @format A data frame with 1122 rows and 10 variables:
#' \describe{
#'   \item{SubjectID}{Identifier for the subject from whom the sample was taken.}
#'   \item{DateZero}{The number of days passed from the first specimen collected in the study.}
#'   \item{DateSubject}{The number of days passed from the first specimen collected from the subject.}
#'   \item{IRIS}{Indicator of insulin resistance or sensitivity.}
#'   \item{CL4_1, CL4_2, CL1, CL2, CL3, CL4}{Various indicators of the health status of the patient at the time of the specimen.}
#' }
#' @source Human Microbiome Project Phase 2
#' @usage data(info_sample_HMP2)
"info_sample_HMP2"

#' Lineage Information for OTUs in HMP2 Gut Dataset
#'
#' This dataset contains the lineage information for each operational taxonomic 
#' unit (OTU) identified in the `otu_HMP2` dataset from the Human Microbiome 
#' Project Phase 2 (HMP2), focusing on gut microbiome samples. 
#' Each row represents an OTU (1953 in total), and each column represents a 
#' different level of taxonomic classification, from phylum to OTU. 
#' The taxa IDs, which are equivalent to the OTU names in this dataset 
#' (but can differ in general, such as with NCBI where names and IDs are distinct), 
#' are provided as the row names.
#'
#' @format A data frame with 1953 rows and columns for each taxonomic rank from phylum to OTU. The row names are the taxa IDs.
#' @source Human Microbiome Project Phase 2
#' @usage data(info_taxa_HMP2)
"info_taxa_HMP2"
