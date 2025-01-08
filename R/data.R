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
#' @usage data(meta_HMP2)
"meta_HMP2"

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
#' @usage data(taxa_HMP2)
"taxa_HMP2"

#' MGNet Object from HMP2 Gut Dataset
#'
#' This MGNet object is constructed from the OTU, sample metadata, and taxa lineage information from the Human Microbiome 
#' Project Phase 2 (HMP2) focusing on gut microbiome samples. It integrates these datasets into a structured object suitable 
#' for metagenomic network analysis, providing a comprehensive framework for analyzing microbial communities and their interactions.
#'
#' The object includes:
#' \itemize{
#'   \item \strong{abun}: A numeric matrix representing the raw abundance data of OTUs across samples.
#'   \item \strong{rela}: A numeric matrix representing the relative abundance of each OTU within each sample.
#'   \item \strong{norm}: A numeric matrix of normalized abundance data, suitable for various types of analyses.
#'   \item \strong{meta}: A data.frame containing metadata for samples, which includes details such as subject ID, insulin resistance status, and other health indicators.
#'   \item \strong{taxa}: A data.frame providing taxonomic classification for each OTU from phylum to species level.
#'   \item \strong{netw}: An `igraph` object representing the network of taxa interactions, based on co-occurrence and other metrics.
#'   \item \strong{comm}: An object storing community detection results, facilitating analysis of microbial community structure.
#' }
#'
#' This object allows researchers to perform advanced metagenomic analyses, such as identifying key taxa, exploring community structures, and understanding the ecological roles of different microbes within the gut microbiome.
#'
#' @format An S4 object of class `mgnet` with multiple slots containing matrices and data.frames structured as described above.
#' @source Human Microbiome Project Phase 2
#' @usage data(HMP2)
"HMP2"

#' MGNet Object with Network Data for Subject HMP2
#'
#' This `mgnet` object contains metagenomic data for a specific subject from the Human Microbiome Project Phase 2 (HMP2),
#' including normalized abundance data, sample metadata, taxonomic information, network interactions, and community
#' detection results. It is designed for comprehensive analysis of microbial communities within the gut microbiome of the subject.
#'
#' @format An S4 object of class `mgnet` with:
#' \itemize{
#'   \item \strong{abun}: abundance matrix for 376 taxa across 51 samples.
#'   \item \strong{rela}: Associated relative abundances.
#'   \item \strong{norm}: CLR transformed abundance matrix.
#'   \item \strong{meta}: Sample metadata including SubjectID, DateZero, DateSubject, IRIS, etc.
#'   \item \strong{taxa}: Detailed taxonomic information for each taxon (domain, phylum, class, order, etc.).
#'   \item \strong{netw}: An `igraph` object representing microbial interactions with 2128 edges and a 3.02% density.
#'   \item \strong{comm}: Community detection results identifying 12 communities with detailed community sizes and 121 isolated nodes.
#' }
#' 
#' @source Human Microbiome Project Phase 2
#' @usage data(subject_HMP2_netw)
"mg"

#' MGNet List with Network Data for Selected Subjects from HMP2 Dataset
#'
#' This `mgnetList` object contains `mgnet` instances for two selected subjects from the Human Microbiome Project Phase 2 (HMP2),
#' focusing on their gut microbiome samples with network data and community analyses. The list facilitates comparative
#' metagenomic studies, highlighting interactions and community structures unique to each subject.
#'
#' @format A `mgnetList` object containing 2 `mgnet` instances:
#' \describe{
#'   \item{69-001}{An `mgnet` object with 185 samples, 377 taxa, a network of 3460 edges (~4.88% density), and 7 detected communities (21.75% isolated nodes).}
#'   \item{69-053}{An `mgnet` object with 40 samples, 375 taxa, a network of 174 edges (~0.25% density), and 18 detected communities (74.40% isolated nodes).}
#' }
#' 
#' @details Each `mgnet` object within the list is structured to include normalized abundance data, comprehensive taxa and sample metadata,
#' networks with their respective density and detected communities indicating how many taxa are isolated or grouped. This setup provides
#' an in-depth look at the microbial ecosystem of each subject under study.
#'
#' @source Human Microbiome Project Phase 2
#' @usage data(mgl)
"mgl"
