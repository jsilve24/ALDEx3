#' Gut Crohn's microbiome dataset (list)
#'
#' This dataset contains a matrix and a data frame: genus-level microbiome profiles and
#' corresponding sample metadata from a Crohn's disease case–control cohort.
#' The dataset is used in examples and vignettes throughout the package.
#'
#' @format A named list of one matrix and a dataframe:
#' \describe{
#'   \item{\code{counts}}{A matrix with read counts with 195 rows and 95 columns)}
#'   \item{\code{metadata}}{A data frame with 95 rows and 7 columns containing subject-level covariates.}
#' }
#'
#' @details
#' The \code{counts} matrix has one row per sample and one column per genus.
#' The \code{metadata} data frame has one row per sample with metadata, critically Health.status either CD or Control, and Average cell count per gram frozen feces.
#'
#' @source
#' Vandeputte D, Falony G, Vieira-Silva S, Wang J, Sailer M, Theis S, Raes J (2017).
#' "Quantitative microbiome profiling links gut community variation to microbial load".
#' Nature, 551, 507--511. \doi{10.1038/nature24460}
"gut_crohns_data"

#' Oral microbiome perturbation dataset (list)
#'
#' This dataset contains a matrix and a data frame: genus-level microbiome profiles and
#' corresponding sample metadata from an oral microbiome perturbation study.
#' 28 participants' oral microbiomes were measured before, 15 minutes after, and 2 hours
#' after perturbation with either a water control, antiseptic mouthwash, alchohol-free mouthwash,
#' or soda. The dataset is used in examples and vignettes throughout the package.
#'
#' @format A named list of one matrix and a dataframe:
#' \describe{
#'   \item{\code{counts}}{A matrix with read counts with 116 rows and 81 columns)}
#'   \item{\code{metadata}}{A data frame with 81 rows and 38 columns containing sample-level covariates.}
#' }
#'
#' @details
#' The \code{counts} matrix has one row per sample and one column per genus.
#' The \code{metadata} data frame has one row per sample with metadata, critically the participant ID, flow cytometry average (and replicate) cells, time_c (time points), and treat (the perturbations )
#'
#' @source
#' Marotz C, Morton JT, Navarro P, Coker J, Belda-Ferre P, Knight R (2021).
#' "Quantifying live microbial load in human saliva samples over time reveals stable composition and dynamic load".
#' mSystems, 6(3), e01182-21. \doi{10.1128/mSystems.01182-20}
"oral_mouthwash_data"
