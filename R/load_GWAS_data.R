#' Load GWAS data
#'
#' @param files The file or files containing GWAS data as a vector
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @param input_type single, two, TASSEL
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @export
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_GWAS_data(demo_association_file, demo_effects_file)
load_GWAS_data <- function(files,
                           trait = "Trait",
                           marker = "Marker",
                           locus = "Locus",
                           site = "Site",
                           p = "p",
                           marker_R2 = "marker_R2",
                           effect = "Effect",
                           input_type) {
  if (length(files) == 1) {
    gwas_file = files[[1]]
    gwas_data = load_single_generic(gwas_file,
                                    trait,
                                    marker,
                                    locus,
                                    site,
                                    p,
                                    marker_R2,
                                    effect)
    
  } else if (length(files) == 2) {
    association_file = files[[1]]
    effects_file = files[[2]]
    if (input_type == "TASSEL") {
      gwas_data = load_TASSEL(files[[1]], 
                              files[[2]],
                              trait,
                              marker,
                              locus,
                              site,
                              p,
                              marker_R2,
                              effect)
    }
  }
  gwas_data
}

#' Load single file data
#'
#' @param gwas_file The file containing GWAS data
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_gwas_file = system.file("extdata", "single_gwas.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_single_generic(c(demo_gwas_file))
load_single_generic <- function(gwas_file,
                                trait = "Trait",
                                marker = "Marker",
                                locus = "Locus",
                                site = "Site",
                                p = "p",
                                marker_R2 = "marker_R2",
                                effect = "Effect") {
  
  gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2),
                  effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2,
                  .data$effect)
  
  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- gwas_data %>%
    dplyr::group_by(.data$Marker_original) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  gwas_data <-
    gwas_data %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  
  # Remove all NaN data to prevent math with NaN
  gwas_data <- gwas_data %>% dplyr::filter(.data$marker_R2 != "NaN")
  gwas_data
}

#' Load two-file generic
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_two_generic(c(demo_association_file, demo_effects_file))
load_two_generic <- function(association_file,
                             effects_file,
                             trait = "Trait",
                             marker = "Marker",
                             locus = "Locus",
                             site = "Site",
                             p = "p",
                             marker_R2 = "marker_R2",
                             effect = "Effect") {
  
  stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)
  
  effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$Effect)
  
  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker_original) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  stats <-
    stats %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  
  # Remove all NaN data to prevent math with NaN
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")
  
  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker_original") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Marker_original,
      .data$Chr,
      .data$Pos,
      .data$p,
      .data$marker_R2,
      .data$Effect.x,
      .data$Effect.y
    )
  all_data
}


#' Load TASSEL data
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_TASSEL(c(demo_association_file, demo_effects_file))
load_TASSEL <- function(association_file,
                        effects_file,
                        trait = "Trait",
                        marker = "Marker",
                        locus = "Locus",
                        site = "Site",
                        p = "p",
                        marker_R2 = "marker_R2",
                        effect = "Effect") {
  
  stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)
  
  effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$Effect)
  
  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker_original) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  stats <-
    stats %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  
  # Remove all NaN data to prevent math with NaN
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")
  
  # Split effects into even and odd rows and
  # recombine into a single row without duplicate columns
  odd_effects <- effects[seq(1, nrow(effects), by = 2), ]
  even_effects <- effects[seq(2, nrow(effects), by = 2), ]
  effects <- merge(odd_effects, even_effects, by = "Marker_original")
  effects <- dplyr::mutate(
    effects,
    Trait = effects$Trait.x,
    Trait.x = NULL,
    Trait.y = NULL
  )
  
  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker_original") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Marker_original,
      .data$Chr,
      .data$Pos,
      .data$p,
      .data$marker_R2,
      .data$Effect.x,
      .data$Effect.y
    )
  all_data
}
