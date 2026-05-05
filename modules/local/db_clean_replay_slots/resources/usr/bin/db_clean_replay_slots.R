#!/usr/bin/env Rscript

library(optparse)

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--replay_tsv", type = "character"),
    optparse::make_option("--db_root",    type = "character"),
    optparse::make_option("--cv_maf",     type = "double", default = NA_real_,
        help = "Global cv_maf (NA -> per-group ms_maf from marker_set_metadata)"),
    optparse::make_option("--cv_ld",      type = "double", default = 0.8)
)))

R_SOURCE_DIR <- Sys.getenv("R_SOURCE_DIR")
source(file.path(R_SOURCE_DIR, "setup.R"))

#' Read and validate the replay manifest TSV
#'
#' @param path Character. Path to `replay.tsv`.
#' @return A tibble with one row per failed slot.
read_replay_manifest <- function(path) {
    readr::read_tsv(path, col_types = readr::cols(
        session    = readr::col_character(),
        task_hash  = readr::col_character(),
        attempt    = readr::col_integer(),
        max_retries = readr::col_integer(),
        species    = readr::col_character(),
        group      = readr::col_character(),
        rep        = readr::col_integer(),
        maf        = readr::col_double(),
        nqtl       = readr::col_character(),
        effect     = readr::col_character(),
        h2         = readr::col_double(),
        mode       = readr::col_character(),
        type       = readr::col_character(),
        exit       = readr::col_integer()
    ))
}

#' Enumerate the four trait parquet file paths for a trait_id
#'
#' @param db_root Character. Database root directory.
#' @param trait_id Character. Trait identifier hash.
#' @return A character vector of four absolute file paths.
trait_files_for_id <- function(db_root, trait_id) {
    c(
        fs::path(db_root, "traits",                  glue::glue("{trait_id}.parquet")),
        fs::path(db_root, "traits/causal_variants",  glue::glue("{trait_id}_causal.parquet")),
        fs::path(db_root, "traits/causal_genotypes", glue::glue("{trait_id}_causal_geno.parquet")),
        fs::path(db_root, "traits/phenotypes",       glue::glue("{trait_id}_phenotype.parquet"))
    )
}

#' Enumerate all mapping parquet file paths for a trait_id across all mode x pca combos
#'
#' @param db_root Character. Database root directory.
#' @param trait_id Character. Trait identifier hash.
#' @param group Character. Population/group identifier.
#' @return A character vector of eight file paths (4 mapping_ids x 2 files each).
mapping_files_for_id <- function(db_root, trait_id, group) {
    combos <- tibble::tibble(
        mode = c("inbred", "inbred", "loco", "loco"),
        pca  = c(TRUE,     FALSE,    TRUE,   FALSE)
    )

    purrr::pmap(combos, function(mode, pca) {
        mid  <- generate_mapping_id(trait_id, mode, pca)$hash
        slot <- fs::path(db_root, "mappings",
                         glue::glue("population={group}"),
                         glue::glue("mapping_id={mid}"))
        c(fs::path(slot, "data.parquet"), fs::path(slot, "meta.parquet"))
    }) |>
        purrr::list_c()
}

#' Delete a set of files, skipping any that do not exist (idempotent)
#'
#' @param files Character vector of file paths to delete.
#' @return Invisibly returns `files`.
delete_if_exists <- function(files) {
    extant <- files[fs::file_exists(files)]
    purrr::walk(extant, function(f) {
        message("Removing ", f)
        fs::file_delete(f)
    })
    invisible(files)
}

#' Derive trait_id and delete all DB files for one replay-manifest slot
#'
#' @param species Character. Species identifier. Present in replay.tsv since
#'   the species-in-fanout-ids update; not required by read_marker_set_metadata()
#'   (keyed by group/maf) but kept in the signature to surface schema mismatches.
#' @param group Character. Population group.
#' @param rep Integer. Replicate number.
#' @param maf Double. Minor allele frequency.
#' @param nqtl Character. Number of QTLs.
#' @param effect Character. Effect size distribution.
#' @param h2 Double. Heritability.
#' @param db_root Character. Database root directory.
#' @param ... Additional manifest columns (ignored).
#' @return Invisibly returns NULL.
clean_slot <- function(species, group, rep, maf, nqtl, effect, h2, db_root,
                       cv_maf = NA_real_, cv_ld = 0.8, ...) {
    ms_meta <- read_marker_set_metadata(group, maf, base_dir = db_root)
    if (is.null(ms_meta)) {
        stop(glue::glue(
            "No marker_set_metadata for group={group} maf={maf} — cannot derive trait_id"
        ))
    }

    cv_maf_effective <- if (is.na(cv_maf)) maf else cv_maf

    ids <- build_ids_from_params(
        params = list(
            population        = group,
            maf               = maf,
            species           = ms_meta$species,
            vcf_release_id    = ms_meta$vcf_release_id,
            ms_ld             = ms_meta$ms_ld,
            nqtl              = nqtl,
            effect            = effect,
            rep               = rep,
            h2                = h2,
            cv_maf_effective  = cv_maf_effective,
            cv_ld             = cv_ld
        ),
        mode = "inbred", pca = FALSE
    )
    trait_id <- ids$trait_id$hash

    c(
        trait_files_for_id(db_root, trait_id),
        mapping_files_for_id(db_root, trait_id, group)
    ) |>
        delete_if_exists()

    invisible(NULL)
}

# --- main ---

manifest <- read_replay_manifest(opts$replay_tsv)
purrr::pwalk(manifest, clean_slot, db_root = opts$db_root,
             cv_maf = opts$cv_maf, cv_ld = opts$cv_ld)
message(glue::glue("DB_CLEAN_REPLAY_SLOTS: cleanup complete for {nrow(manifest)} slot(s)."))
