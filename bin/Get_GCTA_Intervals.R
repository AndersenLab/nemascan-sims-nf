#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(ggbeeswarm)

# argument information
# 1 - Genetoype matrix
# 2 - Phenotype data
# 3 - Mapping data
# 4 - independent tests - eigen
# 5 - number of qtl simulated (numeric)
# 6 - simulation replicate (numeric)
# 7 - If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)
# 8 - Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)
# 9 - simulated heritability2
# 10 - - String, What significance threshold to use for defining QTL,
#       BF = Bonferroni,
#       EIGEN = Defined by Number of independent tests from Eigen Decomposition of SNV correlation matrix,
#       or a user defined number
# 11 - strain set name (string)
# 12 - minor allele frequency (numeric)
# 13 - effect size range (character)
# 14 - mode

# load arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("complete_0.05_Genotype_Matrix.tsv",  ## TESTING
#           "5_49_0.2_0.05_gamma_complete_sims.phen",
#           "5_49_0.2_0.05_gamma_complete_lmm-exact_inbred.fastGWA",
#           "complete_0.05_total_independent_tests.txt", 5, 49, 1000, 150, 0.2, 0.05,
#           "BF", "complete", 0.05, "gamma", "inbred")

# define the trait name
print("Defining trait name...")

trait_name <- glue::glue("{args[5]}_{args[6]}_{args[9]}")

print(
  glue::glue("Trait name defined as: {trait_name}")
)

# load phenotpe data
print(
  glue::glue("Loading phenotype data from {args[2]}...")
)

phenotype_data <- data.table::fread(args[2], col.names = c("strain", "strain_drop", trait_name)) %>%
  na.omit() %>%
  dplyr::select(-strain_drop) %>%
  as.data.frame()

print(
  glue::glue(
    "Phenotype data loaded with {nrow(phenotype_data)} rows."
  )
)

# load GCTA mapping data

print(
  glue::glue("Loading GCTA mapping data from {args[3]}...")
)

if (args[14] == "inbred" | args[14] == "inbred_pca") {
  print(
    glue::glue("Mapping data is in inbred format, processing accordingly...")
  )

  map_df <- data.table::fread(args[3]) %>%
    dplyr::rename(marker = SNP, CHROM = CHR) %>%
    dplyr::mutate(log10p = -log10(P))

  print(
    glue::glue("Inbred mapping data loaded with {nrow(map_df)} rows.")
  )
} else {
  print(
    glue::glue("Mapping data is in LOCO format, processing accordingly...")
  )

  map_df <- data.table::fread(args[3]) %>%
    dplyr::rename(
      marker = SNP,
      CHROM = Chr,
      POS = bp,
      AF1 = Freq,
      BETA = b,
      SE = se,
      P = p
    ) %>%
    dplyr::mutate(log10p = -log10(P))

  print(
    glue::glue("LOCO mapping data loaded with {nrow(map_df)} rows.")
  )
}

# load genotype matrix
print(
  glue::glue("Loading genotype matrix from {args[1]}...")
)
genotype_matrix <- readr::read_tsv(args[1]) %>%
  na.omit()
print(
  glue::glue(
    "Genotype matrix loaded with {nrow(genotype_matrix)} rows and {ncol(genotype_matrix)} columns."
  )
)

# define method for setting significance threshold
significance_threshold <- args[10]

print(
  glue::glue("Setting significance threshold to: {significance_threshold}")
)

# set significance threshold
if (significance_threshold == "EIGEN") {
  QTL_cutoff <- data.table::fread(args[4]) %>% dplyr::pull(V1) # Independent tests should be fed to calculation
} else if (significance_threshold == "BF") {
  QTL_cutoff <- NA
} else {
  QTL_cutoff <- as.numeric(args[10]) # Specified Threshold
}


# process mapping function
process_mapping_df <- function(mapping_df,
                               phenotype_df,
                               CI_size = as.numeric(args[8]),
                               snp_grouping = as.numeric(args[7]),
                               BF = NA,
                               thresh = significance_threshold,
                               geno = genotype_matrix) {
  print("Processing mapping data...")

  print("Setting up phenotype data...")
  pheno <- phenotype_df

  pheno$trait <- colnames(phenotype_df)[2]

  colnames(pheno) <- c("strain", "value", "trait")
  print("Phenotype data set up.")

  print("Flagging significant markers with QTL_cutoff...")
  # Determine how to make threshold
  if (is.na(QTL_cutoff)) {
    print(
      "Using Bonferroni correction for significance threshold to flag significant markers."
    )
    mapping_df <- mapping_df %>%
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>%
      dplyr::filter(log10p != 0) %>%
      dplyr::mutate(BF = -log10(0.05 / sum(log10p > 0, na.rm = T))) %>%
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  } else if (is.numeric(QTL_cutoff) & thresh == "EIGEN") {
    print(
      "Using EIGEN significance threshold to flag significant markers."
    )
    mapping_df <- mapping_df %>%
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>%
      dplyr::filter(log10p != 0) %>%
      dplyr::mutate(BF = -log10(0.05 / BF)) %>%
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  } else {
    print(
      "Using user-defined significance threshold to flag significant markers."
    )
    mapping_df <- mapping_df %>%
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>%
      dplyr::filter(log10p != 0) %>%
      dplyr::mutate(BF = BF) %>%
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  }

  Processed <- mapping_df %>%
    dplyr::filter(sum(aboveBF, na.rm = T) > 0) %>%
    dplyr::ungroup()

  print("Pulling significant markers for variance explained (VE) calculation...")

  snpsForVE <- Processed %>%
    dplyr::filter(aboveBF == 1) %>%
    dplyr::select(marker, trait)

  print(
    glue::glue("Found {nrow(snpsForVE)} significant markers from {nrow(Processed)} markers.")
  )

  snpsForVE$trait <- as.character(snpsForVE$trait)

  if (nrow(snpsForVE) > nrow(Processed) * 0.15) {
    print(
      glue::glue("Too many significant markers QTL intervals are not flagged.")
    )

    Processed <- mapping_df %>%
      dplyr::mutate(
        strain = NA, value = NA, allele = NA, var.exp = NA,
        startPOS = NA, peakPOS = NA, endPOS = NA,
        peak_id = NA, interval_size = NA
      )
  } else if (nrow(snpsForVE) > 0 && nrow(snpsForVE) < nrow(Processed) * 0.15) {
    print(
      "The number of significant markers is acceptable, proceeding to define QTL and calculate Ve."
    )

    print(
      "Joining phenotype data with significant markers"
    )
    row.names(pheno) <- gsub("-", "\\.", row.names(pheno))

    pheno$trait <- gsub("-", "\\.", pheno$trait)

    rawTr <- pheno %>%
      dplyr::left_join(., snpsForVE, by = "trait")

    rawTr$marker <- as.character(rawTr$marker)
    rawTr$strain <- as.character(rawTr$strain)

    print(
      "Joined phenotypes and significant SNPs"
    )

    print(
      "Joining genotype data with significant markers"
    )

    snp_df <- geno %>%
      dplyr::select(-REF, -ALT)

    gINFO <- snp_df %>%
      dplyr::mutate(marker = paste(CHROM, POS, sep = ":")) %>%
      dplyr::filter(marker %in% snpsForVE$marker) %>%
      tidyr::gather(strain, allele, -marker, -CHROM, -POS)

    gINFO$marker <- as.character(gINFO$marker)
    gINFO <- suppressWarnings(data.frame(gINFO) %>%
      dplyr::left_join(., snpsForVE, by = "marker") %>%
      dplyr::left_join(rawTr, ., by = c("trait", "strain", "marker")))

    print(
      "Joined genotyped data with significant markers and phenotypes"
    )

    print(
      "Calculating variance explained (VE) for each marker"
    )

    cors <- gINFO %>%
      dplyr::group_by(trait, marker) %>%
      dplyr::mutate(var.exp = cor(value, allele,
        use = "pairwise.complete.obs",
        method = "pearson"
      )^2)

    print(
      "Variance explained (VE) calculated for each marker."
    )

    print(
      "Joining variance explained (VE) with all marker data"
    )

    CORmaps <- Processed %>%
      dplyr::left_join(., cors, by = c("trait", "marker", "CHROM", "POS"), copy = TRUE)

    print(
      "Joined all marker data with variance explained (VE)."
    )

    processed_mapping_df <- Processed
    correlation_df <- CORmaps
    phenotypes <- as.character(unique(processed_mapping_df$trait))
    intervals <- list()

    print(
      glue::glue("Processing {length(phenotypes)} phenotypes for QTL intervals...")
    )

    # [ ] - Current pipeline will never have more than one trait, so this could be simplified
    for (i in 1:length(phenotypes)) {
      print(
        glue::glue("Processing phenotype {phenotypes[i]}...")
      )

      print(
        "Creating PeakDF"
      )

      # [ ] - Several possibly redundant steps here? Why group create the index multiple times?
      PeakDF <- processed_mapping_df %>%
        dplyr::filter(trait == phenotypes[i]) %>%
        # group makers by chromosome
        dplyr::group_by(CHROM, trait) %>%
        # create an index for SNP based on row number
        dplyr::mutate(index = 1:n()) %>%
        # for each group count the number of markers above sig threshold
        dplyr::mutate(peaks = cumsum(aboveBF)) %>%
        # filter to just sig. markers
        dplyr::filter(aboveBF == 1) %>%
        # group after by filtering
        dplyr::group_by(CHROM, trait) %>%
        # for each group calculate the number of significant markers
        dplyr::mutate(nBF = n()) %>%
        # group again - why?
        dplyr::group_by(CHROM, trait) %>%
        # arrange by chromosome and position
        dplyr::arrange(CHROM, POS)

      print(
        glue::glue("PeakDF created with {nrow(PeakDF)} rows for phenotype {phenotypes[i]}.")
      )
      print(
        "PeakDF structure:"
      )
      print(PeakDF)

      print(
        "Creating SNP index for all markers"
      )

      SNPindex <- processed_mapping_df %>%
        dplyr::filter(trait == phenotypes[i]) %>%
        dplyr::group_by(CHROM, trait) %>%
        dplyr::mutate(index = 1:n()) %>%
        dplyr::distinct(CHROM, POS, .keep_all = T) %>%
        dplyr::select(CHROM, POS, index) %>%
        dplyr::group_by(CHROM) %>% # add because we want the min and max position for each chrom
        dplyr::filter(POS == min(POS) | POS == max(POS))

      print(
        glue::glue("SNP index created with {nrow(SNPindex)} rows for phenotype {phenotypes[i]}.")
      )

      print(
        "Finding peaks"
      )


      findPks <- PeakDF %>%
        dplyr::filter(trait == phenotypes[i]) %>%
        dplyr::group_by(CHROM) %>%
        dplyr::arrange(CHROM, POS)

      print(
        glue::glue("Created findPks with {nrow(findPks)} rows")
      )

      print(
        "Processing peaks to define intervals..."
      )

      if (findPks$nBF == 1 & length(unique(findPks$CHROM)) == 1) {
        print(
          "Only one peak found, creating single interval for this peak."
        )

        print(
          "Assigning peak ID"
        )
        findPks$pID <- 1

        print(
          glue::glue(
            "Creating confidence intervals for the peak of size {CI_size} to marker"
          )
        )

        findPks <- findPks %>%
          dplyr::group_by(CHROM, pID, trait) %>%
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)


        print(
          "Confidence intervals created."
        )

        print(
          "Adjusting start and end positions so they do not exceed SNP index limits."
        )

        for (k in 1:nrow(findPks)) {
          # filter to taget SNPs
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])
          if (findPks$start[k] < min(tSNPs$index)) {
            print(
              glue::glue("Adjusting start position for peak {k} on chromosome {findPks$CHROM[k]} to minimum SNP index.")
            )

            findPks$start[k] <- min(tSNPs$index)

            print(
              glue::glue("Start position for peak {k} on chromosome {findPks$CHROM[k]} adjusted to {findPks$start[k]}.")
            )
          }
          if (findPks$end[k] > max(tSNPs$index)) {
            print(
              glue::glue(
                "Adjusting end position for peak {k} on chromosome {findPks$CHROM[k]} to maximum SNP index."
              )
            )

            findPks$end[k] <- max(tSNPs$index)

            print(
              glue::glue(
                "End position for peak {k} on chromosome {findPks$CHROM[k]} adjusted to {findPks$end[k]}."
              )
            )
          }

          print(
            glue::glue("Finished processing for peak {k} of nrow {nrow(findPks)} peaks on chromosome {findPks$CHROM[k]}.")
          )
        }

        print(
          "Adding peak information to intervals."
        )

        intervals[[i]] <- findPks %>% dplyr::ungroup()

        print(
          "Peak information added to intervals."
        )
      } else {
        print(
          glue::glue("Multiple peaks found, processing {nrow(findPks)} to create intervals.")
        )
        findPks$pID <- 1
        for (j in 2:nrow(findPks)) {
          print(
            glue::glue("Processing peak {j} of {nrow(findPks)} for chromosome {findPks$CHROM[j]}...")
          )
          # if the distance between the current peak and the previous peak is less than the snp_grouping
          # and they are on the same chromosome, assign the same peak ID
          # otherwise, assign a new peak ID
          print(
            glue::glue("Checking if peak {j} is within {snp_grouping} markers of previous peak {j - 1}...")
          )

          findPks$pID[j] <- ifelse(abs(findPks$index[j] - findPks$index[j - 1]) < snp_grouping & findPks$CHROM[j] == findPks$CHROM[j - 1],
            findPks$pID[j - 1],
            findPks$pID[j - 1] + 1
          )
        }
        print(
          glue::glue("Peak ID {findPks$pID[j]} assigned for peak {j}")
        )

        print(
          "Creating confidence intervals for with updated peaks..."
        )

        findPks <- findPks %>%
          dplyr::group_by(CHROM, pID, trait) %>%
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)

        print(
          "Confidence intervals updated."
        )

        print(
          "Adjusting start and end positions so they do not exceed SNP index limits."
        )

        for (k in 1:nrow(findPks)) {
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])

          if (findPks$start[k] < min(tSNPs$index)) {
            print(
              glue::glue("Adjusting start position for peak {k} on chromosome {findPks$CHROM[k]} to minimum SNP index.")
            )


            findPks$start[k] <- min(tSNPs$index)

            print(
              glue::glue("Start position for peak {k} on chromosome {findPks$CHROM[k]} adjusted to {findPks$start[k]}.")
            )
          }
          if (findPks$end[k] > max(tSNPs$index)) {
            print(
              glue::glue(
                "Adjusting end position for peak {k} on chromosome {findPks$CHROM[k]} to maximum SNP index."
              )
            )


            findPks$end[k] <- max(tSNPs$index)

            print(
              glue::glue(
                "End position for peak {k} on chromosome {findPks$CHROM[k]} adjusted to {findPks$end[k]}."
              )
            )
          }
          print(
            glue::glue("Finished processing for peak {k} of nrow {nrow(findPks)} peaks on chromosome {findPks$CHROM[k]}.")
          )
        }
        print(
          "Finished processing peaks."
        )
      }

      print(
        "Adding peak information to intervals."
      )


      intervals[[i]] <- findPks %>% dplyr::ungroup()

      print(
        "Peak information added to intervals."
      )
    }

    print(
      "Combining all intervals into a single data frame."
    )

    intervalDF <- data.table::rbindlist(intervals)

    print(
      glue::glue("Combined intervals data frame created with {nrow(intervalDF)} rows.")
    )

    peak_df <- intervalDF
    peak_list <- intervals

    print(
      "Creating position index reference."
    )

    Pos_Index_Reference <- processed_mapping_df %>%
      dplyr::group_by(CHROM, trait) %>%
      dplyr::mutate(index = 1:n()) %>%
      dplyr::mutate(peaks = cumsum(aboveBF)) %>%
      dplyr::select(trait, CHROM, POS, index) %>%
      dplyr::filter(index %in% c(unique(peak_df$start), unique(peak_df$end))) %>%
      dplyr::ungroup()

    Pos_Index_Reference$trait <- as.character(Pos_Index_Reference$trait)


    print(
      glue::glue("Position index reference created with {nrow(Pos_Index_Reference)} rows.")
    )


    print(
      "Creating interval positions for each peak..."
    )
    interval_positions <- list()
    for (i in 1:length(peak_list)) {
      print(paste(100 * signif(i / length(peak_list), 3),
        "%",
        sep = ""
      ))
      peak_list[[i]]$trait <- as.character(peak_list[[i]]$trait)
      peak_list[[i]] <- peak_list[[i]] %>%
        dplyr::arrange(desc(log10p)) %>%
        dplyr::distinct(pID, .keep_all = T)
      trait_i <- unique(peak_list[[i]]$trait)
      index_i <- c(peak_list[[i]]$start, peak_list[[i]]$end)
      CHROM_i <- peak_list[[i]]$CHROM
      PKpos <- data.frame(Pos_Index_Reference) %>%
        dplyr::filter(trait == trait_i & index %in% index_i & CHROM %in% CHROM_i) %>%
        dplyr::left_join(., peak_list[[i]], by = c("trait", "CHROM")) %>%
        dplyr::mutate(issues = ifelse(start == index.x | end == index.x, 1, 0)) %>%
        dplyr::filter(issues != 0) %>%
        dplyr::select(trait, CHROM, POS.x, POS.y, pID, log10p, index.x, index.y, start, end) %>%
        dplyr::group_by(CHROM, pID) %>%
        dplyr::mutate(startPOS = min(POS.x), peakPOS = POS.y, endPOS = max(POS.x)) %>%
        dplyr::distinct(trait, CHROM, pID, peakPOS, .keep_all = T) %>%
        dplyr::select(trait, CHROM, POS = POS.y, startPOS, peakPOS, endPOS, peak_id = pID)
      interval_positions[[i]] <- PKpos
    }
    print(
      "Interval positions created for each peak."
    )

    print(
      glue::glue("Combining interval positions into a single data frame.")
    )

    interval_pos_df <- data.frame(data.table::rbindlist(interval_positions)) %>%
      dplyr::mutate(interval_size = endPOS - startPOS)

    print(
      glue::glue("Combined interval positions data frame created with {nrow(interval_pos_df)} rows.")
    )
    print(
      "Joining Ve correlation data with interval positions..."
    )

    Processed <- suppressWarnings(dplyr::left_join(correlation_df,
      interval_pos_df,
      by = c("trait", "CHROM", "POS"),
      copy = TRUE
    ))

    print(
      "Joined Ve correlation data with interval positions."
    )
  } else {
    print(
      "No significant peaks found. Returning empty data frame."
    )

    Processed <- mapping_df %>%
      dplyr::mutate(
        strain = NA, value = NA, allele = NA, var.exp = NA,
        startPOS = NA, peakPOS = NA, endPOS = NA,
        peak_id = NA, interval_size = NA
      )
  }

  return(Processed)
}


# process mapping data, define QTL

print("Running process_mapping_df function to process mapping data...")

processed_mapping <- process_mapping_df(
  mapping_df = map_df,
  phenotype_df = phenotype_data,
  CI_size = as.numeric(args[8]),
  snp_grouping = as.numeric(args[7]),
  BF = QTL_cutoff,
  thresh = significance_threshold,
  geno = genotype_matrix
)

print("Finished processing mapping data.")

# save processed mapping data

print("Saving processed mapping data...")

label <- glue::glue("LMM-EXACT_{args[14]}")
readr::write_tsv(processed_mapping,
  file = glue::glue("{trait_name}_{args[12]}_{args[13]}_{args[11]}_processed_{label}_mapping.tsv"),
  col_names = T
)

print("Saved processed mapping data.")

# extract interval information
print("Extracting interval information...")

qtl_region <- processed_mapping %>%
  na.omit() %>%
  dplyr::distinct(CHROM, marker, trait, startPOS, peakPOS, endPOS, peak_id) %>%
  dplyr::mutate(algorithm = args[15])

print(
  glue::glue("Extracted interval information with {nrow(qtl_region)} rows.")
)

# save processed mapping data
print("Saving QTL region data...")

readr::write_tsv(qtl_region,
  file = glue::glue("{trait_name}_{args[12]}_{args[13]}_{args[11]}_{label}_qtl_region.tsv"),
  col_names = T
)

print("Saved QTL region data.")

print("Finished.")

# ## LD ###
# interval_pos_df_LD <- interval_pos_df %>%
#   dplyr::mutate(marker = paste(CHROM, POS, sep = ":"))
# if(nrow(interval_pos_df_LD) == 1){
#
#   marker.LD <- data.frame(interval_pos_df_LD$marker, interval_pos_df_LD$marker, NA, unique(interval_pos_df_LD$trait)) %>%
#     `colnames<-`(c("marker1","marker2","r2","trait"))
#
# } else {
#   QTLcombos <- data.frame(t(combn(x = unique(gINFO$marker), m = 2))) %>%
#     `colnames<-`(c("marker1","marker2"))
#   LD <- list()
#
#   for(q in 1:nrow(QTLcombos)){
#     markers.of.interest <- c(as.character(QTLcombos[q,]$marker1),
#                              as.character(QTLcombos[q,]$marker2))
#     haps <- suppressMessages(gINFO %>%
#                                dplyr::select(strain, allele, marker) %>%
#                                dplyr::filter(marker %in% markers.of.interest) %>%
#                                dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
#                                dplyr::mutate(marker = as.factor(marker)) %>%
#                                tidyr::pivot_wider(names_from = marker, values_from = allele) %>%
#                                `colnames<-`(c("strain","A","B")) %>%
#                                tidyr::unite("hap", c(A,B), sep = "_", remove = F))
#
#     P <- suppressMessages(haps %>%
#                             dplyr::group_by(hap) %>%
#                             dplyr::summarise(n()/nrow(haps)) %>%
#                             `colnames<-`(c("P","freq")) %>%
#                             tidyr::pivot_wider(names_from = P, values_from = freq))
#
#     n <- suppressMessages(gINFO %>%
#                             dplyr::select(marker) %>%
#                             dplyr::mutate(marker = as.factor(marker)) %>%
#                             dplyr::group_by(marker) %>%
#                             dplyr::summarise(n()) %>%
#                             `colnames<-`(c("marker_id","total")))
#
#     p <- suppressMessages(gINFO %>%
#                             dplyr::select(strain, allele, marker) %>%
#                             dplyr::filter(marker %in% markers.of.interest) %>%
#                             dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
#                             dplyr::mutate(marker = as.factor(marker)) %>%
#                             dplyr::group_by(allele, marker) %>%
#                             dplyr::summarise(n()) %>%
#                             `colnames<-`(c("p","marker_id","n")) %>%
#                             dplyr::left_join(.,n) %>%
#                             dplyr::mutate(freq = n/total) %>%
#                             tidyr::unite("allele", c(p,marker_id), sep = "_") %>%
#                             dplyr::ungroup() %>%
#                             dplyr::select(allele, freq) %>%
#                             tidyr::pivot_wider(names_from = allele, values_from = freq))
#
#     pApB <- p %>%
#       dplyr::select(contains("REF")) %>%
#       c(.) %>%
#       unlist() %>%
#       prod()
#
#     pApBpapb <- p %>%
#       unlist() %>%
#       prod()
#
#     D.AB <- P$REF_REF - pApB
#     r2 <- (D.AB^2)/pApBpapb
#     LD[[q]] <- gINFO %>%
#       dplyr::select(marker) %>%
#       dplyr::filter(marker %in% markers.of.interest) %>%
#       dplyr::distinct()  %>%
#       tidyr::pivot_wider(names_from = marker, values_from = marker) %>%
#       `colnames<-`(c("marker1","marker2")) %>%
#       dplyr::mutate(r2 = r2,
#                     trait = unique(gINFO$trait))
#   }
#
#   marker.LD <- Reduce(rbind, LD)
# }

## save mapping LD data

# readr::write_tsv(processed_mapping[[2]],
#                  file = glue::glue("{trait_name}_{args[12]}_{args[13]}_{args[11]}_{args[14]}_qtl_LD.tsv"),
#                  col_names = T)
