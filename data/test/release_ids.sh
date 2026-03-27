#!/usr/bin/env bash
#
# Canonical CaeNDR release IDs for nemascan-sims-nf test VCFs.
#
# These are the release dates currently pinned in the test strainfiles:
#   data/test/test_strains.txt                 → test_ce.CE_VCF_RELEASE
#   data/test/test_strains_three_species.txt   → test_ce.CE_VCF_RELEASE
#                                                 test_cb.CB_VCF_RELEASE
#                                                 test_ct.CT_VCF_RELEASE
#
# When regenerating test VCFs from a new CaeNDR release, update all three
# variables here, re-run the generate scripts, update the strainfiles with
# the new dated symlink names, and commit the new symlinks.
#
# Source VCF URL pattern (replace YYYYMMDD and {species}):
#   https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/{species}/YYYYMMDD/variation/WI.YYYYMMDD.hard-filter.isotype.vcf.gz
#
# Usage (optional — source this file to export the variables):
#   source data/test/release_ids.sh

CE_VCF_RELEASE=20250625   # C. elegans  — test_ce.20250625.vcf.gz
CB_VCF_RELEASE=20250626   # C. briggsae — test_cb.20250626.vcf.gz
CT_VCF_RELEASE=20250627   # C. tropicalis — test_ct.20250627.vcf.gz
