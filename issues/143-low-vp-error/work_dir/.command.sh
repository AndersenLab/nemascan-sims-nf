#!/bin/bash -ue
if [[ inbred == "inbred" ]]; then
    COMMAND='--fastGWA-mlm-exact'
    GWA_THREADS=1  # pinned: BLAS reduction order must be deterministic

    # Inbred: create sparse GRM approximation (supported by fastGWA-mlm-exact)
    gcta64 --grm TO_SIMS_1_43_0.4_0.05_gamma_ct.full_gcta_grm_inbred \
        --make-bK-sparse 0.05 \
        --out 1_43_0.4_0.05_gamma_ct.full_sparse_grm_inbred \
        --thread-num 1

    GRM_OPTION='--grm-sparse'
    GRM_PREFIX=1_43_0.4_0.05_gamma_ct.full_sparse_grm_inbred
else
    COMMAND="--mlma-loco"
    GWA_THREADS=4  # mlma-loco: per-chromosome refits are independent

    # LOCO: must use the full (non-sparse) GRM — --mlma-loco decomposes per-chromosome
    # covariance from the full GRM; a sparse approximation causes REML divergence.
    GRM_OPTION="--grm"
    GRM_PREFIX=TO_SIMS_1_43_0.4_0.05_gamma_ct.full_gcta_grm_inbred
fi

if [[ pca == "pca" ]]; then

    gcta64 --grm TO_SIMS_1_43_0.4_0.05_gamma_ct.full_gcta_grm_inbred \
        --pca 1 \
        --out 1_43_0.4_0.05_gamma_ct.full_pca \
        --thread-num 1  # pinned: BLAS reduction order must be deterministic

    COVAR="--qcovar 1_43_0.4_0.05_gamma_ct.full_pca.eigenvec"
else
    COVAR=""
fi

awk '{print $2}' TO_SIMS_1_43_0.4_0.05_gamma_ct.full.bim > plink_snplist.txt

gcta64 ${COMMAND} \
    --bfile TO_SIMS_1_43_0.4_0.05_gamma_ct.full \
    ${GRM_OPTION} ${GRM_PREFIX} \
    ${COVAR} \
    --out 1_43_0.4_0.05_gamma_ct.full_lmm-exact_inbred_pca \
    --pheno 1_43_0.4_0.05_gamma_ct.full_sims.pheno \
    --extract plink_snplist.txt \
    --thread-num ${GWA_THREADS}

if [[ inbred == "loco" ]]; then
    mv "1_43_0.4_0.05_gamma_ct.full_lmm-exact_inbred_pca.loco.mlma" \
       "1_43_0.4_0.05_gamma_ct.full_lmm-exact_inbred_pca.mlma"
fi

cat <<-END_VERSIONS > versions.yml
"GCTA_PERFORM_GWA":
    GCTA: $( gcta64 --version |& grep version |& cut -f 3 )
END_VERSIONS
