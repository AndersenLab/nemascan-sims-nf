# Bug Report: `write_gwa_to_db.R` — missing `cv_maf_effective` argument to `generate_trait_id()`

**Date:** 2026-03-21
**Nextflow version:** 24.10.4
**Executor:** local

---

## Description

Running the `test_cv_pool` profile causes the `DB_MIGRATION_WRITE_GWA_TO_DB` process to fail. The R script `write_gwa_to_db.R` calls `generate_trait_id()` with 5 positional arguments, but the function signature now requires 7 (`cv_maf_effective` and `cv_ld` were added in v=2 of the trait hash scheme). Neither argument has a default value, so R halts immediately.

The call at `modules/db_migration/write_gwa_to_db/resources/usr/bin/write_gwa_to_db.R:86` is:

```r
trait <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
```

All other DB migration scripts (`analyze_qtl.R`, `assess_sims.R`, `write_trait_data.R`) were updated to pass the two new arguments; `write_gwa_to_db.R` was missed. The module's `main.nf` also does not pass `--cv_maf_effective` or `--cv_ld` to the script.

## Command Used

```bash
NXF_VER=24.10.4 nextflow run main.nf -profile test_cv_pool,docker
```

## Error Output

```
Error in generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep,  :
  argument "cv_maf_effective" is missing, with no default
Calls: generate_trait_id -> paste0 -> sprintf
Execution halted
```

## Config / Profile

```groovy
// test_cv_pool profile — triggers a CV pool larger than the marker pool,
// exercising the non-marker causal variant code path
profiles {
    test_cv_pool {
        // (see nextflow.config for full definition)
    }
}
```

## Additional Context

- `generate_trait_id()` signature is defined in `R/database.R:479`:
  `generate_trait_id(marker_set_hash, nqtl, effect, rep, h2, cv_maf_effective, cv_ld)`
- The two new parameters (`cv_maf_effective`, `cv_ld`) were introduced when the trait hash was bumped to v=2 to encode causal-variant pool metadata.
- Fix requires: (1) adding `--cv_maf_effective` and `--cv_ld` CLI options to `write_gwa_to_db.R`, (2) passing those values through to the `generate_trait_id()` call, and (3) forwarding both arguments from `modules/db_migration/write_gwa_to_db/main.nf`.

---

*Generated with the Quarto GitHub Issue Framework. Submit via `gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
