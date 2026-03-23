Encounted an error when testing the `feature/causal-variant-pool-tests` branch when running local end to end tests with the new `test_cv_pool` profile. The error was:

The error occured in the NF processs `DB_MIGRATION_WRITE_GWA_TO_DB` which runs the script `write_gwa_to_db.R`. The error is that the argument `cv_maf_effective` is missing when calling the function `generate_trait_id`. This argument is required for the function to work properly, and it does not have a default value, which is why the error is occurring.

Error message:

```
Error in generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep,  : 
  argument "cv_maf_effective" is missing, with no default
Calls: generate_trait_id -> paste0 -> sprintf
Execution halted
```

Command run:

```bash
NXF_VER=24.10.4 nextflow run main.nf -profile test_cv_pool,docker
```