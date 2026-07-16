
## 1. Did all 9 WRITE_MARKER_SET tasks complete?

```bash
nextflow log crazy_lalande -f name,status,exit,workdir | grep -i "WRITE_MARKER_SET"
```

Returns:

```log
DB_MIGRATION_WRITE_MARKER_SET (ce.hpc_0.05)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/c7/2dfe09ab38aedf6d28252b5bd9798a
DB_MIGRATION_WRITE_MARKER_SET (ct.hpc.popB_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/69/5e0dcd726e0d674a825a8879148796
DB_MIGRATION_WRITE_MARKER_SET (ce.hpc.popB_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/de/0efb84a670c78319d99e3ccf6f454d
DB_MIGRATION_WRITE_MARKER_SET (cb.hpc.popA_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/0d/bba43bc190457778add4361294bd19
DB_MIGRATION_WRITE_MARKER_SET (ct.hpc_0.05)     FAILED  135     /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/ff/269b244c25131b215209edd37bb2c4
DB_MIGRATION_WRITE_MARKER_SET (ct.hpc.popA_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/d6/25f571fc0cc36b0768439ea4121c42
DB_MIGRATION_WRITE_MARKER_SET (cb.hpc_0.05)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/dc/f535ea680930ace5be9966963bf7fb
DB_MIGRATION_WRITE_MARKER_SET (cb.hpc.popB_0.05)        FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/a6/86b67b4d0f62a628bf19ec88c4f080
DB_MIGRATION_WRITE_MARKER_SET (ce.hpc.popA_0.05)        FAILED  135     /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/f4/af069cb3e2fdf63023225f3f95c365
DB_MIGRATION_WRITE_MARKER_SET (ct.hpc_0.05)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/47/65bc22250ea8f24b57b8dae872a4d0
DB_MIGRATION_WRITE_MARKER_SET (ce.hpc.popA_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/88/7e40daf415e37e84aedc2f35399201
DB_MIGRATION_WRITE_MARKER_SET (cb.hpc.popB_0.05)        FAILED  135     /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/6b/b3bc661c045daf6540729c6393c338
DB_MIGRATION_WRITE_MARKER_SET (cb.hpc.popB_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/6a/5bab1adac6eabf93e7f0c4b8c36f16
```

##  2. How many rows are in the surviving metadata file?

```bash
singularity exec --bind /vast/eande106 --bind "$PWD" \
/vast/eande106/singularity/andersenlab-migrate-sims-20260202.img \
Rscript -e 'arrow::read_parquet("/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db/marker_set_metadata.parquet") |> dplyr::select(population, maf, marker_set_id) |> print()'
```

Returns:

```log
# A tibble: 7 × 3
  population    maf marker_set_id       
  <chr>       <dbl> <chr>               
1 ce.hpc       0.05 87748c22b16e25eac738
2 cb.hpc       0.05 468f0b478060666ce9b4
3 ct.hpc.popA  0.05 ec28451850cf520635da
4 cb.hpc.popA  0.05 9d641394742797c6b8d2
5 ct.hpc       0.05 078617e8396cb6f7c514
6 ce.hpc.popA  0.05 0b532dd017fd6c58bb6a
7 cb.hpc.popB  0.05 43010511affb5579d3f8
```

## 3. Do all 9 per-population marker files exist? (These are written to separate files and should be fine.)

```bash
ls Analysis_Results-20260326/db/marker_sets/*_markers.parquet | wc -l
```

Returns:

```log
ls: cannot access 'Analysis_Results-20260326/db/marker_sets/*_markers.parquet': No such file or directory
0
```

But 

```bash
ls Analysis_Results-20260326/db/markers/marker_sets/*_markers.parquet | wc -l
```
Returns:

```log
9
```