
## Summary 

Failure during `DB_MIGRATION_ANALYZE_QTL` due to I/O errror when trying to read mappings file from parquet DB. 

Failure impacted many processes with 2024 failures of the DB_MIGRATION_ANALYZE_QTL module.

All upstream processes completed successfully so the issue is doesn't impact the database content, but it does prevent the analysis from completing.

## Error message

```
2026-04-09 01:26:37 [INFO] Analyzing mapping: 2e59b3684db2b4a2e256 with threshold: BF      
  Error: rapi_prepare: Failed to prepare query   CREATE VIEW mappings AS                     
    SELECT * FROM read_parquet('/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sim
s-nf/Sims_ce-cb-ct_nqtl1-50_h2grid_50reps/db/mappings/**/data.parquet', hive_partitioning = t
rue, union_by_name = true)                                                                   
  Error: IO Error: Cannot open file "/vast/eande106/projects/Ryan/simulation_pipeline/nemasca
n-sims-nf/Sims_ce-cb-ct_nqtl1-50_h2grid_50reps/db/mappings/population=ce.full/mapping_id=4012
d3f36a69c83d744a/data.parquet": Too many open files in system                                
  In addition: Warning messages:                                                             
  1: In normalizePath("~") :
    path[1]="/home/rmckeow1": No such file or directory                                     
  2: In normalizePath("~") :
    path[1]="/home/rmckeow1": No such file or directory                                     
  Execution halted
  Warning message:
  Connection is garbage-collected, use dbDisconnect() to avoid this.                        
```
