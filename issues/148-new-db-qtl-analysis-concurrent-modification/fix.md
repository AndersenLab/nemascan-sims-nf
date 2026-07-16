## Solutions

### Solution Applied: Simplify the DB analysis input channel creation process

The channel creation strategy in the current code requires the creation of two independent channels that are merged together, which adds unnecessary complexity and the potential for errors if a single channel is not properly synchronized.

Instead, we can directly create a single channel that emits the required inputs. These channels are sent to the `DB_MIGRATION_ANALYZE_QTL` and `R_GET_GCTA_INTERVAL` which are each part of parallele analysis paths. 

The DB_MIGRATION_ANALYZE_QTL and R_GET_GCTA_INTERVALS processes receive multiple queue-channel inputs that are paired by emission index. The current code constructs each channel independently using a wrap-combine-unwrap pattern to duplicate emissions for the BF/EIGEN threshold expansion. This has two problems:

1. ConcurrentModificationException — The unwrap (.map { it -> it[0]}) returns the same underlying ArrayList for both BF and EIGEN emissions. When both land in the same SLURM job array batch, two threads iterate the shared list simultaneously in BashWrapperBuilder.createContainerBuilder (see
bug-report.md).
2. Fragile emission-order pairing — Independent channels rely on implicit index matching. Any operator that reorders one channel but not another silently breaks the pairing.

As a fix, we can merge all upsteam channels into a single channel before the threshold expansion, perform all transforms (adding BF and EIGEN streams) on the merged channel, then split into process inputs with `multiMap` at the end. This ensures each emission is self-contained and eliminates the risk of concurrent modification. It also makes the pairing explicit and robust to reordering.

### Alternative A: Force a new ArrayList after the combine unwrap so each emission owns its own list

```groovy
 // Before (shared reference):
 .map { it -> it[0] }

 // After (independent copy):
 .map { it -> it[0].collect() }
```

Collect creates a shallow copy of the list. Which is a quick and dirty way to ensure that each emission has its own list reference and is not shared between threads.

Its probably better to get rid of the warp, unwrap anti-pattern and just create the channel with the correct structure from the start. 

### Alternative B: `.combine` with sig threshold without wrapping

```groovy
  ch_db_analysis_pheno = ch_gwa_pheno_for_analysis
      .combine(ch_db_sthresh)
      .combine(ch_db_analysis_barrier)
      .map { pheno, par, _thresh, _barrier ->
          tuple(pheno, par)
      }
```
