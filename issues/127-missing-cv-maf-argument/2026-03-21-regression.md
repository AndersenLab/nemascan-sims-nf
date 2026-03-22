## Update — `2026-03-21`

**Phase:** RCA
**Status:** Resolved — `c1fcfe0`

---

## Root Cause Analysis

### Summary

The Option A fix landed in commit `067413f` introduced a new regression: the `.map`
closure that constructs `ch_gwa_db_inputs` was written against an outdated view of
`ch_marker_set_params_for_gwa`, which does not account for the `strains` and
`strainfile` elements added by the strainfile-provenance feature. Groovy receives a
15-element `LinkedList` but the closure declares 13 parameters; the arity mismatch
causes a fatal "Invalid method invocation `call`" error before `write_gwa_to_db.R` ever
runs.

### Diagnostic Steps

1. Reproduced the error by running
   `NXF_VER=24.10.4 nextflow run main.nf -profile test_cv_pool,docker` after
   applying commit `067413f`.
2. Read the Nextflow error report — the failing closure is `_closure52` in
   `main.nf`, not in an R script, ruling out any R-level argument issue.
3. Enumerated the elements in the LinkedList reported by the error:

   | Position | Value (representative) | Field |
   |---|---|---|
   | 0 | `ce.test.200strains` | group |
   | 1 | `0.05` | maf |
   | 2 | `5` | nqtl |
   | 3 | `gamma` | effect |
   | 4 | `1` | rep |
   | 5 | `0.8` | h2 |
   | 6 | `inbred` | mode |
   | 7 | `fastGWA` | suffix |
   | 8 | `pca` | type |
   | 9 | `c_elegans` | species |
   | 10 | `20210901` | vcf_release_id |
   | 11 | `0.8` | ms_ld |
   | 12 | `AB1,BRC20067,...` | strains (comma-separated) |
   | 13 | `/…/test_strains.txt` | strainfile path |
   | 14 | `0.01` | cv_maf_eff |

   15 elements, not 13 — confirming that two unexpected elements appear between
   `ms_ld` (position 11) and `cv_maf_eff` (position 14).

4. Cross-referenced `ch_sf.marker_set_params` definition (`main.nf:213`):

   ```groovy
   marker_set_params: [meta.id, ms_maf, species, extractVcfReleaseId(vcf), ms_ld, strains, strainfile]
   ```

   This 7-element tuple was extended by the strainfile-provenance feature to carry
   `strains` (index 5) and `strainfile` (index 6). All three downstream taps
   (`ch_marker_set_params_for_ms`, `ch_marker_set_params_for_gm`,
   `ch_marker_set_params_for_gwa`) receive this 7-element tuple.

5. Traced `ch_gwa_db_inputs` construction (`main.nf:634-641`):

   - `.combine(ch_marker_set_params_for_gwa, by: [0, 1])` joins on `(group, maf)`;
     remaining non-key elements from `ch_marker_set_params_for_gwa` are
     `species, vcf_release_id, ms_ld, strains, strainfile` (5 fields, not 3).
   - `.combine(ch_cv_maf_keyed_for_gwa_write, by: 0)` appends `cv_maf_eff`.
   - **Actual result:** 15-element tuple.

6. Compared against the closure declaration at `main.nf:637-638`:

   ```groovy
   .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
          species, vcf_release_id, ms_ld, cv_maf_eff ->
   ```

   13 parameters. Groovy cannot spread a 15-element list into 13 parameters; it
   passes the entire list as a single argument → `(LinkedList)` invocation failure.

7. Confirmed the fix document (`2026-03-21-rca-fix.md`) was written before the
   strainfile-provenance feature was merged. Its schema comment
   (`// Result: tuple(…, ms_ld)`) omits `strains` and `strainfile`, confirming
   the fix was authored against the old 5-element non-key assumption.

### Key Evidence

**Nextflow error output:**

```
ERROR ~ Invalid method invocation `call` with arguments: [ce.test.200strains, 0.05, 5,
gamma, 1, 0.8, inbred, fastGWA, pca, c_elegans, 20210901, 0.8,
AB1,BRC20067,...,XZ2018, /Users/ryanmckeown/Desktop/nemascan-sims-nf/data/test/test_strains.txt,
0.01] (java.util.LinkedList) on _closure52 type
```

**Nextflow error report:**

```
No signature of method: Script_a71711320ed3590c$_runScript_closure1$_closure4$_closure52.call()
is applicable for argument types: (LinkedList) values: [[ce.test.200strains, …]]
Possible solutions: any(), any(), any(groovy.lang.Closure), each(groovy.lang.Closure), …
```

**`main.nf:213` — `marker_set_params` schema (current):**

```groovy
marker_set_params: [meta.id, ms_maf, species, extractVcfReleaseId(vcf), ms_ld, strains, strainfile]
//                   [0]      [1]    [2]       [3]                       [4]    [5]      [6]
```

**`main.nf:634-641` — `ch_gwa_db_inputs` construction (current, broken):**

```groovy
ch_gwa_db_inputs = ch_db_params
    .combine(ch_marker_set_params_for_gwa, by: [0, 1])
    // Actual result: (group, maf, nqtl, effect, rep, h2, mode, suffix, type,
    //                 species, vcf_release_id, ms_ld, strains, strainfile)  ← 14 elements
    .combine(ch_cv_maf_keyed_for_gwa_write, by: 0)
    // Actual result: (group, maf, nqtl, effect, rep, h2, mode, suffix, type,
    //                 species, vcf_release_id, ms_ld, strains, strainfile, cv_maf_eff) ← 15 elements
    .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
           species, vcf_release_id, ms_ld, cv_maf_eff ->        // ← 13 params (missing strains, strainfile)
        tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
              species, vcf_release_id, ms_ld, cv_maf_eff, cv_ld)
    }
```

### Confirmed Root Cause

The `.map` closure added by commit `067413f` declares 13 parameters but receives a
15-element tuple because the fix was written against the pre-strainfile-provenance
schema of `ch_marker_set_params_for_gwa`. The two extra elements — `strains`
(comma-separated strain list) and `strainfile` (path to the strainfile) — are present
in positions 12–13 of the combined tuple but absent from the closure parameter list.
Groovy cannot destructure a `LinkedList` of 15 into 13 named parameters and raises an
"Invalid method invocation" error on `_closure52`.

---

## Solution

### Approach

Add `strains` and `strainfile` (prefixed with `_` to signal they are unused by
`write_gwa_to_db`) to the `.map` closure parameter list so that Groovy can destructure
the 15-element tuple. Drop them from the emitted tuple — `DB_MIGRATION_WRITE_GWA_TO_DB`
does not need strain provenance data. Update the schema comment to reflect the full
15-element intermediate and the 14-element output.

Only `main.nf` requires a change; the module and R script are correct as committed.

### Changes Made

**`main.nf` — fix `.map` closure parameter list and schema comment (lines 633–642):**

```groovy
// before
    // combine(by:0) adds cv_maf_eff keyed by group; cv_ld is a scalar folded in via map
    ch_gwa_db_inputs = ch_db_params
        .combine(ch_marker_set_params_for_gwa, by: [0, 1])
        .combine(ch_cv_maf_keyed_for_gwa_write, by: 0)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
               species, vcf_release_id, ms_ld, cv_maf_eff ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                  species, vcf_release_id, ms_ld, cv_maf_eff, cv_ld)
        }
    // Result: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species, vcf_release_id, ms_ld, cv_maf_effective, cv_ld)

// after
    // combine(by:[0,1]): ch_marker_set_params_for_gwa has 7 elements
    //   [group, maf, species, vcf_release_id, ms_ld, strains, strainfile]
    //   → intermediate: (group, maf, nqtl, effect, rep, h2, mode, suffix, type,
    //                    species, vcf_release_id, ms_ld, strains, strainfile)  14 elements
    // combine(by:0) appends cv_maf_eff; cv_ld is a pipeline scalar folded in via map
    //   → intermediate: (…, strains, strainfile, cv_maf_eff)  15 elements
    // map: drop strains/strainfile (not needed by write_gwa_to_db)
    ch_gwa_db_inputs = ch_db_params
        .combine(ch_marker_set_params_for_gwa, by: [0, 1])
        .combine(ch_cv_maf_keyed_for_gwa_write, by: 0)
        .map { group, maf, nqtl, effect, rep, h2, mode, suffix, type,
               species, vcf_release_id, ms_ld, _strains, _strainfile, cv_maf_eff ->
            tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type,
                  species, vcf_release_id, ms_ld, cv_maf_eff, cv_ld)
        }
    // Result: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type, species, vcf_release_id, ms_ld, cv_maf_effective, cv_ld)
```

### Rationale

`write_gwa_to_db.R` writes GWA statistics and constructs the trait hash; it does not
write marker-set provenance. Strain list and strainfile path are only needed by
`write_marker_set.R`. Dropping `_strains` and `_strainfile` in the map keeps the
downstream module input tuple compact and consistent with the comment in
`2026-03-21-rca-fix.md`. Using `_` prefixes on the unused parameters follows the
existing convention in `main.nf` (e.g., `_barrier`).

---

## Validation

### Test Run

**Command:**

```bash
NXF_VER=24.10.4 nextflow run main.nf -profile test_cv_pool,docker
```

**Executor:** local
**Resume from cache:** no (clean run)

### Pass / Fail

Fix applied in commit `c1fcfe0`. Full local pipeline (`-profile test_cv_pool,docker`) ran end-to-end and completed successfully.

- [x] Pipeline completed without error
- [x] `db_simulation_assessment_results.tsv` present and non-empty
- [x] Trait IDs in `mappings_metadata.parquet` match those in `traits/` (v=2 hash, 20-char hex)
- [x] No regression on `test` profile (standard marker-pool run)
- [x] `Rscript tests/run_tests.R` passes

### Remaining Concerns

- The comment at `main.nf:578` (`// Result: tuple(group, maf, genotype_matrix, species,
  vcf_release_id, ms_ld)`) is also outdated — it too omits `strains` and `strainfile`.
  It should be updated to reflect the 8-element actual tuple, though that path has no
  arity mismatch (the module input consumes all elements).
- Any future extension of `marker_set_params` (e.g., adding a new provenance field)
  will require updating all `.map` closures that consume `ch_marker_set_params_for_gwa`,
  `ch_marker_set_params_for_ms`, and `ch_marker_set_params_for_gm`. Consider adding a
  schema comment at the multiMap site (`main.nf:213`) and at each fan-out `.map` to
  make arity drift visible.

---

*Comment generated from `comment.qmd` template · nemascan-sims-nf issue workflow*
