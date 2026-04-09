# 2026-04-08: Test low VP alert
Create test parameters to trigger low VP alert

```bash
echo "1" > tmp/test_nqtl.csv
echo -e "0.05\n0.1\n0.2" > tmp/test_h2.csv
echo "gamma" > tmp/test_effect.csv
```

Run the pipeline

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile docker \
  --strainfile data/test/test_strains.txt \
  --nqtl_file tmp/test_nqtl.csv \
  --h2_file tmp/test_h2.csv \
  --effect_file tmp/test_effect.csv \
  --reps 50
```

Use the task hashes to check the logs for the low VP alert.

Nextflow work dirs are two levels deep: `work/XX/YYYYYY.../`. The log
shows truncated hashes like `[ab/cdef12]`, so a glob suffix (`*`) is
needed to resolve the full directory path.

```bash
#!/usr/bin/env bash
# ── Per-task REML loop output (hash-based) ──────────────────────────
# Extract GCTA_MAKE_GRM task hashes from the Nextflow log and print
# every REML-loop message from each task's .command.out.
grep 'GCTA_MAKE_GRM' .nextflow.log | grep 'Submitted' | \
  sed -n 's/.*\[\([a-f0-9]*\/[a-f0-9]*\)\].*/\1/p' | \
  while read hash; do
    dir=(work/${hash}*)
    [[ -d "${dir[0]}" ]] || continue
    echo "=== ${dir[0]} ==="
    grep -E "Vp=|scaling x1000|WARNING|REML failed" "${dir[0]}/.command.out" 2>/dev/null
  done

# ── Summary across all work directories ─────────────────────────────
# Patterns match the exact echo statements in modules/gcta/make_grm/main.nf

echo ""
echo "=== Traits that needed scaling ==="
grep -rl "scaling x1000" work/*/*/.command.out 2>/dev/null | \
  while read f; do grep "scaling x1000" "$f"; done

echo ""
echo "=== Traits that hit max rounds ==="
grep -rl "WARNING: Vp=" work/*/*/.command.out 2>/dev/null | \
  while read f; do grep "WARNING:" "$f"; done

echo ""
echo "=== Traits where REML failed (fallback Vp=1.0) ==="
grep -rl "REML failed" work/*/*/.command.out 2>/dev/null | \
  while read f; do grep "REML failed" "$f"; done

echo ""
echo "=== Traits that passed on round 1 (no scaling needed) ==="
grep -rl "at round 1; done" work/*/*/.command.out 2>/dev/null | \
  while read f; do grep "at round 1; done" "$f"; done

# ── Counts ──────────────────────────────────────────────────────────
echo ""
echo "=== Summary counts ==="
total=$(grep -rl "Vp=" work/*/*/.command.out 2>/dev/null | wc -l | tr -d ' ')
passed_r1=$(grep -rl "at round 1; done" work/*/*/.command.out 2>/dev/null | wc -l | tr -d ' ')
scaled=$(grep -rl "scaling x1000" work/*/*/.command.out 2>/dev/null | wc -l | tr -d ' ')
max_rounds=$(grep -rl "WARNING: Vp=" work/*/*/.command.out 2>/dev/null | wc -l | tr -d ' ')
reml_fail=$(grep -rl "REML failed" work/*/*/.command.out 2>/dev/null | wc -l | tr -d ' ')
echo "Total REML tasks:        ${total}"
echo "Passed round 1:          ${passed_r1}"
echo "Needed scaling:          ${scaled}"
echo "Hit max rounds (4):      ${max_rounds}"
echo "REML failed (fallback):  ${reml_fail}"
```
