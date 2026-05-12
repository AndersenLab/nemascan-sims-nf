# Implementation plan: pre-transfer manifest + archival workflow

See [`design.md`](design.md) for the architecture and rationale. This
file is the checklist.

## Deliverables

- [x] `pretransfer_manifest.sh` at the repo root (~220 lines of bash, GNU
  tar + zstd + sha256sum).
- [x] `issues/pretransfer-manifest-archive/` folder with `feature_request.md`,
  `design.md`, and this `plan.md`.
- [ ] `docs/transfer-and-archival.qmd` — user-facing Quarto usage page,
  modeled on `docs/testing-hpc.qmd` with a mermaid flowchart and three
  phase sections.
- [ ] `docs/_quarto.yml` — add the new page to the Development submenu.
- [ ] `2026-04-03-production_three_species_sims.qmd` — rewrite Phase 8
  into 8a (stage) / 8b (transfer) / 8c (verify + extract), add a
  `callout-warning` in Phase 9 guarding the `rm -rf /scratch4/...`
  cleanup on Phase 8c success, and expand the top-of-file checklist.

## Verification

Local-runnable:

- [ ] **Unit smoke test.** Build a fake run directory with a `db/` subdir
  and a stub TSV, run `bash pretransfer_manifest.sh /tmp/fake_run
  /tmp/staging`. Confirm all four sidecars appear. `sha256sum -c` against
  both hash files exits 0.
- [ ] **Reproducibility check.** Run the script twice on unchanged source.
  `sha256sum` the two `.tar.zst` outputs; they must be identical.
- [ ] **Idempotency check.** Re-run without `--force` after source is
  untouched; script should print "staging is current" and exit 0 without
  rewriting. Then re-run with `--force` and confirm regeneration.
- [ ] **Failure injection.** Corrupt one byte of the archive after
  staging; confirm `sha256sum -c <RUN>.tar.zst.sha256` exits non-zero
  and the script-printed verify step catches it before extraction.

Requires HPC access:

- [ ] **Round-trip on real HPC data.** On Rockfish, against a completed
  `test_hpc`-profile output (~1 GB): stage → rsync to local → extract →
  `sha256sum -c` the inner manifest. Every file reports `OK`.

Requires Quarto toolchain:

- [ ] **Docs render.** `quarto render docs/` — confirm the new page
  renders, the navbar entry shows up, and the mermaid flowchart draws.
- [ ] **Production notes render.** `quarto render
  2026-04-03-production_three_species_sims.qmd` — confirm no broken
  syntax after the Phase 8/9 rewrite. (The file is `format: gfm`, not
  html, so the mermaid in Quarto-html doesn't apply here — just check
  that the new `callout-warning` block is valid.)

## Rollback

All changes are additive except the `2026-04-03-production_three_species_sims.qmd`
edit (which modifies existing Phase 8/9 content). To roll back fully:

1. `rm pretransfer_manifest.sh`
2. `rm -rf issues/pretransfer-manifest-archive/`
3. `rm docs/transfer-and-archival.qmd`
4. Revert the `docs/_quarto.yml` navbar edit.
5. Revert the `2026-04-03-production_three_species_sims.qmd` Phase 8/9
   rewrite.

No data files, Parquet schemas, Nextflow modules, or containers are
touched, so there is no migration or database concern on rollback.
