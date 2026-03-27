# Sparse GRM Incompatibility with LOCO GWA Mapping

## Background

### What is LOCO?

LOCO (Leave-One-Chromosome-Out) is a refinement of standard mixed-model association (MLMA) that
improves power by preventing the kinship matrix from absorbing signal on the chromosome being
tested. In standard MLMA, the GRM captures genome-wide relatedness — including from the same
chromosome as the SNP under test. This creates a circularity: the kinship estimate is correlated
with the SNP being tested, which can deflate test statistics and reduce power.

LOCO's solution is to fit a separate mixed model for each chromosome using kinship estimated from
**all other chromosomes only**. For chromosome *k*, the leave-out kinship matrix is:

```
K_(-k) = K_total - (M_k / M) * K_k
```

where `K_total` is the full genome-wide GRM, `K_k` is the GRM contribution from chromosome *k*'s
SNPs (computed on-the-fly from the PLINK binary files), `M_k` is the number of markers on
chromosome *k*, and `M` is the total number of markers.

### What is a sparse GRM?

GCTA's `--make-bK-sparse` command creates a compressed approximation of the GRM (`.grm.sp`) by
zeroing out any pairwise relatedness values below a threshold (e.g., `|r_ij| < 0.05`). This
reduces storage and computation time and is appropriate for `--fastGWA-mlm-exact` (inbred mode),
which uses the sparse GRM directly as the covariance structure without any arithmetic on it.

---

## The Bug

### Root cause

The sparse GRM approximation is **incompatible with the LOCO subtraction step**.

`K_k` is computed directly from raw genotype data in the PLINK binaries — every pair gets a real,
nonzero value. `K_total_sparse` is the stored sparse GRM, which has small pairwise values set to
exactly zero. When GCTA computes `K_(-k) = K_total_sparse - (M_k / M) * K_k`, any pair where
`K_total_sparse[i,j] = 0` but `K_k[i,j] > 0` produces a **negative off-diagonal kinship value**.

Negative pairwise kinship is physically meaningless (you cannot be less related than unrelated).
A kinship matrix containing negative off-diagonal elements is not positive semi-definite. The REML
optimizer, tasked with fitting variance components under this impossible covariance structure,
cannot find a valid solution and runs off to the boundary of the parameter space.

### Where this step occurs in the mapping workflow

Within a single `GCTA_PERFORM_GWA` task running in LOCO mode, GCTA iterates over each chromosome
in sequence. For each chromosome *k* it performs two sub-steps before moving to the next:

1. **Variance component estimation (REML)** — fits a null mixed model to the phenotype using the
   leave-out kinship matrix K_(-k). The goal is to estimate what fraction of phenotypic variance
   is explained by genome-wide genetic relatedness (V(G)) vs. unexplained noise (V(e)). No SNPs
   are tested yet at this stage.

2. **Association testing** — using the V(G) and V(e) estimates from step 1, constructs the
   phenotypic variance-covariance matrix `V = V(G) * K_(-k) + V(e) * I` and inverts it. Each SNP
   on chromosome *k* is then tested by asking whether adding that SNP as a fixed effect
   significantly improves the model fit, given the genetic background captured by V. This produces
   the per-SNP p-values in the `.mlma` output file.

The REML step in (1) is a prerequisite — if it fails, association testing in (2) cannot proceed
and the entire chromosome's results are lost. This is why the error terminates the task rather
than producing partial output.

#### What the AI-REML optimizer is doing

REML (Restricted Maximum Likelihood) finds the values of V(G) and V(e) that make the observed
phenotype data most probable under the mixed model, after projecting out fixed effects (intercept
and, if applicable, the first PC as a covariate). "Restricted" means the likelihood is computed
on residuals from the fixed-effect fit rather than on raw phenotypes, which removes upward bias in
variance estimates that would otherwise occur with small sample sizes.

GCTA uses the **AI-REML** (Average Information) variant, a Newton-like algorithm where each
iteration computes an update step using the average of the observed and expected information
matrices (the matrix of second derivatives of the log-likelihood with respect to V(G) and V(e)).
This is more computationally efficient than full Newton-Raphson and converges quickly for
well-specified models.

Each iteration: compute the gradient and information matrix → solve for the update step → apply
it to the current V(G), V(e) estimates → check for convergence (change in logL < threshold).
The "information matrix is not invertible" error occurs when the information matrix becomes
singular, meaning the optimizer cannot compute the update step — the log-likelihood surface has
become so flat or pathological in the region being explored that no descent direction exists.

### Observed symptom

The REML optimizer diverged on chromosome III while chromosome II converged normally in the same
run. The outputs below are from the same GCTA process (same mapping task, same samples, same
sparse GRM), captured from `.nextflow.log`.

#### Column definitions

Each REML iteration reports four values:

| Column | Meaning |
|--------|---------|
| `Iter.` | Iteration number of the AI-REML algorithm |
| `logL` | Log-likelihood of the model at the current parameter values — should increase toward a maximum at convergence |
| `V(G)` | Estimated genetic variance — the proportion of phenotypic variance attributed to genome-wide relatedness from the leave-out chromosomes |
| `V(e)` | Estimated residual (environmental/error) variance — everything not explained by the genetic component |

`(1 component(s) constrained)` means V(G) has been forced to zero because the optimizer reached
its lower bound. In a well-specified model this is a valid outcome (some traits have no detectable
heritability); in a mis-specified model it can mean the genetic component simply cannot be fit and
all variance is being pushed into V(e).

#### Chromosome II — normal convergence

```
#Chr 2:
1655 SNPs on chromosome 2 are included in the analysis.

Performing REML analysis ... (Note: may take hours depending on sample size).
13 observations, 2 fixed effect(s), and 2 variance component(s)(including residual variance).
Prior values of variance components:  0.00000 10.13857
logL: -16.28978
Running AI-REML algorithm ...
Iter.  logL    V(G)     V(e)
1      -16.68  0.00000  12.54204   (1 component(s) constrained)
2      -17.03  0.00000  15.66727   (1 component(s) constrained)
3      -17.88  0.00000  16.64007   (1 component(s) constrained)
4      -18.13  0.00000  16.70423   (1 component(s) constrained)
5      -18.14  0.00000  16.70448   (1 component(s) constrained)
6      -18.14  0.00000  16.70448   (1 component(s) constrained)
Log-likelihood ratio converged.

Running association tests for 1655 SNPs ...
```

V(e) increases gradually from the prior (10.14) and stabilises at 16.70 by iteration 5.
logL decreases steadily toward a stable maximum. Association testing proceeds normally.

#### Chromosome III — divergence and failure

Note that GCTA carries the converged V(e) from chromosome II (16.70448) forward as the prior for
chromosome III. Despite starting from a valid prior, the optimizer immediately escapes to an
implausible value on iteration 1 and then diverges exponentially:

```
#Chr 3:
1175 SNPs on chromosome 3 are included in the analysis.

Performing REML analysis ... (Note: may take hours depending on sample size).
13 observations, 2 fixed effect(s), and 2 variance component(s)(including residual variance).
Prior values of variance components:  0.00000 16.70448
logL: -18.14486
Running AI-REML algorithm ...
Iter.  logL     V(G)     V(e)
1      -21.35   0.00000  94.15479                (1 component(s) constrained)
2      -26.52   0.00000  1391.99609              (1 component(s) constrained)
3      -41.11   0.00000  265874.10680            (1 component(s) constrained)
4      -69.98   0.00000  9598528154.52585        (1 component(s) constrained)
5      -127.70  0.00000  12509437511561904128.00000                                    (1 component(s) constrained)
6      -243.13  0.00000  21247315130549679265321089811208994816.00000                  (1 component(s) constrained)
7      -474.00  0.00000  61296632154175416873972050354089382167194763460234399866129176153559138304.00000  (1 component(s) constrained)
8      -935.74  0.00000  510154500476796989690264502143601286972537336331836277351193086180083958472323999376496350824771995757497428787879997613186320542481173765641207808.00000  (1 component(s) constrained)
Error: the information matrix is not invertible.
An error occurs, please check the options or data
```

V(e) roughly squares each iteration (94 → 1,392 → 265,874 → 9.6×10^9 → 1.25×10^19 → ...),
reaching ~10^149 before the information matrix — the matrix of second derivatives of the
log-likelihood used to compute the parameter update step — becomes numerically singular and GCTA
exits with an error. logL also diverges negatively at the same rate, indicating the model fit is
getting worse with every iteration rather than better.

The contrast between the two chromosomes is the key diagnostic signal: same run, same data, same
sparse GRM — but chromosome III's leave-out kinship matrix, after subtracting chromosome III's
contribution from the sparse total, contains sufficiently many negative off-diagonal elements to
make the covariance structure unsolvable.

---

## Empirical Evidence for This Root Cause

Three observations together support the sparse GRM subtraction mismatch as the root cause:

**1. Failures were isolated to individual chromosomes, not whole runs.**
The sparse GRM was used as the starting point for all chromosomes. If the sparse GRM were simply
incorrect or corrupt, every chromosome would fail. Chromosome-specific failures point to the
per-chromosome subtraction step — the only part of the computation that varies between
chromosomes.

**2. Inbred mode worked correctly on the same data.**
Inbred mode (`--fastGWA-mlm-exact`) uses the same sparse GRM on the same samples and produced
valid results throughout. The only difference between inbred and LOCO modes is the subtraction
step. This rules out the samples, phenotypes, and sparse GRM itself as the source of the problem
and points directly at what LOCO does differently.

**3. Switching to the full binary GRM resolved the failures completely.**
After replacing the sparse GRM with the full binary GRM (`.grm.bin`) for LOCO mode, every
chromosome across all 9 populations produced valid variance component estimates and well-formed
p-value distributions. Changing only the GRM format eliminated the failure entirely.

---

## Fix

`modules/gcta/perform_gwa/main.nf` was refactored so that the two modes use distinct GRM formats:

- **Inbred** (`--fastGWA-mlm-exact`): creates sparse GRM with `--make-bK-sparse`, passes it via
  `--grm-sparse`. Correct and efficient — no subtraction step involved.
- **LOCO** (`--mlma-loco`): skips sparse GRM creation entirely, passes the full binary GRM from
  `GCTA_MAKE_GRM` directly via `--grm`. Required for valid per-chromosome covariance decomposition.

---

## Key References

- Yang et al. (2014). "Advantages and pitfalls in the application of mixed-model association
  methods." *Nature Genetics* — introduces MLMA-LOCO and the K_(-k) decomposition.
- Yang et al. (2011). "GCTA: A Tool for Genome-wide Complex Trait Analysis." *AJHG* — describes
  GRM construction and the REML framework.
- GCTA documentation: `--mlma-loco` and `--grm-sparse` command references at
  `cnsgenomics.com/software/gcta`.

> **Note:** No published reference explicitly documents the incompatibility between `--grm-sparse`
> and `--mlma-loco`. The mechanistic explanation above is an inference from the mathematical
> structure of the LOCO algorithm, supported by the three empirical observations described.
