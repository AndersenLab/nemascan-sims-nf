# tests/pytest/test_create_causal_vars.py
#
# pytest tests for bin/create_causal_vars.py via subprocess.run().
# numpy and pandas must be installed in the active Python environment.

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parents[2]
SCRIPT = PROJECT_ROOT / "bin" / "create_causal_vars.py"
FIXTURES = PROJECT_ROOT / "tests" / "fixtures"
BIM = FIXTURES / "test.bim"

pytestmark = [
    pytest.mark.skipif(not SCRIPT.exists(), reason="create_causal_vars.py not found"),
    pytest.mark.skipif(not BIM.exists(), reason="test.bim fixture not found"),
]


def run_script(tmp_path, nqtl, effect, rep):
    """Copy PLINK binaries into tmp_path (as CV_TO_SIMS.*) and run the script.

    create_causal_vars.py calls read_bed_genotypes() unconditionally after writing
    causal_vars.txt, so .bed/.fam must be present or every call raises
    FileNotFoundError before producing any output.
    """
    for ext in (".bim", ".bed", ".fam"):
        src = BIM.with_suffix(ext)
        if src.exists():
            shutil.copy(src, tmp_path / f"CV_TO_SIMS{ext}")
    return subprocess.run(
        [sys.executable, str(SCRIPT), "CV_TO_SIMS.bim", str(nqtl), effect, str(rep)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )


def test_pool_guard_exit1_when_nqtl_exceeds_pool(tmp_path):
    n_pool = sum(1 for _ in BIM.open())
    result = run_script(tmp_path, n_pool + 1, "gamma", 1)
    assert result.returncode == 1
    msg = result.stdout + result.stderr
    assert any(kw in msg for kw in ("--cv_maf", "--cv_ld", "nqtl"))


def test_pool_guard_exit0_at_pool_boundary(tmp_path):
    n_pool = sum(1 for _ in BIM.open())
    result = run_script(tmp_path, n_pool, "gamma", 1)
    assert result.returncode == 0


def test_rng_determinism_same_rep(tmp_path):
    tmp1 = tmp_path / "run1"; tmp1.mkdir()
    tmp2 = tmp_path / "run2"; tmp2.mkdir()
    r1 = run_script(tmp1, 3, "gamma", 7)
    r2 = run_script(tmp2, 3, "gamma", 7)
    assert r1.returncode == 0
    assert r2.returncode == 0
    assert (tmp1 / "causal_vars.txt").read_text() == (tmp2 / "causal_vars.txt").read_text()


def test_rng_advance_different_reps(tmp_path):
    tmp1 = tmp_path / "run1"; tmp1.mkdir()
    tmp2 = tmp_path / "run2"; tmp2.mkdir()
    r1 = run_script(tmp1, 5, "gamma", 1)
    r2 = run_script(tmp2, 5, "gamma", 2)
    assert r1.returncode == 0
    assert r2.returncode == 0
    assert (tmp1 / "causal_vars.txt").read_text() != (tmp2 / "causal_vars.txt").read_text()


def test_output_format_two_columns_no_header(tmp_path):
    result = run_script(tmp_path, 3, "gamma", 1)
    assert result.returncode == 0
    lines = (tmp_path / "causal_vars.txt").read_text().splitlines()
    assert len(lines) == 3
    cols = lines[0].split()
    assert len(cols) == 2
    float(cols[1])  # raises ValueError if not numeric


def test_uniform_effect_range(tmp_path):
    result = run_script(tmp_path, 5, "0.1-0.5", 1)
    assert result.returncode == 0
    lines = (tmp_path / "causal_vars.txt").read_text().splitlines()
    effects = [abs(float(line.split()[1])) for line in lines]
    assert all(0.1 <= e <= 0.5 for e in effects)


@pytest.mark.skipif(
    not (BIM.with_suffix(".bed").exists() and BIM.with_suffix(".fam").exists()),
    reason="test.bed / test.fam not available for genotype extraction",
)
def test_causal_genotypes_tsv_written(tmp_path):
    import pandas as pd

    result = run_script(tmp_path, 3, "gamma", 1)
    assert result.returncode == 0
    gfiles = list(tmp_path.glob("causal_genotypes.sim.3.1.tsv"))
    assert len(gfiles) == 1
    g = pd.read_csv(gfiles[0], sep="\t")
    assert {"QTL", "strain", "allele"}.issubset(g.columns)
    assert len(g) > 0
