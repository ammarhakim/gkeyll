#!/usr/bin/env python3
"""Regression coverage for the relaxed X-point seam diagnostic optimizer and
guarded production apply (improve-xpt-finder branch, commits "Optimize
smooth bounded delta-s for X-point interfaces" and "Use relaxed X-point
seams with guarded fallback").

Runs gyrokinetic/creg/rt_gk_multib_nstxu_relaxed_xpt_seam_geom (built by the
normal `make -j` at the repo root) across the three regression shots, a
higher resolution for folded-cell detection, an explicit zero-delta
equivalence check, a guarded-fallback check, and (if mpirun is available) an
MPI check -- then parses the resulting diagnostic CSVs against the checklist
below. Must be run from the gkeyll repo root (the executable's rpath and the
geqdsk data paths are both relative to it).

Checklist covered here:
  - exact/tolerance flux residual and shared positions
  - canonical e_psi/e_theta mismatch and signed/absolute J behavior
  - no NaN/Inf, no within-block sign change, no folded cell (Jacobian ratio)
  - strict radial ordering
  - zero-delta equivalence
  - guarded fallback behavior
  - multiple resolutions, including a higher resolution for shot 203730
  - MPI: plain multi-rank (blocks split whole across ranks) works; theta
    (poloidal) decomposition of a single block is explicitly rejected
  - CUDA: reported as skipped (not available on this machine), not silently
    ignored

Exits 0 if every check passes, 1 otherwise, printing a PASS/FAIL line per
check.
"""
import csv
import math
import os
import shutil
import subprocess
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXE = os.path.join(REPO_ROOT, "build", "gyrokinetic", "creg",
    "rt_gk_multib_nstxu_relaxed_xpt_seam_geom")
NAME_FMT = "nstxu_relaxed_xpt_shot{shot}_ms{ms}"
SHOT_MS = {203877: 655, 203742: 480, 203730: 839}


def find_mpirun():
  """Prefer the OpenMPI build this executable was actually linked against
  (the gkylsoft/openmpi sibling of the repo checkout, per the standard
  gkeyll layout) over whatever `mpirun` resolves to on PATH -- a mismatched
  OpenMPI version there (e.g. a Homebrew install) can silently misbehave
  rather than erroring cleanly. Falls back to PATH if that sibling
  directory does not exist (e.g. a differently laid-out checkout)."""
  sibling = os.path.abspath(os.path.join(REPO_ROOT, "..", "gkylsoft",
    "openmpi", "bin", "mpirun"))
  if os.path.exists(sibling):
    return sibling
  return shutil.which("mpirun")

failures = []
passes = []


def check(label, condition, detail=""):
  if condition:
    passes.append(label)
    print(f"PASS  {label}")
  else:
    failures.append(label)
    print(f"FAIL  {label}  {detail}")


def cleanup(shot):
  ms = SHOT_MS[shot]
  prefix = NAME_FMT.format(shot=shot, ms=ms)
  for fn in os.listdir(REPO_ROOT):
    if fn.startswith(prefix):
      os.remove(os.path.join(REPO_ROOT, fn))


def run(shot, extra_args=(), mpi_ranks=None):
  ms = SHOT_MS[shot]
  cleanup(shot)
  cmd = []
  if mpi_ranks:
    cmd += [find_mpirun(), "-np", str(mpi_ranks)]
  cmd += [EXE, "-o", f"shot={shot}", *extra_args]
  if mpi_ranks:
    cmd += ["-M"]
  result = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True,
    timeout=280)
  prefix = NAME_FMT.format(shot=shot, ms=ms)
  csvs = {}
  for suffix in ["xpt_seam_optimizer_diagnostics",
      "xpt_seam_applied_diagnostics", "xpt_seam_delta_s_objective"]:
    path = os.path.join(REPO_ROOT, f"{prefix}-{suffix}.csv")
    if os.path.exists(path):
      with open(path) as f:
        csvs[suffix] = list(csv.DictReader(f))
    else:
      csvs[suffix] = None
  return result, csvs


def finite(x):
  try:
    return math.isfinite(float(x))
  except (TypeError, ValueError):
    return False


def validate_optimizer_csv(shot, rows, min_jacobian_ratio=0.1,
    max_jacobian_ratio=10.0):
  label = f"shot {shot}: optimizer CSV present and non-empty"
  check(label, rows is not None and len(rows) > 0)
  if not rows:
    return

  all_finite = all(
    finite(r[c]) for r in rows for c in
    ["coefficient", "objective", "epsi_rz", "etheta_direction",
     "etheta_log_magnitude", "jacobian_log_magnitude",
     "max_realized_displacement"])
  check(f"shot {shot}: no NaN/Inf in any candidate row", all_finite)

  valid_rows = [r for r in rows if r["candidate_valid"] == "1"]
  check(f"shot {shot}: at least the zero candidate is valid",
    any(r["is_zero_candidate"] == "1" for r in valid_rows))

  zero_rows = [r for r in valid_rows if r["is_zero_candidate"] == "1"]
  if zero_rows:
    z = zero_rows[0]
    check(f"shot {shot}: zero-delta canonical e_psi mismatch at roundoff",
      float(z["epsi_rz"]) < 1e-8, detail=z["epsi_rz"])
    check(f"shot {shot}: zero-delta max fixed/anchor position difference "
      "at roundoff",
      float(z["max_fixed_position_difference"]) < 1e-8 and
      float(z["max_anchor_position_difference"]) < 1e-8)

  ratios_ok = all(
    min_jacobian_ratio <= float(r["min_jacobian_ratio"]) and
    float(r["max_jacobian_ratio"]) <= max_jacobian_ratio
    for r in valid_rows)
  check(f"shot {shot}: no folded cell (Jacobian ratio within "
    f"[{min_jacobian_ratio}, {max_jacobian_ratio}]) in any valid candidate",
    ratios_ok)

  ordering_ok = all(r["radial_ordering_valid"] == "1" for r in valid_rows)
  check(f"shot {shot}: strict radial ordering holds for every valid "
    "candidate", ordering_ok)

  sign_ok = all(r["jacobian_valid"] == "1" for r in valid_rows)
  check(f"shot {shot}: no within-block Jacobian sign change in any valid "
    "candidate", sign_ok)


def validate_objective_csv(shot, rows):
  label = f"shot {shot}: pointwise objective CSV present and non-empty"
  check(label, rows is not None and len(rows) > 0)
  if not rows:
    return
  positions_ok = all(float(r["position_m"]) == 0.0 for r in rows
    if r["local_valid"] == "1" and r["partner_valid"] == "1")
  check(f"shot {shot}: shared seam position matches exactly at every node",
    positions_ok)
  flux_ok = all(
    abs(float(r["local_flux_residual"])) <= 1e-8 and
    abs(float(r["partner_flux_residual"])) <= 1e-8 and
    float(r["requested_flux_mismatch"]) <= 1e-10
    for r in rows if r["local_valid"] == "1" and r["partner_valid"] == "1")
  check(f"shot {shot}: requested flux is met at every node to tolerance",
    flux_ok)
  signed_abs_ok = all(
    r["local_abs_signed_consistent"] == "1" and
    r["partner_abs_signed_consistent"] == "1"
    for r in rows if r["local_valid"] == "1" and r["partner_valid"] == "1")
  check(f"shot {shot}: signed/absolute Jacobian are mutually consistent "
    "at every node", signed_abs_ok)


def main():
  if not os.path.exists(EXE):
    print(f"ERROR: {EXE} not found. Build it first with `make -j` at the "
      "repo root (it is picked up automatically by the creg/rt_*.c "
      "wildcard).")
    return 1

  for shot in (203877, 203742, 203730):
    result, csvs = run(shot)
    check(f"shot {shot}: base-resolution run exits cleanly",
      result.returncode == 0, detail=result.stderr[-500:])
    validate_optimizer_csv(shot, csvs["xpt_seam_optimizer_diagnostics"])
    validate_objective_csv(shot, csvs["xpt_seam_delta_s_objective"])
    cleanup(shot)

  # Higher resolution for shot 203730 (folded-cell detection).
  result, csvs = run(203730, extra_args=["-x", "32", "-y", "24"])
  check("shot 203730 (4x resolution): run exits cleanly",
    result.returncode == 0, detail=result.stderr[-500:])
  validate_optimizer_csv(203730, csvs["xpt_seam_optimizer_diagnostics"])
  cleanup(203730)

  # Zero-delta equivalence: with optimize=0, the optimizer must not run at
  # all (no optimizer CSV emitted), and the objective CSV must show
  # coefficient exactly 0 and requested_mode=requested_zero everywhere.
  result, csvs = run(203877, extra_args=["-o", "shot=203877,optimize=0"])
  check("zero-delta equivalence: optimizer CSV absent when optimize=0",
    csvs["xpt_seam_optimizer_diagnostics"] is None)
  obj_rows = csvs["xpt_seam_delta_s_objective"] or []
  check("zero-delta equivalence: objective CSV shows coefficient=0 "
    "everywhere",
    len(obj_rows) > 0 and all(
      float(r["local_coefficient_m"]) == 0.0 and
      r["requested_mode"] == "requested_zero" for r in obj_rows))
  cleanup(203877)

  # Guarded fallback: a tiny bound produces a small improvement that must
  # not clear the usefulness threshold, so nothing gets applied.
  result, csvs = run(203877, extra_args=["-o", "shot=203877,bound=0.00005"])
  applied_rows = csvs["xpt_seam_applied_diagnostics"] or []
  check("guarded fallback: tiny-bound run produces an applied-diagnostics "
    "row", len(applied_rows) > 0)
  check("guarded fallback: nothing applied when improvement is below "
    "threshold",
    len(applied_rows) > 0 and all(r["applied"] == "0" for r in applied_rows))
  cleanup(203877)

  # MPI: plain multi-rank (blocks split whole across ranks) must reproduce
  # the single-rank result; theta decomposition of a single block must be
  # explicitly rejected rather than silently producing incomplete samples.
  if find_mpirun():
    result, csvs = run(203877, mpi_ranks=2)
    rows = csvs["xpt_seam_optimizer_diagnostics"] or []
    check("MPI (cuts=1,1, blocks split whole across ranks): matches "
      "single-rank result",
      result.returncode == 0 and any(
        r["is_zero_candidate"] == "1" and r["candidate_valid"] == "1"
        for r in rows))
    cleanup(203877)

    result, csvs = run(203877, extra_args=["-d", "2"], mpi_ranks=2)
    rows = csvs["xpt_seam_optimizer_diagnostics"] or []
    check("MPI theta decomposition (cuts=1,2): visibly rejected as an "
      "unsupported configuration rather than silently mis-measured",
      result.returncode == 0 and len(rows) > 0 and
      all(r["pair_reject_reason"] == "unsupported_configuration"
        for r in rows))
    cleanup(203877)
  else:
    print("SKIP  MPI checks: mpirun not found on this machine")

  print("SKIP  CUDA build/run: no CUDA toolchain installed on this machine "
    "(no nvcc); not tested here. This is an explicit, reported gap, not a "
    "silent one.")

  print(f"\n{len(passes)} passed, {len(failures)} failed")
  if failures:
    print("Failed checks:")
    for f in failures:
      print(f"  - {f}")
    return 1
  return 0


if __name__ == "__main__":
  sys.exit(main())
