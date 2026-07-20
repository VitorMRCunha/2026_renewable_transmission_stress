#!/usr/bin/env python3
"""Validate principal publication values in the frozen CSV package."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FAILURES = []

def check_close(label, actual, expected, tol=1e-3):
    if abs(float(actual) - float(expected)) > tol:
        FAILURES.append(f"{label}: expected {expected}, got {actual}")

def check_equal(label, actual, expected):
    if actual != expected:
        FAILURES.append(f"{label}: expected {expected}, got {actual}")

# Intact convergence
intact = ROOT / "results/frozen/intact"
rep = pd.read_csv(intact / "convergence_representative_branches.csv")
final = rep[(rep["N_requested"] == 2000) &
            (rep["from_bus"] == 26) &
            (rep["to_bus"] == 30)]
expected_lcp = {
    "S1_Low_Conc": 5.5166666667,
    "S2_High_Conc": 10.6104166667,
    "S3_Low_Dist": 2.5291666667,
    "S4_High_Dist": 6.2979166667,
}
for scenario, expected in expected_lcp.items():
    row = final[final["scenario"] == scenario]
    if row.empty:
        FAILURES.append(f"Missing final branch 26-30 row for {scenario}")
    else:
        check_close(f"{scenario} branch 26-30 LCP80", row.iloc[0]["LCP_pct"], expected, 1e-6)

inter = pd.read_csv(intact / "convergence_factorial_interaction.csv")
row = inter[inter["N_requested"] == 2000].iloc[0]
check_close("Interaction estimate", row["interaction_pp"], 1.325, 1e-9)
check_close("Interaction CI low", row["ci_low_pp"], 0.51875, 1e-9)
check_close("Interaction CI high", row["ci_high_pp"], 2.121875, 1e-9)

opf = pd.read_csv(intact / "opf_convergence_summary.csv")
check_equal("Intact successful OPFs", int(opf["successful_OPFs"].sum()), 192024)
check_equal("Intact failed OPFs", int(opf["failed_OPFs"].sum()), 0)

# Deterministic N-1
det = ROOT / "results/frozen/deterministic_n1"
summary = pd.read_csv(det / "n1_contingency_summary_all.csv")
check_equal("Deterministic outages", len(summary), 186)
check_equal("Islanding outages", int(summary["is_islanding"].sum()), 9)
check_equal("Connected outages", int((summary["is_islanding"] == 0).sum()), 177)

# Probabilistic N-1
n1 = ROOT / "results/frozen/probabilistic_n1"
n1opf = pd.read_csv(n1 / "n1_prob_opf_summary.csv")
check_equal("Probabilistic N-1 successful OPFs", int(n1opf["successful_OPFs"].sum()), 230400)
check_equal("Probabilistic N-1 failed OPFs", int(n1opf["failed_OPFs"].sum()), 0)

target = pd.read_csv(n1 / "n1_prob_target_branch_metrics.csv")
for outage, base, res, delta in [
    ("94-100", 0.0, 12.5520833333, 12.5520833333),
    ("94-96", 8.3333333333, 19.9583333333, 11.625),
]:
    q = target[(target["scenario"] == "S2_High_Conc") &
               (target["target_branch"] == "26-30") &
               (target["outage_branch"] == outage)]
    if q.empty:
        FAILURES.append(f"Missing probabilistic N-1 target row for outage {outage}")
    else:
        r = q.iloc[0]
        check_close(f"{outage} baseline LCP80", r["base_LCP80_pct"], base, 1e-6)
        check_close(f"{outage} renewable LCP80", r["LCP80_pct"], res, 1e-6)
        check_close(f"{outage} Delta LCP80", r["delta_LCP80_pp"], delta, 1e-6)

# Peak sensitivity
peak = ROOT / "results/frozen/peak_sensitivity"
popf = pd.read_csv(peak / "peak_sensitivity_opf_summary.csv")
check_equal("Peak successful OPFs", int(popf["successful_opfs"].sum()), 192120)
check_equal("Peak failed OPFs", int(popf["failed_opfs"].sum()), 0)

pt = pd.read_csv(peak / "peak_sensitivity_target_branches.csv")
expected_peak = {
    0.60: (25.0, 31.34375),
    0.65: (16.6666666667, 23.9375),
    0.70: (0.0, 10.375),
    0.75: (0.0, 5.28125),
    0.80: (0.0, 2.9791666667),
}
for factor, (base, s2) in expected_peak.items():
    q = pt[(pt["peak_frac"] == factor) &
           (pt["from_bus"] == 26) &
           (pt["to_bus"] == 30) &
           (pt["scenario"] == "S2_High_Conc")]
    if q.empty:
        FAILURES.append(f"Missing peak target row at factor {factor}")
    else:
        r = q.iloc[0]
        check_close(f"Peak {factor} base", r["base_LCP80_pct"], base, 1e-6)
        check_close(f"Peak {factor} S2", r["res_LCP80_pct"], s2, 1e-6)

if FAILURES:
    print("VALIDATION FAILED")
    for failure in FAILURES:
        print(" -", failure)
    sys.exit(1)

print("VALIDATION PASSED")
print("All principal solver counts and numerical spot checks match the frozen package.")
