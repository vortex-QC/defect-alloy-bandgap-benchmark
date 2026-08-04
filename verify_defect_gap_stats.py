#!/usr/bin/env python3
"""Defect/Alloy band-gap endpoint-interpolation benchmark - reproduction script.

Reproduces all numbers in:
Qin, C. (2026). Systematic Underestimation of Defect and Alloy Band Gaps by DFT
Endpoint Interpolation. Zenodo: 10.5281/zenodo.21788464 (v1.1).

Inputs (in this directory):
  defect_gap_pipeline_v2.json  - 274 matched + 275 unmatched analysis set
Outputs: prints all reported statistics.
"""
import json
import re
import statistics

def main():
    with open("defect_gap_pipeline_v2.json") as f:
        data = json.load(f)
    results = data["results"]
    unmatched = data["unmatched"]

    devs = [r["dev"] for r in results]
    print(f"Matched systems: {len(results)}/{data['n_total']} ({100*len(results)/data['n_total']:.0f}%)")
    print(f"Types: {dict(__import__('collections').Counter(r['type'] for r in results))}")
    print(f"Mean deviation:  {statistics.mean(devs):+.3f} eV")
    print(f"Median deviation:{statistics.median(devs):+.3f} eV")
    print(f"MAE:             {statistics.mean(abs(d) for d in devs):.3f} eV")
    und = 100 * sum(1 for d in devs if d > 0) / len(devs)
    print(f"Underestimation rate: {und:.0f}%")

    print("\nBand-gap dependence:")
    for lo, hi in [(0, .5), (.5, 1), (1, 2), (2, 4)]:
        v = [r["dev"] for r in results if lo <= r["expt"] < hi]
        if v:
            rate = 100 * sum(1 for d in v if d > 0) / len(v)
            print(f"  [{lo},{hi}) eV: n={len(v)} mean={statistics.mean(v):+.3f} under={rate:.0f}%")

    # Endpoint-quality classes
    PHYS_ZERO = {"HgTe", "BiTe", "BiSe", "PbTe", "SnTe", "HgSe", "HgS"}
    def parse_endpoints(ep):
        return [(f, float(g)) for f, g in re.findall(r"([A-Za-z0-9]+)\(([\d.]+)\)", ep)]
    phys, calc = [], []
    for r in results:
        eps = parse_endpoints(r["endpoints"])
        if len(eps) < 2:
            continue
        zero = [f for f, g in eps if g < 0.01]
        if any(f in PHYS_ZERO for f in zero):
            phys.append(r["dev"])
        elif any(g < 0.01 for f, g in eps):
            calc.append(r["dev"])
    print(f"\nPhysically gapless endpoints:  n={len(phys)} mean={statistics.mean(phys):+.3f}")
    print(f"Computationally gapless:      n={len(calc)} mean={statistics.mean(calc):+.3f}")

    # Correction model
    import numpy as np
    rows = []
    for r in results:
        eps = parse_endpoints(r["endpoints"])
        if len(eps) < 2:
            continue
        zero = [f for f, g in eps if g < 0.01]
        rows.append({"bg": r["expt"], "dev": r["dev"],
                     "phys": 1 if any(f in PHYS_ZERO for f in zero) else 0,
                     "calc": 1 if any(g < 0.01 for f, g in eps) and not any(f in PHYS_ZERO for f in zero) else 0})
    X = np.array([[r["bg"], r["phys"], r["calc"], 1.0] for r in rows])
    Y = np.array([r["dev"] for r in rows])
    coef = np.linalg.lstsq(X, Y, rcond=None)[0]
    pred = X @ coef
    mae = np.mean(np.abs(pred - Y))
    r2 = 1 - np.sum((Y - pred) ** 2) / np.sum((Y - np.mean(Y)) ** 2)
    print(f"\nCorrection model: dev = {coef[0]:.3f}*gap + {coef[1]:+.3f}*phys + {coef[2]:+.3f}*calc {coef[3]:+.3f}")
    print(f"  MAE={mae:.3f} eV  R2={r2:.3f}")

    print(f"\nUnmatched systems: {len(unmatched)} (forward predictions in paper)")

if __name__ == "__main__":
    main()
