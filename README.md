# Defect/Alloy Band-Gap Endpoint-Interpolation Benchmark

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21788464.svg)](https://doi.org/10.5281/zenodo.21788464)

**Systematic Underestimation of Defect and Alloy Band Gaps by DFT Endpoint Interpolation: A 549-Entry Experimental Benchmark with a Hypothesis-Driven Correction Model**

**Author**: Chao Qin | Xingyi Juexiao Information Consulting Center | ORCID 0009-0006-2000-5644
**Version**: v1.1 (2026-08-04) — Zenodo record [10.5281/zenodo.21788556](https://zenodo.org/records/21788556); the versioned DOI above always points to the latest release.

## Summary

- **81% of 274 gold-standard defect/alloy systems are underestimated** by DFT endpoint interpolation (median +0.80 eV, MAE 0.87 eV).
- Deviation grows monotonically with the alloy's experimental band gap: underestimation rate 38% at [0, 0.5) eV → **100% at [2, 4) eV**.
- Endpoint quality matters: physically gapless endpoints (HgTe, BiTe, PbTe, SnTe, HgSe) give +0.45 eV; computationally gapless endpoints (DFT predicts 0 for non-semimetals) give +1.19 eV.
- **Mainstream DFT databases contain zero defective systems** (verified: OQMD, JARVIS-3D, fused V9) — a structural blind spot for infrared, thermoelectric and optoelectronic materials.
- A falsifiable hypothesis framework survived in-sample falsification of half its predictions and yields a correction model:
  `dev = 0.714×gap + 0.232×phys-gapless + 0.401×comput-gapless − 0.402` (MAE 0.28 eV, R² = 0.79).
- Forward predictions for 275 unresolved systems are provided for external testing.

## Reproduce

```bash
python3 verify_defect_gap_stats.py   # prints all reported statistics
```

Requires: Python 3, numpy. Input data (`defect_gap_pipeline_v1.json`) is included.

## Data

| File | Description |
|---|---|
| `defect_gap_pipeline_v1.json` | 274 matched (formula, expt gap, endpoints, prediction, deviation) + 275 unmatched |
| `verify_defect_gap_stats.py` | Reproduction script |

## Data sources

- Experimental gaps: matminer `expt_gap` (literature-verified, 6,354 entries; 549 non-stoichiometric subset)
- DFT endpoint gaps: JARVIS-DFT 3D (OptB88vdW; NIST public domain)
- License: manuscript CC BY 4.0; code MIT
