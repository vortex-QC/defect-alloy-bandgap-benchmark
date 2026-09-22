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

---

## 2026-09-22 Update — Channel-dependent bias: a testable diagnostic

A challenge-shift re-read of this benchmark (asking *why* the bias is where it is, not just how large) yields a refined, falsifiable claim:

> **Systems containing a 4f channel carry an extra band-gap deviation of +0.51 eV (after controlling for the gap), which d-electron systems do not show.**

| Test | Result |
|---|---|
| `dev ~ gap + has_f` (has_f = contains lanthanide/actinide) | coef **+0.511**, p=1.6e-7; ΔR² = **+0.0255** |
| Gap-matched subset (gap ∈ [0.35, 2.61] eV) | with-f (n=16) dev=+1.429 vs without (n=216) +0.867; **+0.562 eV** (p=1.6e-2) |
| Permutation test (10,000 label shuffles) | real ΔR² above the 95th percentile by an order of magnitude; **empirical p = 0.0000** |
| Leave-one-out (drop each f sample) | coefficient stays in [+0.48, +0.56], all positive; worst p = 1.3e-6 |
| Confound check | r(gap, has_f proxy) = +0.17 — the f-term is not a gap proxy |

**d-electron systems show no extra term** (has_d ΔR² = +0.0004, p = 0.50; 3d/4d/5d tested separately, none significant).

**Mechanistic reading** (candidate, supported by the above, not proven): a DFT calculation is a *projection* — Kohn-Sham single-determinant + functional. What such a projection preserves is instantaneous structure (density, band shape); what it discards is multi-configurational structure (near-degeneracies, spin-state competition). 4f electrons are the most localized, hence the most multi-configurational, hence the most affected — the bias is **channel-dependent by construction**, not a matter of functional tuning.

**Falsifiable next steps** (we welcome anyone running these):
1. **Independent-dataset replication** — the f-term should survive on other defect-alloy datasets; if it vanishes, our claim is dead.
2. **Mechanism discrimination** — recompute f-containing systems with a multi-configurational method (DMFT / multi-reference / hybrid+U); the mechanism predicts the f-gap bias should shrink markedly, while non-f systems improve far less.
3. **Direct correlation** — test whether a multi-configurational metric (near-degeneracy density, f-occupation fluctuation) correlates with the bias more strongly than the gap does.

No priority is claimed on this — the goal is verification, not credit. Use, test, refute freely.

**Files** (added 2026-09-22): `prereg_Bprime_v0.1.md`, `prereg_Bprime_v0.2.md` (criteria locked before running), `test_Bprime_v0.1.py`, `test_Bprime_v0.2.py`, `result_Bprime_v0.1.json`, `result_Bprime_v0.2.json`, `RESULT_Bprime_v0.2.md`.

**Honest limits**: 16 f-containing systems (double robustness handled, but small); composition-level proxy (not electronic-structure-level); conclusions limited to this 304-system dataset; the absence of a d-signal is not yet explained.
