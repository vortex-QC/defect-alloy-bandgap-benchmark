# Systematic Underestimation of Defect and Alloy Band Gaps by DFT Endpoint Interpolation: A 549-Entry Experimental Benchmark with a Hypothesis-Driven Correction Model

**Author**: Chao Qin
**Affiliation**: Xingyi Juexiao Information Consulting Center, Xingyi 562400, Guizhou, China
**ORCID**: 0009-0006-2000-5644
**Version**: v1.2 preprint (2026-08-04)
**Status**: All numbers reproducible from included data pipeline

---

## Abstract

Density functional theory (DFT) systematically underestimates band gaps, and a recent result (arXiv:2026-04) proved that DFT+U eigen-spectrum gaps are *not* valid for defective or isolated systems. For defect and alloy systems (non-stoichiometric compounds, solid solutions, doped semiconductors), no systematic quantification of DFT-side prediction error exists — in part because mainstream DFT databases contain no defective entries at all (we verify 0 non-stoichiometric entries in OQMD, JARVIS-3D, and the fused V9 library). Here we quantify the error of the standard practical approach — endpoint (virtual-crystal) interpolation of DFT endpoint band gaps — against 549 experimental band gaps of non-stoichiometric systems from the matminer expt_gap gold standard. Of 304 systems with resolvable endpoints, **82% are underestimated**, with median deviation +0.84 eV (MAE 0.87 eV). The deviation is strongly structured: it grows monotonically with the alloy's experimental band gap (underestimation rate 38% at [0,0.5) eV vs 100% at [2,4) eV) and depends on endpoint quality (endpoints that are physically gapless give deviations of +0.45 eV; endpoints that DFT incorrectly predicts as gapless give +1.19 eV). Under a hypothesis-driven framework ("gap as formation difficulty of an interstitial object"), we formulated six falsifiable predictions; three were falsified in-sample, and the surviving model — deviation = 0.70×band gap + endpoint-quality terms — explains 78% of the variance (5-fold out-of-sample MAE 0.27 eV, R2 = 0.78). Forward predictions for the 245 unresolved systems are provided for external testing. The work identifies a structural blind spot of computational materials screening: the systems most relevant to infrared, thermoelectric and optoelectronic applications are precisely those with no DFT labels and the largest interpolation errors.

## 1. Introduction

Density functional theory (DFT) in its local and semi-local forms underestimates band gaps (Perdew & Levy 1983; Sham & Schlüter 1983). Corrective approaches (hybrid functionals, GW, DFT+U) improve this at high computational cost, and a recent rigorous result (arXiv 2026-04-28) established that even DFT+U eigen-spectrum gaps, while valid for pristine periodic systems, **are not valid for defective or isolated systems**. Defect and alloy systems — non-stoichiometric compounds, solid solutions, doped semiconductors — are at the heart of infrared detectors (HgCdTe), thermoelectrics (PbSnTe), microelectronics (SiGe) and optoelectronics (II-VI and III-V alloys).

For such systems, the practical community standard is *endpoint interpolation*: predict the alloy band gap from the DFT band gaps of its endpoints (virtual-crystal approximation, VCA). To our knowledge, **no systematic, gold-standard quantification of this error exists**. Two structural reasons:

1. **DFT databases contain no defective systems.** We verify that the OQMD compositions table (integer stoichiometry), the fused V9 library (1,407,278 entries) and JARVIS-3D (93,902 entries) contain **zero** non-stoichiometric formulas. Defective systems require supercells; mainstream databases do not carry them. Consequently, machine-learning screening pipelines (e.g., 200,000-scale property databases) have no defect/alloy labels at all — a blind spot flagged by the recent materials-AI benchmark critique (arXiv 2026-07-22: benchmarks based mainly on computational labels).
2. **Experimental data exist but are unpaired.** The matminer expt_gap gold standard contains 6,354 experimental gaps, of which 549 (8.6%) are non-stoichiometric — i.e., precisely the systems with no DFT labels.

This paper quantifies the endpoint-interpolation error against these 549 experimental gaps, documents a strongly structured deviation, and — under an explicitly falsifiable hypothesis — derives a correction model explaining 79% of the variance, with forward predictions for the unresolved systems.

## 2. Data and Methods

### 2.1 Data

| Dataset | Source | Use |
|---|---|---|
| Experimental band gaps | matminer expt_gap (6,354 entries, literature-verified) | Gold standard |
| Non-stoichiometric subset | 549 entries (8.6%) — ternary solid solutions 319, quaternary 191, binary 39 | Analysis set |
| DFT endpoint gaps | JARVIS-DFT 3D (93,902 entries, OptB88vdW band gaps; 65,395 unique formulas) | Endpoint values |
| DFT database coverage check | OQMD compositions table; fused V9 library (1,407,278); JARVIS-3D | Blind-spot verification |

### 2.2 Non-stoichiometry detection and endpoint decomposition

A formula is classified as non-stoichiometric if it contains a fractional coefficient (e.g., Hg0.7Cd0.3Te). Endpoint decomposition rules (cations = metallic elements; anions = non-metallic):

- **Two cations, one anion** (A_xB_1-x)C → endpoints AC, BC (solid solution)
- **One cation, two anions** A(C_xD_1-x) → endpoints AC, AD (anion solid solution)
- **Elemental alloys** (e.g., GeSi) → endpoints = elements
- **Defect compounds** (e.g., Pt0.97S2) → endpoint = integer parent compound (PtS2); **elemental alloys** (GeSi) → elemental endpoints

Endpoint gaps taken from JARVIS-3D OptB88vdW; interpolation is linear in mole fraction (VCA).

### 2.3 Deviation measure

dev = E_expt − E_pred (positive = DFT endpoint interpolation underestimates the gap).

## 3. Results

### 3.1 Overall

| Metric | Value |
|---|---|
| Matched systems | 304 / 549 (55%) |
| Mean deviation | **+0.745 eV** |
| Median deviation | **+0.837 eV** |
| MAE | 0.865 eV |
| **Underestimation rate** | **82%** |

### 3.2 Structured deviation I: band-gap dependence (monotone)

| Experimental gap range (eV) | n | Mean dev | Underest. rate |
|---|---|---|---|
| [0, 0.5) | 67 | −0.177 | 39% |
| [0.5, 1) | 51 | +0.436 | 78% |
| [1, 2) | 117 | +0.976 | 98% |
| [2, 4) | 69 | +1.476 | **100%** |

### 3.3 Structured deviation II: anion family

| Family | n | Mean dev | Underest. rate |
|---|---|---|---|
| S | 33 | +1.455 | 100% |
| P | 38 | +0.889 | 97% |
| Se | 127 | +0.662 | 72% |
| Te | 72 | +0.465 | 82% |

### 3.4 Structured deviation III: endpoint quality

Endpoints whose DFT gap is ≈0 split into two classes with very different behavior:

| Endpoint class | n | Mean dev |
|---|---|---|
| **Physically gapless** (HgTe, BiTe, PbTe, SnTe, HgSe — intrinsic semimetals) | 70 | **+0.446** |
| **Computationally gapless** (DFT predicts 0 but material is not a semimetal, e.g., EuSe, TaO) | 78 | **+1.186** |

Example: Hg0.7Cd0.3Te (infrared detector) dev +0.19 (HgTe is a physical semimetal); Sr0.5TaO3 dev +2.57 (TaO is a computational-zero endpoint).

## 4. Hypothesis-driven analysis

We adopt an explicitly falsifiable heuristic ("gap as formation difficulty of an interstitial object"): the band gap is the difficulty of forming a delocalized interstitial object between two oscillating field regions; DFT is an axle-based (occupied/unoccupied orbital) description that is accurate when the interstitial object forms easily (narrow gap) and fails when it does not (wide gap). Six falsifiable predictions were derived:

| ID | Prediction | In-sample test | Verdict |
|---|---|---|---|
| P1 | dev grows with **endpoint gap mean** | r = 0.080 | falsified |
| P2 | dev grows with **endpoint gap difference** | r = 0.174 | weak |
| P3 | alloys with **zero-gap endpoints** have small dev | under-est. rate 91% vs 69% | falsified (coarse) |
| P4 | dev grows with **experimental alloy gap** | 38%→100% monotone | **supported** |
| P5 | **physically gapless** endpoints give small dev; **computationally gapless** give large dev | +0.446 vs +1.186 | **supported** |
| P6 | dev grows with **bowing parameter** (interference strength) | r = 0.100 (7 families) | falsified |

Three of six predictions were falsified in-sample — the framework is falsifiable, not post-hoc. The surviving model:

**dev = 0.698 × E_gap(expt) + 0.208 × (phys. gapless endpoint) + 0.456 × (comput. gapless endpoint) − 0.364**

MAE = 0.269 eV, R² = 0.785 (n = 304). 5-fold out-of-sample: MAE = 0.274 eV, R² = 0.776 (band-gap-only baseline: MAE = 0.329 eV).

## 5. Forward predictions (unresolved 245 systems)

| Class | n | Prediction |
|---|---|---|
| Complex multi-cation systems | ~107 | uncertain (to be refined) |
| Binary defect compounds (endpoint lookup) | 65 | endpoint gap directly usable |
| Quaternary (double interference) | 60 | medium-to-large deviation |
| Solid solutions with missing endpoints (semimetal elements) | 39 | small deviation (physically-gapless class) |
| Solid solutions with missing endpoints (other) | 4 | medium deviation |

These are published as testable predictions; external verification is invited.

## 6. Conclusions

1. DFT databases contain **zero** defective systems — the defect/alloy band-gap problem has no computational labels at all.
2. Endpoint interpolation underestimates defect/alloy gaps in **82%** of 304 gold-standard systems; the error grows monotonically with the alloy gap (100% underestimation above 2 eV) and depends on endpoint quality (physically vs computationally gapless endpoints).
3. A falsifiable hypothesis framework survives in-sample falsification of half its predictions and yields a correction model with **out-of-sample R² = 0.78** (MAE 0.27 eV).
4. The systems most relevant to infrared, thermoelectric and optoelectronic applications are precisely those with no DFT labels and the largest interpolation errors — a structural blind spot of computational materials screening.

## Data availability

- Analysis pipeline and results: reproducible from the included scripts (GitHub repository to be linked; DOI to be assigned on Zenodo)
- Experimental data: matminer expt_gap (literature-verified, per-entry sources in the dataset)
- DFT data: JARVIS-DFT 3D (NIST public domain)
- The 549-entry analysis set and the 274 matched results are provided as supplementary JSON.

## References (selected)

- Perdew & Levy 1983, PRL 51, 1884; Sham & Schlüter 1983, PRL 51, 1888
- Validity of DFT+U band gaps in all its known functional forms (arXiv 2026-04-28)
- Generative and multimodal AI for materials prediction and design: Progress, challenges, and perspectives (arXiv 2026-07-22)
- JARVIS-DFT: Choudhary et al., npj Comput. Mater. (2018)
- matminer: Ward et al., Comput. Mater. Sci. (2018)
