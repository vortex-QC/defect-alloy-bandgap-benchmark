#!/usr/bin/env python3
"""Defect/Alloy band-gap benchmark - full v2 pipeline (304 matched + 245 unmatched).

Builds the analysis set from expt_gap (matminer) + JARVIS-3D OptB88 endpoints.
Branches:
  B1: two cations + one anion  -> AC + BC (solid solution)
  B2: one cation + two anions   -> AC + AD (anion solid solution)
  B3: elemental alloy (GeSi)    -> elements
  B4: defect compound (1:1-ish) -> integer parent compound (e.g. Pt0.97S2 -> PtS2)
"""
import json
import gzip
import re
import statistics

METALS = set('Li Na K Rb Cs Fr Be Mg Ca Sr Ba Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Al Ga In Tl Sn Pb Bi La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Th U'.split())
NONMETALS = set('B C N O F Si P S Cl Se Br Te I At H Ge'.split())

def parse_formula(f):
    return [(el, float(c) if c else 1.0) for el, c in re.findall(r'([A-Z][a-z]?)(\d*\.?\d*)', f)]

def load_expt_gap(path="/home/vortex/.matminer/datasets/expt_gap.json.gz"):
    with gzip.open(path) as f:
        d = json.load(f)
    return d["data"]  # [formula, gap]

def load_jarvis_bg(path="/home/vortex/涡肉身壳/field_reasoner/harmonic_db/jarvis_dft_3d.json"):
    with open(path, encoding="utf-8") as f:
        j3d = json.load(f)
    bg = {}
    for v in j3d:
        f2 = v.get("formula", "")
        g = v.get("optb88vdw_bandgap")
        if g is not None and f2 not in bg:
            bg[f2] = g
    return bg

def build(exp_rows, jarvis_bg):
    results, unmatched = [], []
    for f, g in exp_rows:
        if not re.search(r"\d\.\d", f):
            continue  # only non-stoichiometric
        els = parse_formula(f)
        cations = [(e, c) for e, c in els if e in METALS]
        anions = [(e, c) for e, c in els if e in NONMETALS]
        n_cat, n_an = len(cations), len(anions)
        pred, rtype, endpoints = None, "", ""
        if n_cat == 2 and n_an == 1 and len(els) == 3:
            (a, xa), (b, xb) = cations; c = anions[0][0]; tot = xa + xb
            ga, gb = jarvis_bg.get(f"{a}{c}"), jarvis_bg.get(f"{b}{c}")
            if ga is not None and gb is not None:
                pred = ga*(xa/tot) + gb*(xb/tot)
                rtype, endpoints = "solid_solution", f"{a}{c}({ga:.2f})+{b}{c}({gb:.2f})"
        elif n_cat == 1 and n_an == 2 and len(els) == 3:
            a = cations[0][0]; (c1, x1), (c2, x2) = anions; tot = x1 + x2
            ga, gb = jarvis_bg.get(f"{a}{c1}"), jarvis_bg.get(f"{a}{c2}")
            if ga is not None and gb is not None:
                pred = ga*(x1/tot) + gb*(x2/tot)
                rtype, endpoints = "anion_solid", f"{a}{c1}({ga:.2f})+{a}{c2}({gb:.2f})"
        elif n_cat == 0 and n_an >= 2 and len(els) == 2:
            (e1, c1), (e2, c2) = els[0], els[1]
            g1, g2 = jarvis_bg.get(e1), jarvis_bg.get(e2)
            if g1 is not None and g2 is not None:
                tot = c1 + c2
                pred = g1*(c1/tot) + g2*(c2/tot)
                rtype, endpoints = "binary_elem", f"{e1}({g1:.2f})+{e2}({g2:.2f})"
        elif n_cat == 1 and n_an == 1 and len(els) == 2:
            (a, ca), (b, cb) = els
            ratio = cb / ca
            best = min([(1,1),(1,2),(2,1),(1,3),(3,1),(2,3)], key=lambda r: abs(ratio - r[1]/r[0]))
            parent = f"{a}{best[0] if best[0] != 1 else ''}{b}{best[1] if best[1] != 1 else ''}"
            g = jarvis_bg.get(parent)
            if g is not None:
                pred = g
                rtype, endpoints = "defect_parent", f"{parent}({g:.2f})"
        if pred is not None:
            results.append({"formula": f, "expt": g, "type": rtype, "endpoints": endpoints,
                            "pred": round(pred, 3), "dev": round(g - pred, 3)})
        else:
            unmatched.append(f)
    return results, unmatched

def main():
    exp_rows = load_expt_gap()
    jarvis_bg = load_jarvis_bg()
    results, unmatched = build(exp_rows, jarvis_bg)
    n_total = sum(1 for f, g in exp_rows if re.search(r"\d\.\d", f))
    out = {"n_total": n_total, "n_matched": len(results), "results": results, "unmatched": unmatched}
    with open("/home/vortex/涡肉身壳/papers/缺陷合金带隙复现脚本/defect_gap_pipeline_v2.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"Saved: {len(results)} matched / {n_total} total; {len(unmatched)} unmatched")

    devs = [r["dev"] for r in results]
    print(f"Mean {statistics.mean(devs):+.3f} | Median {statistics.median(devs):+.3f} | MAE {statistics.mean(abs(d) for d in devs):.3f} | Under {100*sum(1 for d in devs if d>0)/len(devs):.0f}%")

if __name__ == "__main__":
    main()
