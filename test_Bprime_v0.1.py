# -*- coding: utf-8 -*-
"""难题 B 重读 · 预言 B′ 数据检验（预注册 v0.1 判据锁定执行）
数据：papers/缺陷合金带隙复现脚本/defect_gap_pipeline_v2.json（304 条）
判据：主=dev~expt+MC 的 MC 系数>0 且 p<0.05 且 ΔR²≥0.02；副=偏相关>0 显著；副2=含 f 组控制带隙后正差
纪律：零改动原 pipeline；结果含否定可能，如实落盘
"""
import json
import re
import numpy as np
from scipy import stats

SRC = "papers/缺陷合金带隙复现脚本/defect_gap_pipeline_v2.json"
OUT = "papers/难题B重读_预言B_prime_结果.json"

# ---- 元素权重表（预注册锁定：f=1.0 / 3d=0.5 / 4d,5d=0.3 / sp=0）----
F_ELEMS = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
           "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm"]
D3 = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
D45 = ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
       "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]


def elem_weights(entries):
    """解析化学式 → {elem: 计量数}"""
    out = {}
    for sym, num in entries:
        n = float(num) if num else 1.0
        out[sym] = out.get(sym, 0.0) + n
    return out


def parse(formula):
    return re.findall(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)", formula)


def mc_of(formula):
    """加权局域电子分 + 二值标记"""
    ents = [(s, n) for s, n in parse(formula) if s]
    comp = elem_weights(ents)
    mc = 0.0
    for e, n in comp.items():
        if e in F_ELEMS:
            mc += 1.0 * n
        elif e in D3:
            mc += 0.5 * n
        elif e in D45:
            mc += 0.3 * n
    has_f = 1.0 if any(e in F_ELEMS for e in comp) else 0.0
    has_d = 1.0 if any((e in D3 or e in D45) for e in comp) else 0.0
    return mc, has_f, has_d


def ols(y, X):
    """最小二乘 + 系数标准误/t/p（含截距由调用方加入）"""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n, k = X.shape
    dof = n - k
    s2 = (resid @ resid) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * (1 - stats.t.cdf(np.abs(t), dof))
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - (resid @ resid) / ss_tot
    return beta, se, t, p, r2, resid


def partial_corr(x, y, z):
    """corr(x, y | z)：双方对 z 回归取残差后的相关"""
    Xz = np.column_stack([np.ones_like(z), z])
    rx = x - Xz @ np.linalg.lstsq(Xz, x, rcond=None)[0]
    ry = y - Xz @ np.linalg.lstsq(Xz, y, rcond=None)[0]
    return stats.pearsonr(rx, ry)


def main():
    d = json.load(open(SRC, encoding="utf-8"))
    recs = d["results"]
    formulas = [x["formula"] for x in recs]
    dev = np.array([x["dev"] for x in recs], dtype=float)
    expt = np.array([x["expt"] for x in recs], dtype=float)
    mc = np.zeros(len(recs))
    hf = np.zeros(len(recs))
    hd = np.zeros(len(recs))
    for i, f in enumerate(formulas):
        mc[i], hf[i], hd[i] = mc_of(f)

    res = {"n": len(recs), "mc_stats": {}, "corr": {}, "models": {}, "groups": {}, "judgement": {}}

    # ---- 描述统计 ----
    res["mc_stats"] = {
        "mc_mean": float(mc.mean()), "mc_std": float(mc.std(ddof=1)),
        "mc_max": float(mc.max()), "n_mc_pos": int((mc > 0).sum()),
        "n_has_f": int(hf.sum()), "n_has_d": int(hd.sum()),
    }

    # ---- 相关矩阵（含共线性检查）----
    r_de, p_de = stats.pearsonr(dev, expt)
    r_dm, p_dm = stats.pearsonr(dev, mc)
    r_em, p_em = stats.pearsonr(expt, mc)
    res["corr"] = {
        "dev~expt": {"r": float(r_de), "p": float(p_de)},
        "dev~mc": {"r": float(r_dm), "p": float(p_dm)},
        "expt~mc": {"r": float(r_em), "p": float(p_em)},
        "collinearity_warning": bool(abs(r_em) > 0.7),
    }

    # ---- 模型 1: dev ~ expt ----
    X1 = np.column_stack([np.ones_like(expt), expt])
    b1, se1, t1, p1, r2_1, _ = ols(dev, X1)

    # ---- 模型 2: dev ~ expt + mc ----
    X2 = np.column_stack([np.ones_like(expt), expt, mc])
    b2, se2, t2, p2, r2_2, _ = ols(dev, X2)

    res["models"] = {
        "M1_dev~expt": {"beta": b1.tolist(), "se": se1.tolist(), "t": t1.tolist(),
                        "p": p1.tolist(), "R2": float(r2_1)},
        "M2_dev~expt+mc": {"beta": b2.tolist(), "se": se2.tolist(), "t": t2.tolist(),
                           "p": p2.tolist(), "R2": float(r2_2),
                           "delta_R2": float(r2_2 - r2_1),
                           "mc_coef": float(b2[2]), "mc_p": float(p2[2])},
    }

    # ---- 副判据 1：偏相关 ----
    rp, pp = partial_corr(dev, mc, expt)
    res["partial_corr_dev_mc_given_expt"] = {"r": float(rp), "p": float(pp)}

    # ---- 副判据 2：分组（含 f vs 不含 f），并做带隙匹配的粗略控制 ----
    g_f = dev[hf == 1]
    g_nf = dev[hf == 0]
    e_f = expt[hf == 1]
    e_nf = expt[hf == 0]
    res["groups"] = {
        "has_f": {"n": int(len(g_f)), "dev_mean": float(g_f.mean()), "dev_std": float(g_f.std(ddof=1)),
                  "expt_mean": float(e_f.mean())},
        "no_f": {"n": int(len(g_nf)), "dev_mean": float(g_nf.mean()), "dev_std": float(g_nf.std(ddof=1)),
                 "expt_mean": float(e_nf.mean())},
        "diff_dev": float(g_f.mean() - g_nf.mean()),
        "diff_expt": float(e_f.mean() - e_nf.mean()),
        "ttest_dev": [float(stats.ttest_ind(g_f, g_nf, equal_var=False).statistic),
                      float(stats.ttest_ind(g_f, g_nf, equal_var=False).pvalue)],
        # 带隙匹配子集：expt 落在重叠区间内的样本
        "matched": None,
    }
    lo, hi = max(e_f.min(), e_nf.min()), min(e_f.max(), e_nf.max())
    m_f = dev[(hf == 1) & (expt >= lo) & (expt <= hi)]
    m_nf = dev[(hf == 0) & (expt >= lo) & (expt <= hi)]
    if len(m_f) >= 3 and len(m_nf) >= 3:
        tt = stats.ttest_ind(m_f, m_nf, equal_var=False)
        res["groups"]["matched"] = {
            "band": [float(lo), float(hi)],
            "has_f": {"n": int(len(m_f)), "dev_mean": float(m_f.mean())},
            "no_f": {"n": int(len(m_nf)), "dev_mean": float(m_nf.mean())},
            "diff": float(m_f.mean() - m_nf.mean()),
            "ttest": [float(tt.statistic), float(tt.pvalue)],
        }

    # ---- 判定 ----
    j = res["judgement"]
    j["criterion_main"] = {"mc_coef>0": bool(b2[2] > 0), "mc_p<0.05": bool(p2[2] < 0.05),
                           "delta_R2>=0.02": bool((r2_2 - r2_1) >= 0.02)}
    j["main_pass"] = bool((b2[2] > 0) and (p2[2] < 0.05) and ((r2_2 - r2_1) >= 0.02))
    j["criterion_sub1"] = {"partial_r>0": bool(rp > 0), "p<0.05": bool(pp < 0.05)}
    j["sub1_pass"] = bool((rp > 0) and (pp < 0.05))
    j["Verdict"] = "B' 成立" if j["main_pass"] else ("B' 部分成立（副判据过主判据未过）" if j["sub1_pass"] else "B' 未成立（证伪条件）")

    json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    # ---- 打印 ----
    print(f"n={res['n']}  MC>0: {res['mc_stats']['n_mc_pos']}  含f: {res['mc_stats']['n_has_f']}  含d: {res['mc_stats']['n_has_d']}")
    print(f"\n共线性检查 r(expt, MC) = {r_em:+.3f} (p={p_em:.2e})  {'⚠️共线' if abs(r_em)>0.7 else 'OK'}")
    print(f"\n相关: dev~expt {r_de:+.3f} (p={p_de:.2e}) | dev~MC {r_dm:+.3f} (p={p_dm:.2e})")
    print(f"\nM1 dev~expt      : R²={r2_1:.4f}")
    print(f"M2 dev~expt+MC   : R²={r2_2:.4f}  ΔR²={r2_2-r2_1:+.4f}")
    print(f"   MC 系数={b2[2]:+.3f} (se={se2[2]:.3f}, t={t2[2]:+.2f}, p={p2[2]:.3e})")
    print(f"   expt 系数={b2[1]:+.3f} (p={p2[1]:.3e})")
    print(f"\n偏相关 corr(dev, MC | expt) = {rp:+.3f} (p={pp:.3e})")
    print(f"\n分组: 含f n={len(g_f)} dev={g_f.mean():+.3f} (expt={e_f.mean():.2f}) | 不含f n={len(g_nf)} dev={g_nf.mean():+.3f} (expt={e_nf.mean():.2f})")
    print(f"   diff_dev={g_f.mean()-g_nf.mean():+.3f}  Welch t p={res['groups']['ttest_dev'][1]:.3e}")
    if res["groups"]["matched"]:
        m = res["groups"]["matched"]
        print(f"   带隙匹配子集 [{m['band'][0]:.2f},{m['band'][1]:.2f}]: 含f n={m['has_f']['n']} dev={m['has_f']['dev_mean']:+.3f} | 不含 n={m['no_f']['n']} dev={m['no_f']['dev_mean']:+.3f} | diff={m['diff']:+.3f} p={m['ttest'][1]:.3e}")
    print(f"\n★判定: {j['Verdict']}")
    print(f"   主判据: {j['criterion_main']} -> {j['main_pass']}")
    print(f"   副判据1: {j['criterion_sub1']} -> {j['sub1_pass']}")
    print(f"\n结果落盘: {OUT}")


if __name__ == "__main__":
    main()
