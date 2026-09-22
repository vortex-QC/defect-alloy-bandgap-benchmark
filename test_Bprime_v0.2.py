# -*- coding: utf-8 -*-
"""难题 B 重读 · 预言 B′ 复跑 v0.2（预注册 v0.2 判据锁定执行）
主代理=has_f（4f 通道存在性，机理选定）；副代理=MC_f（f 元素计量数）；对照=MC_v0.1
稳健性：置换检验 10000 次 + LOFF 留一法 + v0.1 对照复算
纪律：不改变 v0.1 主判定；结果无论过否如实落盘
"""
import json, re
import numpy as np
from scipy import stats

SRC = "papers/缺陷合金带隙复现脚本/defect_gap_pipeline_v2.json"
OUT = "papers/难题B重读_预言B_prime_结果_v0.2.json"

F_ELEMS = ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Ac","Th","Pa","U","Np","Pu","Am","Cm"]
D3 = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
D45 = ["Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]


def parse(formula):
    return [(s, float(n) if n else 1.0) for s, n in re.findall(r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)", formula) if s]


def proxies(formula):
    ents = parse(formula)
    elems = [e for e, _ in ents]
    has_f = 1.0 if any(e in F_ELEMS for e in elems) else 0.0
    mc_f = sum(n for e, n in ents if e in F_ELEMS)
    mc_v01 = sum(n * (1.0 if e in F_ELEMS else 0.5 if e in D3 else 0.3 if e in D45 else 0.0) for e, n in ents)
    return has_f, mc_f, mc_v01


def ols(y, *cols):
    X = np.column_stack([np.ones_like(y)] + list(cols))
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    n, k = X.shape
    dof = n - k
    s2 = (r @ r) / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = b / se
    p = 2 * (1 - stats.t.cdf(np.abs(t), dof))
    r2 = 1 - (r @ r) / ((y - y.mean()) ** 2).sum()
    return r2, b, se, t, p


def main():
    d = json.load(open(SRC, encoding="utf-8"))
    recs = d["results"]
    n = len(recs)
    dev = np.array([x["dev"] for x in recs]); expt = np.array([x["expt"] for x in recs])
    hf = np.zeros(n); mcf = np.zeros(n); mcv = np.zeros(n)
    for i, x in enumerate(recs):
        hf[i], mcf[i], mcv[i] = proxies(x["formula"])

    res = {"n": n, "n_has_f": int(hf.sum()), "models": {}, "robustness": {}, "judgement": {}}

    r2_1, b1, se1, t1, p1 = ols(dev, expt)
    res["models"]["M1_dev~expt"] = {"R2": float(r2_1), "beta": b1.tolist(), "p": p1.tolist()}

    # 主判据：has_f
    r2_hf, b_hf, se_hf, t_hf, p_hf = ols(dev, expt, hf)
    res["models"]["M2_dev~expt+has_f"] = {
        "R2": float(r2_hf), "delta_R2": float(r2_hf - r2_1),
        "has_f_coef": float(b_hf[2]), "has_f_se": float(se_hf[2]),
        "has_f_t": float(t_hf[2]), "has_f_p": float(p_hf[2]),
        "expt_coef": float(b_hf[1]), "expt_p": float(p_hf[1])}

    # 副判据 A：MC_f
    r2_mf, b_mf, se_mf, t_mf, p_mf = ols(dev, expt, mcf)
    res["models"]["M3_dev~expt+MC_f"] = {
        "R2": float(r2_mf), "delta_R2": float(r2_mf - r2_1),
        "MC_f_coef": float(b_mf[2]), "MC_f_p": float(p_mf[2])}

    # 对照：MC_v0.1（复现一致性）
    r2_mv, b_mv, se_mv, t_mv, p_mv = ols(dev, expt, mcv)
    res["models"]["M4_dev~expt+MC_v0.1(对照)"] = {
        "R2": float(r2_mv), "delta_R2": float(r2_mv - r2_1),
        "MC_coef": float(b_mv[2]), "MC_p": float(p_mv[2])}

    # 副判据 B：带隙匹配子集
    lo, hi = max(expt[hf == 1].min(), expt[hf == 0].min()), min(expt[hf == 1].max(), expt[hf == 0].max())
    m_f = dev[(hf == 1) & (expt >= lo) & (expt <= hi)]
    m_nf = dev[(hf == 0) & (expt >= lo) & (expt <= hi)]
    tt = stats.ttest_ind(m_f, m_nf, equal_var=False)
    res["matched_subset"] = {"band": [float(lo), float(hi)],
                             "n_f": int(len(m_f)), "n_nf": int(len(m_nf)),
                             "dev_f": float(m_f.mean()), "dev_nf": float(m_nf.mean()),
                             "diff": float(m_f.mean() - m_nf.mean()),
                             "t": float(tt.statistic), "p": float(tt.pvalue)}

    # 稳健性 1：置换检验（打乱 has_f 标签 10000 次）
    rng = np.random.default_rng(20260922)
    perm = np.empty(10000)
    for i in range(10000):
        hp = rng.permutation(hf)
        r2p, *_ = ols(dev, expt, hp)
        perm[i] = r2p - r2_1
    emp_p = float((perm >= (r2_hf - r2_1)).mean())
    res["robustness"]["permutation"] = {"n_perm": 10000, "real_deltaR2": float(r2_hf - r2_1),
                                        "perm_mean": float(perm.mean()), "perm_p95": float(np.percentile(perm, 95)),
                                        "empirical_p": emp_p}

    # 稳健性 2：LOO（逐条剔除 f 样本）
    idx_f = np.where(hf == 1)[0]
    loo = []
    for j in idx_f:
        keep = np.ones(n, bool); keep[j] = False
        r2l, bl, sel, tl, pl = ols(dev[keep], expt[keep], hf[keep])
        loo.append({"dropped": recs[j]["formula"], "coef": float(bl[2]), "p": float(pl[2]),
                    "delta_R2": float(r2l - ols(dev[keep], expt[keep])[0])})
    coefs = [x["coef"] for x in loo]
    res["robustness"]["LOO"] = {"n_dropped": len(loo), "coef_min": float(min(coefs)), "coef_max": float(max(coefs)),
                                "sign_all_positive": bool(all(c > 0 for c in coefs)),
                                "p_max": float(max(x["p"] for x in loo)),
                                "p_all_below_0.05": bool(all(x["p"] < 0.05 for x in loo)),
                                "detail": loo}

    # 判定
    j = res["judgement"]
    j["main"] = {"coef>0": bool(b_hf[2] > 0), "p<0.05": bool(p_hf[2] < 0.05),
                 "delta_R2>=0.02": bool((r2_hf - r2_1) >= 0.02),
                 "pass": bool(b_hf[2] > 0 and p_hf[2] < 0.05 and (r2_hf - r2_1) >= 0.02)}
    j["subA"] = {"pass": bool(b_mf[2] > 0 and p_mf[2] < 0.05 and (r2_mf - r2_1) >= 0.02)}
    j["subB"] = {"pass": bool(m_f.mean() - m_nf.mean() > 0 and tt.pvalue < 0.05)}
    j["robust_perm"] = {"pass": bool(emp_p < 0.05)}
    j["robust_loo"] = {"pass": bool(all(c > 0 for c in coefs))}
    j["verdict"] = ("B″ 获判定支持（限本数据集）"
                    if (j["main"]["pass"] and j["subA"]["pass"] and j["robust_perm"]["pass"] and j["robust_loo"]["pass"])
                    else "B″ 未获完整支持（见各项）")

    json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    print(f"n={n} 含f={int(hf.sum())}")
    print(f"\nM1 dev~expt           R²={r2_1:.4f}")
    print(f"M2 dev~expt+has_f     R²={r2_hf:.4f}  ΔR²={r2_hf-r2_1:+.4f} | has_f={b_hf[2]:+.3f} (t={t_hf[2]:+.2f}, p={p_hf[2]:.2e})")
    print(f"M3 dev~expt+MC_f      R²={r2_mf:.4f}  ΔR²={r2_mf-r2_1:+.4f} | MC_f={b_mf[2]:+.3f} (p={p_mf[2]:.2e})")
    print(f"M4 dev~expt+MC_v0.1   R²={r2_mv:.4f}  ΔR²={r2_mv-r2_1:+.4f} (对照，应≈+0.0098)")
    print(f"\n带隙匹配 [{lo:.2f},{hi:.2f}]: 含f n={len(m_f)} dev={m_f.mean():+.3f} | 不含 n={len(m_nf)} dev={m_nf.mean():+.3f} | diff={m_f.mean()-m_nf.mean():+.3f} p={tt.pvalue:.2e}")
    print(f"\n置换检验: 真实ΔR²={r2_hf-r2_1:+.4f} | 置换均值={perm.mean():+.5f} | 95分位={np.percentile(perm,95):+.5f} | 经验p={emp_p:.4f}")
    print(f"LOO: 系数范围 [{min(coefs):+.3f}, {max(coefs):+.3f}] 全正={all(c>0 for c in coefs)} | 最差 p={max(x['p'] for x in loo):.2e} 全<0.05={all(x['p']<0.05 for x in loo)}")
    print(f"\n★判定: {j['verdict']}")
    print(f"   主判据 {j['main']}")
    print(f"   副A(MC_f) {j['subA']['pass']} | 副B(匹配) {j['subB']['pass']} | 置换 {j['robust_perm']['pass']} | LOO {j['robust_loo']['pass']}")
    print(f"\n落盘: {OUT}")


if __name__ == "__main__":
    main()
