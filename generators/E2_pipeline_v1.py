# -*- coding: utf-8 -*-
# ============================================================
# E2_pipeline_v1_5.py  定案线(数据源 float32 → 筛查级效力)
# 修订史:
#   v1     初版
#   v1_2   D-1 钉扎修订(numpy 2.2 / scipy 1.16 / sympy 1.14)
#          A10 float32 识别门 1e-4
#          E2-C1 参考席补球心 O(13 席对 13 顶点)
#   v1_3   E2-C3 签名 inf 污染修复 + isfinite 闸
#          E2-C4 尺度估计改中位数比(顺序无关)
#   v1_4   E2-C5 落盘纪律变更: 报告平铺写入 repo\text\,
#          文件名自带时间戳, 取消日期子文件夹
#   v1_5   E2-C6 dp 全谱对账表索引错位修复:
#          按模型顶点 k 配对席位 ci[k], dp_m[k] vs dp_r[ci[k]]
#          (v1_4 报告中 1.45e+00 量级 |Δ| 为打印伪影, 判废)
# 消费: E1_blend_extract_v1_2_*.json 的 primary(+Z闭合构型)
# 纪律: 只对完整运行成功的输出做记录; 定案效力由闭式/sympy 链另行出具
# ============================================================
import sys, json, math
from datetime import datetime
from pathlib import Path
import numpy as np
import scipy
import sympy as sp
from scipy.optimize import linear_sum_assignment

OUT_ROOT = Path(r"D:\理论体系\repo\text")   # E2-C5: 平铺, 无子文件夹
CODE, VER = "E2_pipeline", "v1_5"

PINS = {"python": (3, 11), "numpy": (2, 2), "scipy": (1, 16), "sympy": (1, 14)}  # D-1
GATE      = 1e-4                       # A10: float32 源识别放行门
SCALE_FLG = 1e-4                       # float32 级尺度偏离阈值
BINS2     = {"l1": 1.0, "l2": 32.0/27.0, "l3": 25.0/19.0}
TOL_LOOSE = 5e-3
TOL_ALIGN = 1e-4
DELTA_EX  = 74*math.sqrt(6)/1377
N_REF     = np.array([7/9, 0.0, -4*math.sqrt(2)/9])

S2, S3, S6 = math.sqrt(2), math.sqrt(3), math.sqrt(6)
REF = {                                # 13 席(E2-C1: 含球心 O)
 "O":  (0.0, 0.0, 0.0),
 "x2": (0.0, 1.0, 0.0),
 "x3": (0.0, 0.5, S3/2),
 "x4": (-S6/3, 0.5, S3/6),
 "x5": (-2*S6/9, 0.5, -7*S3/18),
 "x6": (5*S6/27, 0.5, -23*S3/54),
 "x7": (28*S6/81, 0.5, 17*S3/162),
 "x8": (14*S6/57, -11/38, 49*S3/114),
 "x9": (-2*S6/9, -1/3, 4*S3/9),
 "x10":(-10*S6/27, -1/3, -4*S3/27),
 "x11":(-2*S6/81, -1/3, -44*S3/81),
 "x12":(86*S6/243, -1/3, -52*S3/243),
 "P":  (8*S6/459, -295/306, 139*S3/918),
}
KEYS   = list(REF.keys())
KEYS_A = np.array(KEYS)
X_REF  = np.array([REF[k] for k in KEYS])

STEM = None
LOG = []
def log(s=""): print(s); LOG.append(str(s))

def make_stem():
    global STEM
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    today = datetime.now().strftime("%Y%m%d")
    seq = 1 + len(list(OUT_ROOT.glob(CODE + "_" + VER + "_*_" + today + "-*")))
    STEM = (CODE + "_" + VER + "_%03d_" % seq
            + datetime.now().strftime("%Y%m%d-%H%M%S"))

def write_report(status, extra=None):
    payload = dict(code=CODE, ver=VER, stem=STEM, status=status,
                   created=datetime.now().isoformat(timespec="seconds"),
                   pins={k: ".".join(map(str, v)) for k, v in PINS.items()},
                   gate_float32=GATE, delta_exact="74*sqrt(6)/1377",
                   n_ref=[7/9, 0.0, -4*math.sqrt(2)/9],
                   log="\n".join(LOG), **(extra or {}))
    pj = OUT_ROOT / (STEM + ".json")
    pj.write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=float),
                  encoding="utf-8")
    (OUT_ROOT / (STEM + ".txt")).write_text("\n".join(LOG) + "\n", encoding="utf-8")
    log("[E2] 报告落盘: " + str(pj))

def failfast(msg):
    log("[FAIL-FAST] " + msg); write_report("FAIL_FAST", {"msg": msg}); sys.exit(1)

def pinning_check():
    import platform
    pv = tuple(int(x) for x in platform.python_version().split(".")[:2])
    libs = {"python": pv,
            "numpy": tuple(int(x) for x in np.__version__.split(".")[:2]),
            "scipy": tuple(int(x) for x in scipy.__version__.split(".")[:2]),
            "sympy": tuple(int(x) for x in sp.__version__.split(".")[:2])}
    bad = [k + "=" + ".".join(map(str, v)) + "≠钉扎" + ".".join(map(str, PINS[k]))
           for k, v in libs.items() if v != PINS[k]]
    if bad: failfast("D.1 版本钉扎不符 → 本次运行不产生定案效力: " + "; ".join(bad))
    log("[D.1] 钉扎核验通过(D-1 修订表): " +
        ", ".join(k + "=" + ".".join(map(str, v)) for k, v in libs.items()))

def signatures(X):
    # E2-C3: 升序排序后丢弃末位 inf, 只保留 12 个真实距离
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
    np.fill_diagonal(D, np.inf)
    return np.sort(D, axis=1)[:, :12]

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc, float(np.linalg.det(R))

def main():
    make_stem()
    pinning_check()

    # ---- 输入 ----
    cands = sorted(OUT_ROOT.rglob("E1_blend_extract_v1_2_*.json"))
    if not cands: failfast("未找到 E1 v1_2 输出")
    e1 = cands[-1]; log("[E2] 输入: " + str(e1))
    data = json.loads(e1.read_text(encoding="utf-8"))
    rec = {r["object"]: r for r in data["records"]}
    prim = rec.get("+Z闭合构型")
    if prim is None or prim.get("role") != "primary": failfast("缺 primary 对象")
    Xm = np.array(prim["vertices"], dtype=np.float64)
    if Xm.shape != (13, 3): failfast("primary 形状 " + str(Xm.shape) + " ≠ 13×3")

    # ---- 尺度侦察(E2-C4: 中位数比, 与顶点顺序无关) ----
    dm = np.linalg.norm(Xm - Xm.mean(0), axis=1)
    dr = np.linalg.norm(X_REF - X_REF.mean(0), axis=1)
    scale = float(np.median(dm) / np.median(dr))
    log("[尺度] scale = %.12f  偏离1 = %.3e" % (scale, abs(scale - 1)))
    Xn = Xm / scale
    if abs(scale - 1) > SCALE_FLG:
        log("[尺度][FLAG] 尺度偏离超 float32 阈值, 按登记值归一, 不作定案吸收")

    # ---- 签名指派 ----
    cost = np.linalg.norm(signatures(Xn)[:, None, :]
                          - signatures(X_REF)[None, :, :], axis=2)
    if not np.all(np.isfinite(cost)): failfast("签名代价矩阵含非有限值")
    ri, ci = linear_sum_assignment(cost)
    log("[对应] 匈牙利签名代价和 = %.3e" % cost[ri, ci].sum())
    log("[对应] 席位映射(参考席→模型顶点): " +
        ", ".join(str(KEYS_A[ci[i]]) + "→v" + str(int(ri[i]) + 1) for i in range(13)))
    Q = X_REF[ci]

    # ---- Procrustes ----
    R, T, detR = kabsch(Xn, Q)
    Xa = Xn @ R.T + T
    res = np.linalg.norm(Xa - Q, axis=1)
    tag = "保向" if detR > 0 else "镜像!FLAG"
    log("[对齐] det(R) = %+.6f(%s)  残差 max = %.3e  mean = %.3e"
        % (detR, tag, res.max(), res.mean()))
    if res.max() > TOL_ALIGN: log("[对齐][FLAG] 残差超 1e-4, 结果降级筛查级")

    # ---- 支判认 + 全谱 dp 对账(E2-C6 修复: 按模型顶点 k 配对席位 ci[k]) ----
    i_P = int(np.flatnonzero(KEYS_A[ci] == "P")[0])
    n = N_REF / np.linalg.norm(N_REF)
    dp_m, dp_r = Xa @ n, X_REF @ n
    chi = float(dp_m[i_P])
    log("[支判认] ── dp 全谱对账(模型实测 vs 参考闭式) ──")
    worst = 0.0
    for k in range(13):
        i = int(ci[k])
        d = abs(dp_m[k] - dp_r[i]); worst = max(worst, d)
        log("  %4s: dp_model = %+.7f   dp_ref = %+.7f   |Δ| = %.2e"
            % (KEYS[i], dp_m[k], dp_r[i], d))
    br_tag = "X₋(−Z)" if chi < 0 else "X₊(+Z)"
    log("[支判认] 全谱最大偏差 = %.3e;  χ = n·x_P = %+.9f  → 支 = %s"
        % (worst, chi, br_tag))

    # ---- R19 数值识别(A10 float32 门) ----
    rel = abs(abs(chi) - DELTA_EX) / DELTA_EX
    r19_tag = ("识别命中(float32 裕度达标, 筛查级)" if rel <= GATE
               else "识别失败→退回复核(A.2)")
    log("[R19数值] |χ| = %.9f  vs 74√6/1377 = %.9f  相对残差 = %.3e"
        % (abs(chi), DELTA_EX, rel))
    log("[R19数值] " + r19_tag + "  [注] 1e-12 定案门仅对闭式/sympy 链有效")

    # ---- R17 反射检测 ----
    Xr = Xa - 2.0 * np.outer(dp_m, n)
    reg = [i for i in range(13) if i != i_P]
    Dall = np.linalg.norm(Xr[:, None, :] - Xa[None, :, :], axis=2)
    d_reg = float(np.max(np.min(Dall[reg], axis=1)))
    d_P   = float(np.min(Dall[i_P]))
    dpp   = float(Xr[i_P] @ n)
    ok17 = d_reg < TOL_ALIGN and abs(dpp - DELTA_EX) < 1e-4 and d_P > 0.1
    r17_tag = "通过——两支结构复现" if ok17 else "未通过→登记"
    log("[R17] 正则席反射自映射最大偏差 = %.3e;  Rπ(P) 镜距 = %+.9f (预期 +δ, 偏差 %.2e);  Rπ(P) 距最近席 = %.6f(新位置,非席)"
        % (d_reg, dpp, abs(dpp - DELTA_EX), d_P))
    log("[R17] " + r17_tag)

    # ---- 棱长谱分档(78 全对距) ----
    def bins_of(X):
        D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        d2 = D[np.triu_indices(13, 1)] ** 2
        return {k: int((np.abs(d2 - t) < TOL_LOOSE * t).sum())
                for k, t in BINS2.items()}
    bm, br = bins_of(Xa), bins_of(X_REF)
    cov = sum(bm.values())
    spec_tag = "一致" if bm == br else "不一致[FLAG]"
    log("[棱谱] 模型分档: %s  覆盖 %d/78  谱外 %d(含塌缩合并消失的网格棱)"
        % (str(bm), cov, 78 - cov))
    log("[棱谱] 参考闭式同口径: %s  → %s" % (str(br), spec_tag))

    # ---- aux 四面体旁证(不入定案) ----
    aux_ok = None
    a, b = rec.get("+Z闭合构型δ空间"), rec.get("+Z闭合构型δ镜像空间")
    if a and b:
        A1 = np.array(a["vertices"]); A2 = np.array(b["vertices"])
        f1 = np.sort(np.linalg.norm(A1[:, None, :] - A1[None, :, :], axis=2).ravel())
        f2 = np.sort(np.linalg.norm(A2[:, None, :] - A2[None, :, :], axis=2).ravel())
        aux_ok = bool(np.allclose(f1, f2, rtol=1e-5, atol=1e-6))
        aux_tag = "一致 → 互为镜像(筛查级旁证)" if aux_ok else "不一致[FLAG]"
        log("[aux] 两四面体边长谱" + aux_tag + "  [注] 可视化辅助, 不入定案链")

    # ---- 终态 ----
    log("[溯源] %s | Blender %s | %s"
        % (e1.name, data.get("blender", "?"), data.get("precision_note", "")))
    ok = rel <= GATE and ok17 and res.max() <= TOL_ALIGN and bm == br
    write_report("DONE" if ok else "DONE_FLAGGED",
                 extra=dict(scale=float(scale), chi=chi,
                            branch="X₋" if chi < 0 else "X₊",
                            delta_rel_residual=float(rel), R17_pass=bool(ok17),
                            bins_model=bm, bins_ref=br, aux_mirror=bool(aux_ok),
                            procrustes_max_res=float(res.max()),
                            mapping={KEYS[ci[k]]: int(ri[k]) + 1
                                     for k in range(13)}))

main()
