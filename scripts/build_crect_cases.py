#!/usr/bin/env python
"""build_crect_cases.py
C 心矩形磁亚晶格 (centered rectangular net) 磁体: 超胞 + 共线磁序枚举 + J/J2/J3 能量映射系数。

为什么不能复用 build_rectangular_cases.py —— 两点都是**物理**差异, 不是风格差异:

  (1) 磁亚晶格是 **2 个磁性原子/胞** 的 C 心网。以 CrOCl 单层为例: Cr 在 (0,0,z) 与
      (1/2,1/2,z'), 于是同胞两个 Cr 必须能取**反平行** —— 那正是 C 心网的 Néel 态。
      矩形脚本的磁序只由 (n0,n1) 决定, 两个 basis 原子被强行赋同一个 σ,
      「亚晶格反平行」这个态根本表达不出来。
  (2) 最近两个 Cr-Cr 壳 (CrOCl: 3.035 A 与 3.198 A) 只差 0.163 A, 而矩形脚本的
      `magnetic_shells(tol=0.2)` 会把它们**静默合并**成一个 J1。两个物理上完全
      不同的交换路径被算成一个, 而 J 反解照样收敛、不报错 —— 这是最危险的一类错。
      本脚本因此把壳容差默认设为 0.05, 并显式打印相邻壳间距; 间距 < 0.3 A 时告警。

距离壳 (CrOCl 单层, a=3.1980 b=3.8603 A, 2 Cr/胞):
  J1 = 3.0347 A  x4   --- (1/2,1/2) 键, 即 C 心矢量, 连接**同胞的两个 Cr**
  J2 = 3.1980 A  x2   --- 沿 a, 同一亚晶格
  J3 = 3.8603 A  x2   --- 沿 b, 同一亚晶格
⚠ 顺序**不是**按「近邻/次近邻」的直觉: 最短的键反而是把两个 basis 连起来的 C 心矢量。

模型: H = -J1 S·S - J2 S·S - J3 S·S  (J>0 = 铁磁)
能量映射 (每超胞): E = e0 + c1*J1~ + c2*J2~ + c3*J3~,  J~ = J*S^2
  c_k = -0.5 * sum_{第 k 壳有向键} sigma_i * sigma_j

磁序 (basis = 同胞两个 Cr 的亚晶格符号):
  fm        : sigma = +1                         (全同向)
  afm_c     : sigma = basis                      (C 心网 Néel: 同胞两 Cr 反平行)
  afm_a     : sigma = (-1)^n0                    (沿 b 的条纹)
  afm_b     : sigma = (-1)^n1                    (沿 a 的条纹)
  afm_diag  : sigma = (-1)^(n0+n1)               (棋盘)
  afm_c_a   : sigma = basis * (-1)^n0            (对照: 与 afm_a 简并)
  afm_c_diag: sigma = basis * (-1)^(n0+n1)       (对照: 与 afm_diag 简并)

⚠ 核心结构事实 (也是必须钉住的回归不变式): C 心矢量 (1/2,1/2) 连接的是**同一胞内**
  的两个 basis, 所以 basis 符号对一切只依赖 (n0,n1) 奇偶的磁序**完全不可见**:
      afm_c_a == afm_a,  afm_c_b == afm_b,  afm_c_diag == afm_diag
  于是像 (fm, afm_c, afm_a, afm_c_a) 这样的组合**奇异** (det=0) —— 不能选。
  反过来 fm 与 afm_c 只在 c1 上不同 (c2/c3 完全相同), 所以 c1 的分辨完全靠这一对。
  selected = fm / afm_c / afm_a / afm_b, 4x4x1 下 cond(A)=64.02, det=524288 (非奇异)。
  ⚠ 注意这里**不需要**周期 4 的序: 与方形晶格不同, 周期 2 的四个序已经给出
    互不相同的 (c1,c2,c3), J3 与 e0 不简并。

键的重数: J1 每胞 4 根、J2/J3 各 2 根 (fm 行的 per-mag 系数 = (-2,-1,-1))。
  J1 的分数位移恰好是 (1/2,1/2) = 两方向的半周期, 所以与矩形/方形/蜂窝同理
  用显式平移搜索 (nshift=1) 保留像的重数; 实测本几何下 min_image 与
  shift_search 给出相同系数 (因为 ±(1/2,1/2) 两个像对应**不同**的邻位, 不像
  方形 J3=2a 那样折叠到同一位点), 但显式搜索是更稳的选择。

弱铁磁倾斜族 cant (P2′, 系数层; 只有显式给 --cant-theta 才生成, 不给时默认路径逐位不变):
  基序 (默认 afm_c, 可用 --cant-base 改) 的每个原子绕 z 轴倾斜 ±θ、面内方位 φ:
      S_i = (sinθ·cosφ, sinθ·sinφ, σ_i·cosθ),   σ_i = 基序的共线符号
  于是每个键的点积只依赖 θ 的**二倍角**:
      S_i·S_j = sin²θ + σ_iσ_j·cos²θ
        · 同亚晶格 (σ_iσ_j=+1): sin²θ+cos²θ = 1    => 该壳系数**与 θ 无关** (恒等于共线值)
        · 跨亚晶格 (σ_iσ_j=-1): sin²θ-cos²θ = -cos2θ => 该壳系数按 **-cos2θ** 缩放
  对本晶格 (CrOCl 4x4x1, 每胞 2 个 Cr): J1 = C 心键 (1/2,1/2) **跨**亚晶格 (连同胞两个 basis),
  J2/J3 = 同一亚晶格内的 a/b 键, 系数与 θ 无关。每胞系数因此满足**仿射恒等式**
      c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(base)
  等价写法 (注意符号陷阱): c_k(θ) = -sin²θ·(z_k·n_mag_cell/2) + cos²θ·c_k(base), 其中
  z_k = 每格点键数 = (4,2,2) **是正数**, 而 c_k(fm) = -z_k·n_mag_cell/2 是**负数** ——
  两套约定差一个负号, 故报告里两种量都落盘 (per_site_bonds / coeff_fm_per_primitive_cell)。
  为什么**不能**拿 fm 当基序: fm 只有一个自旋亚晶格, ±θ 倾斜对全同向自旋就是一次整体旋转,
  任意键的点积恒为 +1 => 系数行对一切 θ 都等于 fm 行, 该族不携带任何信息 (θ 标度自检还会
  "假通过")。故 --cant-base 只接受 afm_c/afm_a/afm_b, 给 fm 直接报错。

  ⚠ 三条**诚实边界** (P2′ 最重要的结论, 同时写进 cant_report 每一项的 note):
    (1) θ 的**符号不可识别**: cos2θ 是偶函数 => cant_p8 与 cant_m8 的系数**逐位相同**;
        总能量只能定出 cos2θ (进而 |θ|), 定不出倾斜方向。
    (2) 方位 φ 由几何**不定**: φ 不进任何点积 => 系数与 φ 完全无关 (自检专门查这一条)。
        φ 只能靠 SOC/磁各向异性 (真 DFT) 定; 本脚本只透传并落盘。
    (3) 倾斜族**不增加设计矩阵的秩**: 每条 cant 行都是 [1, c(fm)] 与 [1, c(base)] 的仿射组合
        => 它让你换一组构型拟合**同一套** {e0, J1, J2, J3}, 但**不能**分离 DMI 项 D, 也
        **不能**把 J1 与 J1′ (同距离不同亚晶格对) 分开。要分离它们必须有真正独立的自旋
        构型或 SOC 计算 (自旋螺旋/能带) —— 这一条是"needs a real SOC DFT run"的边界。

  cant 序命名规则 (稳定、可解析; 下游按名字找目录):
      "cant_" + ("p" if θ>0 else "m") + "%.12g" % |θ|
      θ=+8 -> cant_p8;  θ=-8 -> cant_m8;  θ=±8.5 -> cant_p8.5 / cant_m8.5;
      θ=0 -> cant_p0 / cant_m0 (退化: 系数与基序逐位相同, 自检 (b) 就用它)
      即: 整数值不带小数点, 其余按 12 位有效数字格式化; 若两个角度格式化成同一个名字,
      或与已有共线序重名, **硬报错** (否则后者的目录会静默覆盖前者)。
  每个输入值 t 生成**两个**序: cant_p|t| (θ=+|t|) 与 cant_m|t| (θ=-|t|); 输入 -8 与 8
  等价 (cos2θ 偶), 只打印提示、不静默改值。
  cant 序一律需要 SOC: mag_vectors.needs_soc("cant*") 为 True => 下游必须 --lsorbit;
  名字也照实落进 magnetic_orders.json 的 soc_orders (共线序一个都不进)。

  JSON 契约 (新增键, 全部向后兼容; 不给 --cant-theta 时全是中性值):
      cant_orders: [序名]             cant_theta_deg: {序名: 带符号 θ}
      cant_azimuth_deg: float         noncollinear: bool (有 cant 序才 True; model_spec 同)
      soc_orders: [需要 SOC 的序名]    cant_report: {序名: {...}}
      cant_report[<序名>] 的键 (solve 端按这些拼): order, base, theta_deg, azimuth_deg,
        coefficients_per_cell (= **每超胞**, 与 report 顶层同名键同口径),
        coeff_per_primitive_cell, coeff_fm_per_primitive_cell, coeff_base_per_primitive_cell,
        prediction_per_primitive_cell, prediction_rule, max_abs_dev, n_cells, n_mag_cell,
        per_site_bonds, intra_bonds_per_cell, inter_bonds_per_cell, theta_dependent_shells,
        theta_independent_shells, soc_required, sign_of_theta_unidentifiable,
        azimuth_not_determined_by_geometry, note
      ⚠ 口径陷阱: 本仓库历史命名里的 `coefficients_per_cell` 实际是**每超胞** (顶层既有键就是
        如此, 自检 (6) 的解析理想值也是每胞 4/2/2 根键的 per-mag 值); cant_report 里同名键
        沿用同口径, **每原胞**的量一律显式叫 `*_per_primitive_cell` (已除以 n_cells),
        闭式预测用的是每原胞。

用法:
  python build_crect_cases.py POSCAR --supercell 4 4 1 --outdir cases \
      --element Cr --nonmag Cl O --moment 3.0 --spin 1.5 [--selftest]
  # 弱铁磁倾斜族 (P2′; 只有显式给 --cant-theta 才生成, MAGMOM 一律三分量):
  python build_crect_cases.py POSCAR --supercell 4 4 1 --outdir cases \
      --element Cr --nonmag Cl O --cant-theta "8 16" [--cant-azimuth 0.0] [--cant-base afm_c]
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from itertools import product

import numpy as np

# 本文件所在 scripts/ 的上一级 = 仓库根; --selftest 省略 POSCAR 时用它去取 tests/fixtures/。
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORDERS = ["fm", "afm_c", "afm_a", "afm_b", "afm_diag", "afm_c_a", "afm_c_diag"]
# 物理最直白且线性无关的 4 个: cond(A)=64.02 @ 4x4x1 (远低于 1e4 闸门), det=524288
SELECTED = ["fm", "afm_c", "afm_a", "afm_b"]

RULES = {
    "fm": "sigma_i = +1 (全同向)",
    "afm_c": "sigma = basis (C 心网 Neel: 同胞两个磁原子反平行)",
    "afm_a": "sigma = (-1)^n0 (沿 b 的条纹)",
    "afm_b": "sigma = (-1)^n1 (沿 a 的条纹)",
    "afm_diag": "sigma = (-1)^(n0+n1) (棋盘)",
    "afm_c_a": "sigma = basis * (-1)^n0 (对照: 与 afm_a 完全简并)",
    "afm_c_diag": "sigma = basis * (-1)^(n0+n1) (对照: 与 afm_diag 完全简并)",
}

# 周期 2 的 basis 无关序 → 它们的 basis 翻倍版本必然与其逐位相同 (见模块 docstring)
BASIS_DEGENERATE = (("afm_c_a", "afm_a"), ("afm_c_diag", "afm_diag"))

SHELL_NAMES = ("J1", "J2", "J3")

# ---------------------------------------------------------------------------
# 弱铁磁倾斜 (cant, P2′) 规格表 —— 声明式, 自检直接引用这张表而不是重抄一遍
# ---------------------------------------------------------------------------
CANT_BASE_DEFAULT = "afm_c"

# 允许当 cant 基序的共线序 + 每个基序**预期**的亚晶格结构。
# 为什么 fm 不在表里: fm 只有一个自旋亚晶格, ±θ 倾斜对全同向的自旋就是一次整体旋转,
#   任意键的点积恒为 +1 => 系数行对一切 θ 都等于 fm 行 (该族零信息, θ 标度自检还会假通过)。
# 为什么默认 afm_c: 只有 C 心网的 Néel 态让 J1 (C 心键 (1/2,1/2)) 恰好**跨**亚晶格、
#   J2/J3 同亚晶格 —— 这是唯一能把 J1 的 θ 标度与 J2/J3 分开的基序。
# pair: σ 图样必须逐位等于哪个声明图样 (守卫用); expect_inter/intra: 预期哪几个壳是
#   跨/同亚晶格 (只对 afm_c 声明: 条纹基序的 C 心键集合里同/跨混在一起, 没有"壳纯"分法,
#   故实测值照实落盘, 不假装有预期)。
CANT_BASES = {
    "afm_c": dict(pair="basis", expect_inter=("J1",), expect_intra=("J2", "J3")),
    "afm_a": dict(pair="stripe_a", expect_inter=None, expect_intra=None),
    "afm_b": dict(pair="stripe_b", expect_inter=None, expect_intra=None),
}

# cant 打开时**预期会变**的既有键 (布尔开关), 以及全部新增键 —— 自检 (a) 的差集比对靠它们。
CANT_DIFF_KEYS = ("noncollinear",)
CANT_NEW_KEYS = ("cant_orders", "cant_theta_deg", "cant_azimuth_deg", "cant_report", "soc_orders")
CANT_NEW_MODEL_KEYS = ("soc_orders", "cant_base")
# 裸 --selftest (没给 --cant-theta) 也要把 cant 代码路径覆盖到: 用这个样本角在子进程里真跑
# 一遍 cant 自检。选 8/16 是因为 cos2θ = 0.9613/0.8480 两个值都远离 0 和 ±1 (θ=0/45/90 是
# 标度关系的退化点, 用它们测不出"到底是 θ 相关还是无关")。
CANT_SELFTEST_SAMPLE = "8 16"

# order_spec / cant_report 里 cant 条目共用的规则串 (一个字都不许含糊, 求解器直接读):
CANT_RULE = ("cant: S_i = (σ_i·sinθ, 0, σ_i·cosθ) (z 反平行的两个自旋亚晶格各自倾斜 ±θ, "
             "面内方位 φ); 同亚晶格键 S_i·S_j = sin²θ+cos²θ = 1 (θ 无关), "
             "跨亚晶格键 = sin²θ-cos²θ = -cos2θ; θ 的**符号**不可辨识 (cos2θ 偶函数), "
             "方位 φ 不由几何决定; 每行 c(θ) = sin²θ·z + cos²θ·c(base) 是 fm 行与基序行的"
             "仿射组合 => 倾斜族不给设计矩阵增加任何新秩")


def sign_pow(base):
    """(-1)^base, 对整数数组。"""
    return np.where(np.asarray(base) % 2 == 0, 1.0, -1.0)


def sigmas_of(n0, n1, basis):
    """七个磁序的 sigma 数组 (±1)。basis 是同胞亚晶格符号 (±1)。"""
    n0 = np.asarray(n0, dtype=int)
    n1 = np.asarray(n1, dtype=int)
    b = np.asarray(basis, dtype=float)
    return {
        "fm": np.ones(len(n0)),
        "afm_c": b.copy(),
        "afm_a": sign_pow(n0),
        "afm_b": sign_pow(n1),
        "afm_diag": sign_pow(n0 + n1),
        "afm_c_a": b * sign_pow(n0),
        "afm_c_diag": b * sign_pow(n0 + n1),
    }


# ---------------------------------------------------------------------------
# cant 参数解析 / 命名 / 亚晶格守卫 (纯函数, 不做 I/O => 守卫在任何计算之前 fail fast)
# ---------------------------------------------------------------------------
def format_cant_theta(theta_deg):
    """把 |θ| 格式化成序名里的角度串 (规则必须稳定, 见模块 docstring「命名规则」)。

    整数值不带小数点 (8 -> "8", 16 -> "16"), 其余按 %.12g (8.5 -> "8.5", 8.25 -> "8.25")。
    12 位有效数字是"够用又能拆穿重名"的折中: 两个角度若差得比这更小, 几乎肯定是输入写错,
    由 cant_order_plan 硬报错, 绝不静默覆盖目录。
    """
    a = abs(float(theta_deg))
    s = "%.12g" % a
    return s[:-2] if s.endswith(".0") else s


def cant_order_name(theta_deg):
    """带符号 θ -> 序名: cant_p8 / cant_m8 / cant_p8.5 (θ 的正负用 p/m 前缀表示)。

    用 copysign 而不是 `t >= 0`: 这样才能把 -0.0 与 +0.0 分成 cant_m0 / cant_p0
    (二者系数与基序逐位相同, 但名字必须不同, 否则 θ=0 的请求会撞名报错)。
    """
    t = float(theta_deg)
    sign = "p" if math.copysign(1.0, t) > 0 else "m"
    return "cant_" + sign + format_cant_theta(t)


def parse_cant_theta(raw):
    """`"8 16"` / `"8,16"` / `"8, 16.5"` -> [8.0, 16.0, 16.5]; 非法值直接 SystemExit。

    刻意**不用** argparse 的 type=float + nargs: 那样空串/逗号形式会得到英文 usage 报错,
    而本仓库要求可执行的中文提示 (哪个值、为什么不行、怎么改)。
    """
    if raw is None or not str(raw).strip():
        raise SystemExit(
            "--cant-theta 是空值; 例: --cant-theta \"8 16\" (空格或逗号分隔的角度, 单位: 度)。\n"
            "  不想生成倾斜族就不要给这个开关 (默认路径与旧版逐位相同)。")
    out = []
    for tok in str(raw).replace(",", " ").split():
        try:
            t = float(tok)
        except ValueError:
            raise SystemExit(
                f"--cant-theta 里的 {tok!r} 不是数; 只接受空格/逗号分隔的角度, "
                "例 --cant-theta \"8 16\"。")
        if not math.isfinite(t):
            raise SystemExit(
                f"--cant-theta 里的 {t} 不是有限值; cant 角度必须是有限实数, 且 |θ| < 90 度。")
        if abs(t) >= 90.0:
            raise SystemExit(
                f"--cant-theta {t} 不合法: 要求 |θ| < 90 度。θ=90° 时 cosθ=0, 两个亚晶格的 z "
                "分量同时归零, 自旋只剩面内分量 —— 「弱铁磁倾斜」退化成面内序, σ 不再出现在任何"
                "点积里 (本族的全部物理都没了); 需要那个态请显式构造面内序。")
        out.append(t)
    return out


def cant_order_plan(thetas, base):
    """[(序名, 带符号 θ)]: 每个输入值生成 +|t| 与 -|t| 两个序, 并做全部**静态**守卫。

    守卫 (硬报错, 绝不静默):
      * 基序不在 CANT_BASES 里 —— 含 fm: 单亚晶格 => ±θ 只是整体自旋旋转, 该族零信息;
      * 序名与已有共线序冲突, 或两个角度格式化成同一个名字 (后者目录会被静默覆盖);
      * 输入角度的符号: 本族按 ±|θ| 成对生成, 输入 -8 与 8 等价 (cos2θ 是偶函数),
        只**打印提示**, 不静默取绝对值。
    """
    if base not in CANT_BASES:
        raise SystemExit(
            f"--cant-base {base!r} 不是允许的 cant 基序; 只能取 {sorted(CANT_BASES)}。\n"
            "  说明: fm (全同向) 只有一个自旋亚晶格, ±θ 倾斜对它只是一次**整体**自旋旋转 —— "
            "任意键的点积恒为 +1, 系数行对一切 θ 都等于 fm 行, 该族不携带信息 (θ 标度自检还会"
            "假通过)。")
    plan, seen = [], {}
    for t in thetas:
        if t < 0:
            print(f"[note] --cant-theta {t}: 本族按 ±|θ| 成对生成, 输入符号没有物理意义 "
                  f"(cos2θ 是偶函数) —— 按 |θ|={abs(t):g} 处理, 不静默改值。")
        for sgn in (1.0, -1.0):
            th = sgn * abs(float(t))
            nm = cant_order_name(th)
            if nm in ORDERS:
                raise SystemExit(
                    f"cant 序名 {nm!r} 与已有共线磁序重名; 该目录的 magmom.txt 会被静默覆盖。"
                    "请改 --cant-theta。")
            if nm in seen:
                raise SystemExit(
                    f"cant 序名冲突: θ={seen[nm]:g} 与 θ={th:g} 都格式化为 {nm!r} "
                    "(命名规则 %.12g); 请把两个角度差开 (>1e-9), 否则后者会静默覆盖前者的目录。")
            seen[nm] = th
            plan.append((nm, th))
    return plan


def cant_sublattice_guard(sig_base, base, basis, n0, n1):
    """核对 cant 的两个自旋亚晶格**就是预期的那一对**; 不符直接 SystemExit。

    这是 P2′ 最容易静默出错的点: 只要 σ 不是"每个原胞恰好一上一下"的图样, spin_cant 照样
    能算出一组漂亮的三分量自旋, 系数照样落盘、反解照样收敛 —— 但倾斜的对象换成了别的配对,
    J1 的 θ 标度就没了。所以用两个**独立判据**:
      (1) σ 必须恰好 2 个取值 (单取值 = 全同向 = 整体旋转, 不是两个亚晶格);
      (2) σ 必须逐位等于基序声明的图样 (afm_c -> 原胞两个 basis; afm_a/b -> 沿 a/b 的条纹)。
    """
    sig = np.asarray(sig_base, dtype=float)
    vals = np.unique(sig)
    if vals.size != 2:
        raise SystemExit(
            f"cant 基序 {base!r} 的 σ 有 {vals.size} 个取值 ({vals.tolist()}), 需要恰好 2 个自旋"
            "亚晶格。单取值 (= 全同向 fm) 时 ±θ 倾斜只是整体自旋旋转, 系数对一切 θ 都等于 fm 行。")
    want = {"basis": np.asarray(basis, dtype=float),
            "stripe_a": sign_pow(n0), "stripe_b": sign_pow(n1)}[CANT_BASES[base]["pair"]]
    if sig.shape != want.shape or not np.array_equal(sig, want):
        raise SystemExit(
            f"cant 基序 {base!r} 的 σ 图样与声明不符 (预期 {CANT_BASES[base]['pair']}); "
            "倾斜的两个自旋亚晶格不是本脚本预期的对手 —— 拒绝生成, 否则 J1 的 θ 标度会静默错位。\n"
            "  提示: 只有「每个原胞恰好一上一下」的图样才能让 C 心键 (1/2,1/2) 跨亚晶格。")
    return sig


def parse_poscar(path):
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    scale = float(lines[1].strip())
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip().lower()
    idx += 1
    coords, atoms = [], []
    for sp, n in zip(species, counts):
        for _ in range(n):
            coords.append([float(x) for x in lines[idx].split()[:3]])
            atoms.append(sp)
            idx += 1
    coords = np.array(coords)
    if mode[:1] in "ck":
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    return lat, atoms, coords, mode


def write_poscar(path, lat, atoms, coords, comment="centered-rectangular supercell"):
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(sp) for sp in order]
    lines = [comment, "1.00000000000000"]
    for row in lat:
        lines.append(" ".join(f"{x:20.12f}" for x in row))
    lines.append("  ".join(order))
    lines.append("  ".join(str(c) for c in counts))
    lines.append("Direct")
    for sp in order:
        for a, c in zip(atoms, coords):
            if a == sp:
                lines.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_supercell(lat, atoms, coords, sc):
    sc = np.asarray(sc, dtype=int)
    new_atoms, new_coords, cell_idx = [], [], []
    for n0 in range(sc[0]):
        for n1 in range(sc[1]):
            for n2 in range(sc[2]):
                for a, c in zip(atoms, coords):
                    fc = (c + np.array([n0, n1, n2], dtype=float)) / sc
                    new_atoms.append(a)
                    new_coords.append(fc)
                    cell_idx.append((n0, n1, n2))
    new_lat = lat.copy()
    new_lat[0] *= sc[0]
    new_lat[1] *= sc[1]
    new_lat[2] *= sc[2]
    return new_lat, new_atoms, np.array(new_coords), np.array(cell_idx, dtype=int)


def magnetic_shells(idxs, coords, lat, max_shells=3, nshift=2, tol=0.05):
    """显式平移搜索识别磁性亚晶格前 max_shells 个距离壳。

    tol 默认 0.05 (不是矩形脚本的 0.2): C 心网的 C 心键与 a 键在 CrOCl 只差 0.163 A,
    用 0.2 会把两个物理上不同的壳合并成一个 J, 而且不会报错。

    距离是 **3D 键长** (含面外起伏)。注意 poscar_to_mc3.py 走的是另一套口径 —— 它把
    POSCAR 投到 2D 再搜索, 所以它把同一个 C 心键报成 2.506 A (面内投影), 而不是 3.0347 A。
    两套口径的**键拓扑相同** (每格点 4/2/2), 只是尺子不同; 别拿一方的距离去对另一方。
    """
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    ds = []
    n = len(idxs)
    for a in range(n):
        ci = coords[idxs[a]]
        for b in range(n):
            base = coords[idxs[b]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                ds.append(d)
    ds = sorted(ds)
    shells = []
    for d in ds:
        if not shells or d - shells[-1] > tol:
            shells.append(d)
    return shells[:max_shells]


def classify_bonds(idxs, coords, lat, sigma, shells, nshift=1, tol=0.05):
    """有向键和按壳归类, 返回与 shells 等长的列表。

    `sigma` 可以是两种口径, **机制完全相同**, 所以 cant 族不需要第二套代码、也不会漂移:
      * 一维 ±1 共线符号数组 (旧口径): 累加 sigma_i * sigma_j;
      * (n,3) **自旋矢量** (非共线口径): 累加 S_i·S_j。
    ⚠ 共线分支的表达式**一个字都没改** —— 这是"不给 --cant-theta 时逐位不变"的硬不变式。
    """
    sig = np.asarray(sigma, dtype=float)
    if sig.ndim == 1:
        if sig.size != len(idxs):
            raise SystemExit(f"σ 数组长度 {sig.size} != 磁原子数 {len(idxs)}; 拒绝静默错位。")
        dot = lambda i, j: float(sig[i] * sig[j])
    elif sig.shape == (len(idxs), 3):
        dot = lambda i, j: float(np.dot(sig[i], sig[j]))
    else:
        raise SystemExit(
            f"classify_bonds 的自旋参数形状 {sig.shape} 既不是 (n,) 也不是 (n,3) (n={len(idxs)});"
            " 共线请给 ±1 符号数组, cant 请给 mag_vectors.spin_cant 的三分量矢量。")
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    sums = [0.0] * len(shells)
    n = len(idxs)
    for i in range(n):
        ci = coords[idxs[i]]
        for j in range(n):
            base = coords[idxs[j]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                k = int(np.argmin([abs(d - s) for s in shells]))
                if abs(d - shells[k]) < tol:
                    sums[k] += dot(i, j)
    return sums


def sublattice_bond_stats(idxs, coords, lat, sigma, shells, nshift=1, tol=0.05):
    """每个壳里「同亚晶格 (σ_iσ_j=+1)」与「跨亚晶格 (σ_iσ_j=-1)」的有向键条数。

    用 σ 的**符号**分类 (而不是 σσ 的数值), 与 classify_bonds 的双重循环完全同构 ——
    同一套平移搜索与容差, 于是统计出来的条数与系数行严格自洽, 并直接给出 θ 标度的判据:
      * 同亚晶格壳 N_intra 条: c_k = -0.5*N_intra, **与 θ 无关** (σσ 恒 +1, sin²+cos²=1);
      * 跨亚晶格壳 N_inter 条: c_k(θ) = -0.5*N_inter*(sin²θ-cos²θ) = +0.5*N_inter*cos2θ
        (壳纯跨时 c_k(θ) = c_k(base)·cos2θ, 因为 c_k(base) = +0.5*N_inter)。
    这个三元组 (N_intra, N_inter) 也是"两个自旋亚晶格是不是预期那一对"的硬证据。
    """
    sig = np.asarray(sigma, dtype=float)
    if sig.ndim != 1 or sig.size != len(idxs):
        raise SystemExit("sublattice_bond_stats 需要与磁原子数等长的一维 σ 数组。")
    sgn = np.sign(sig)
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    intra = [0] * len(shells)
    inter = [0] * len(shells)
    n = len(idxs)
    for i in range(n):
        ci = coords[idxs[i]]
        for j in range(n):
            base = coords[idxs[j]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                k = int(np.argmin([abs(d - s) for s in shells]))
                if abs(d - shells[k]) < tol:
                    if sgn[i] * sgn[j] > 0:
                        intra[k] += 1
                    else:
                        inter[k] += 1
    return intra, inter


def _cross_check(spec_kwargs, coeff, n):
    """与 P0.1 规格层 (magnetic_orders.py) 交叉校验: 两套独立实现必须逐位一致。"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import magnetic_orders as MO

    sh_tol = spec_kwargs.get("shell_tol", 0.05)
    spec = MO.Spec(
        key="crect_selftest", label="C 心矩形 (selftest)",
        structure=MO.StructureSpec(
            poscar=spec_kwargs["poscar"], element=spec_kwargs["element"],
            nonmag=tuple(spec_kwargs["nonmag"]), supercell=tuple(spec_kwargs["supercell"])),
        shells=tuple(MO.Shell(f"J{i + 1}", d, tol=sh_tol)
                     for i, d in enumerate(spec_kwargs["shells"])),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: sigmas_of(
            c.n0, c.n1, np.where(c.prim_idx[c.mag_idx] == spec_kwargs["basis_prim"], 1.0, -1.0))[nm]))
            for nm in ORDERS),
        selected=tuple(SELECTED), spin=spec_kwargs["spin"],
        policy="shift_search", nshift=1,
    )
    res, ctx = MO.compute(spec)
    worst, where = 0.0, ""
    for nm in ORDERS:
        a = np.asarray(coeff[nm], dtype=float)
        b = np.asarray(res[nm].coeffs, dtype=float)
        if a.shape != b.shape:
            return False, f"{nm}: 形状不符 {a.shape} vs {b.shape}"
        err = float(np.max(np.abs(a - b)))
        if err > worst:
            worst, where = err, nm
    return worst <= 1e-9, f"最大逐位偏差 {worst:.3g} ({where})"


def _mag_vectors():
    """惰性 import 同目录的 mag_vectors (自旋矢量的**唯一真源**: spin_cant/render_magmom/...)。

    刻意不在模块顶部 import: 不给共线默认路径增加任何 import 副作用 —— "不给 --cant-theta
    时输出与旧版逐位相同"这条不变式不该依赖 import 顺序。同理也不自己实现 sinθ/cosθ:
    自造一套公式就会与 mag_vectors 的约定 (两亚晶格面内同向、z 反平行) 漂移。
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import mag_vectors as MV
    return MV


def _run_builder_sub(argv, timeout=900):
    """在**子进程**里跑本脚本, 返回 (returncode, 合并输出)。

    守卫与"双建字节比对"必须在子进程里做: 同进程 main() 会把文件写到一半就被 SystemExit 打断,
    也会污染当前 args.outdir。用 sys.executable 而不是字面的 "python" —— 本机 `python` 可能是
    Windows Store 的别名桩 (退出码 9009), 用它会得到与守卫无关的假失败。
    PYTHONIOENCODING=utf-8 + errors=replace: 中文/希腊字母提示在管道 (非控制台) 里不会因
    本地 codepage 编码失败而炸掉子进程。
    """
    cmd = [sys.executable, os.path.abspath(__file__)] + [str(x) for x in argv]
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    pr = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                        errors="replace", env=env, timeout=timeout)
    return pr.returncode, (pr.stdout or "") + (pr.stderr or "")


def _sha12(path):
    """文件 SHA-256 前 12 位 (自检打印的逐字节证据)。"""
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:12]


def _poscar_same_except_comment(p1, p2):
    """两个 POSCAR 除**第 1 行注释**外逐字节相同? (cant 序只该换注释里的序名)"""
    with open(p1, "rb") as f:
        a = f.read().split(b"\n", 1)
    with open(p2, "rb") as f:
        b = f.read().split(b"\n", 1)
    return len(a) == 2 and len(b) == 2 and a[1] == b[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar", nargs="?", default=None,
                    help="POSCAR 路径; 仅 --selftest 允许省略 (省略时用仓库自带 CrOCl 夹具, "
                         "并**显式打印**用了哪个文件)")
    ap.add_argument("--supercell", type=int, nargs=3, default=[4, 4, 1])
    ap.add_argument("--outdir", default=None,
                    help="输出目录 (默认 cases; 仅 `--selftest` 省略 POSCAR 的裸调用默认写到"
                         "临时目录, 并把实际路径打印出来)")
    ap.add_argument("--element", default="Cr")
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--moment", type=float, default=3.0)
    ap.add_argument("--spin", type=float, default=1.5)
    ap.add_argument("--shell-tol", type=float, default=0.05,
                    help="壳容差; 必须小于相邻**不等价**壳的最小间距 (默认 0.05)")
    # ---- 弱铁磁倾斜族 (P2′); 不给 --cant-theta 时一切照旧 ----
    ap.add_argument("--cant-theta", default=None,
                    help="cant 角度列表 (空格/逗号分隔, 度); 每个值生成 cant_p|t| 与 cant_m|t| 两个序。"
                         "不给这个开关就不生成倾斜族 (默认路径逐位不变)")
    ap.add_argument("--cant-azimuth", type=float, default=0.0,
                    help="cant 面内方位 φ (度, [0,360)); 传给 mag_vectors.spin_cant (默认 0.0)")
    ap.add_argument("--cant-base", default=CANT_BASE_DEFAULT,
                    help=f"cant 的共线基序, 取 {sorted(CANT_BASES)} (默认 {CANT_BASE_DEFAULT}); "
                         "fm 被显式拒绝 (单亚晶格 => ±θ 只是整体旋转, 该族零信息)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    # --outdir 的默认值分两种: 给了 POSCAR 就是历史默认 "cases"; 只有"裸 --selftest"
    # (没给 POSCAR) 才落到临时目录 —— 那次调用以前根本跑不起来 (argparse 直接退出),
    # 所以这不改变任何既有调用的落盘位置, 也不会往仓库里塞目录。
    if args.outdir is None:
        args.outdir = (os.path.join(tempfile.gettempdir(), "crect_selftest_cases")
                       if args.poscar is None else "cases")

    # ---- POSCAR 位置参数省略的处理 (绝不静默换输入) ----
    # 只有 --selftest 允许省略: 那时用仓库自带的 CrOCl 夹具, 并在任何计算之前把"到底用了
    # 哪个文件 / 哪个 --element --nonmag"打印出来。非 selftest 的省略一律硬报错。
    if args.poscar is None:
        if not args.selftest:
            raise SystemExit(
                "缺少 POSCAR 位置参数: 用法 `build_crect_cases.py <POSCAR> [选项]`; "
                "只有 --selftest 允许省略 (省略时用仓库自带的 CrOCl 夹具)。")
        fix = os.path.join(_REPO_ROOT, "tests", "fixtures", "CrOCl_POSCAR")
        if not os.path.isfile(fix):
            raise SystemExit(f"--selftest 省略 POSCAR 时需要仓库夹具 {fix}, 但它不存在; "
                             "请显式给出 POSCAR 路径。")
        if not args.nonmag:
            args.nonmag = ["Cl", "O"]
        args.poscar = fix
        print(f"note: --selftest 未给 POSCAR, 使用仓库夹具 {fix} "
              f"(element={args.element}, nonmag={args.nonmag}); 落盘目录 {args.outdir}")

    # ---- cant 参数的**静态**校验: 在任何 I/O 与重计算之前 fail fast ----
    # (守卫自检也在子进程里跑这些路径, 因此这里刻意不读 POSCAR / 不碰壳搜索)
    # 注意用 `is not None` 而不是真值判断: `--cant-theta ""` 必须是**硬报错**, 不能因为空串是
    # 假值就被静默当成"没请求 cant"(那就是静默降级)。
    cant_thetas = parse_cant_theta(args.cant_theta) if args.cant_theta is not None else []
    cant_plan = cant_order_plan(cant_thetas, args.cant_base) if cant_thetas else []
    if not (math.isfinite(args.cant_azimuth) and 0.0 <= args.cant_azimuth < 360.0):
        raise SystemExit(
            f"--cant-azimuth {args.cant_azimuth} 不合法: 要求**有限**且落在 [0, 360) 度内。\n"
            "  方位角是几何量 (面内方位), 刻意不做 mod 360 的静默归一 —— 那会把输入错误藏起来。")

    sc = [int(x) for x in args.supercell]
    if sc[0] < 2 or sc[1] < 2:
        raise SystemExit(
            f"--supercell {sc} 的前两维必须 >= 2: C 心键沿 a/b 都是半周期, 1 宽的超胞里 "
            "沿该方向的平移像会与原件重合。")

    lat, atoms, coords, mode = parse_poscar(args.poscar)  # 已把 Cartesian 归一为 Direct
    n_prim = len(atoms)
    cr_prim = [i for i, a in enumerate(atoms) if a == args.element]
    if cant_plan and len(cr_prim) != 2:
        # cant 专属提示放在通用守卫之前: 倾斜族的物理就是"两个亚晶格反平行 + 各自倾斜 ±θ",
        # 亚晶格数不等于 2 时 spin_cant 会把同一亚晶格的自旋也算成倾斜对 (静默错物理)。
        raise SystemExit(
            f"--cant-theta 要求模型恰好有 2 个磁性原子/胞 (两个自旋亚晶格), 本胞有 {len(cr_prim)} 个 "
            f"{args.element}; ±θ 倾斜的对象不存在。\n"
            "  提示: cant 族只对 C 心矩形 (2 个磁性原子/胞) 定义; 其他晶格请用对应的 builder。")
    if len(cr_prim) < 2:
        raise SystemExit(
            f"胞内只有 {len(cr_prim)} 个 {args.element} 原子; 本脚本要求 "
            "**2 个磁性原子/胞** 的 C 心网 (矩形/方形请改用 build_rectangular_cases.py)")
    if len(cr_prim) != 2:
        raise SystemExit(f"胞内有 {len(cr_prim)} 个 {args.element}, 本脚本只支持恰好 2 个")

    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    mag = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(mag)
    prim_of_mag = np.array([i % n_prim for i in mag], dtype=int)
    basis = np.where(prim_of_mag == cr_prim[0], 1.0, -1.0)

    # n0/n1 提前到这里 (纯整数切片, 数值上与原位置完全一致), 因为 cant 的**静态亚晶格守卫**
    # 必须在任何壳搜索/写盘之前 fail fast: 倾斜的对象错了, 后面的系数再漂亮也是错的。
    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]
    sig_base = None
    if cant_plan:
        sig_base = cant_sublattice_guard(basis, args.cant_base, basis, n0, n1)

    shells = magnetic_shells(mag, coords_sc, lat_sc, max_shells=3, tol=args.shell_tol)
    if len(shells) < 3:
        raise SystemExit(f"只找到 {len(shells)} 个磁距壳, 需要 3 (J1/J2/J3); 请加大 --supercell")
    d_nn, d_nnn, d_3rd = shells
    gaps = [d_nnn - d_nn, d_3rd - d_nnn]

    sigmas = sigmas_of(n0, n1, basis)

    sums = {k: classify_bonds(mag, coords_sc, lat_sc, s, shells,
                              tol=args.shell_tol) for k, s in sigmas.items()}
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums.items()}

    # 每格点的键数 (fm 序的系数 = -per_site/2): (4, 2, 2) 对 CrOCl。
    # 这个三元组是**跨口径不变量**: poscar_to_mc3.py 在 2D 投影里搜索, 距离读数是
    # 2.506/3.198/3.860, 但每格点键数必须同样是 (4, 2, 2) —— 驱动层就拿它做交叉校验。
    per_site_expected = [-2.0 * c / n for c in coeff["fm"]]

    A = np.array([[1.0] + coeff[k] for k in SELECTED])
    detA = float(np.linalg.det(A))
    if abs(detA) < 1e-6:
        raise SystemExit(
            f"选定的 4 个磁序设计矩阵奇异 (det={detA:.3g}); 拒绝反解。\n"
            "  提示: 只依赖 (n0,n1) 奇偶的序与它的 basis 翻倍版本必然逐位相同 "
            "(C 心矢量连接同一胞内的两个 basis), 不能同时选。")
    Ainv = np.linalg.inv(A)
    cond = float(np.linalg.cond(A))
    if cond > 1e4:
        raise SystemExit(
            f"选定的 4 个磁序设计矩阵病态: cond(A)={cond:.4g} > 1e4; 拒绝反解。\n"
            "  提示: 先算系数矩阵条件数再决定提交哪些构型 (LATTICE_ROADMAP §122)。")

    # ---------------- 弱铁磁倾斜族 (P2′): 系数一律从**真实三分量自旋**算 ----------------
    # 物理 (刻意写成闭式, 便于与数值结果对拍):
    #   S_i = (σ_i sinθ, 0, σ_i cosθ) —— 两个 z 反平行的亚晶格在 (x,z) 面内各自倾斜 ±θ,
    #   面内方位由 φ 给出 (φ=0 即上面那个 x 方向)。于是点积只有两种取值:
    #     同亚晶格 (σ_iσ_j=+1): sin²θ + cos²θ = 1       -> **与 θ 无关**
    #     跨亚晶格 (σ_iσ_j=-1): sin²θ - cos²θ = -cos2θ  -> 按 -cos2θ 标度
    #   => c_k(θ) = -0.5*[N_intra_k + N_inter_k*(-cos2θ)]
    #             = sin²θ*(-0.5*N_k) + cos²θ*(-0.5*sums[base]_k)
    #   即 c(θ) = sin²θ·z_k + cos²θ·c_base_k (z_k = fm 行 = 每格点键数)。
    #   这是本族**最重要的诚实结论**: 每一行都是 fm 行与基序行的仿射组合 => 倾斜族不给设计
    #   矩阵增加任何新秩 (能拟合同一组 {e0,J1,J2,J3}, 但分不开 DMI、也分不开 J1 与 J1′)。
    # 数值上刻意**不重写公式**: 复用 classify_bonds, 只是把 σ 换成 (n,3) 自旋矢量 ——
    #   共线分支的代码一个字都没动, 所以不给 --cant-theta 时输出逐位不变。
    cant_names, cant_theta, cant_coeff, cant_report = [], {}, {}, {}
    cant_intra_shells, cant_inter_shells = {}, {}
    if cant_plan:
        MV = _mag_vectors()
        # n_cells 在下面 (化学计量校验处) 才定义 —— 这里用**同一个表达式**提前算一遍, 让 cant
        # 的守卫与 cant_report 的每原胞口径都在同一段里, 不依赖后文定义位置 (值完全相同)。
        n_cells = int(sc[0] * sc[1] * sc[2])
        if n_cells < 1:
            raise SystemExit(f"--supercell {tuple(sc)} 不是合法胞数; cant 的每原胞口径无法定义。")
        SH = ("J1", "J2", "J3")
        base_row = [float(x) for x in coeff[args.cant_base]]
        intra3, inter3 = sublattice_bond_stats(mag, coords_sc, lat_sc, sig_base, shells,
                                              tol=args.shell_tol)
        cant_intra_shells = {k: int(intra3[i]) for i, k in enumerate(SH)}
        cant_inter_shells = {k: int(inter3[i]) for i, k in enumerate(SH)}
        print(f"cant 基序 {args.cant_base}: 同亚晶格壳 (系数与 θ 无关) "
              f"{[k for k in SH if cant_intra_shells[k]]} "
              f"{[cant_intra_shells[k] for k in SH]}, 跨亚晶格壳 (-cos2θ 标度) "
              f"{[k for k in SH if cant_inter_shells[k]]} "
              f"{[cant_inter_shells[k] for k in SH]}  (有向键条数, 与系数行同构统计)")
        for k, i, e in zip(SH, intra3, inter3):
            if i and e:
                raise SystemExit(
                    f"壳 {k} 同时含 {i} 条同亚晶格键与 {e} 条跨亚晶格键 (混合壳): 它的 θ 标度不再"
                    "是单一的 -cos2θ, cant 族在该壳上没有闭式预言, cant_report 的偏差就没法解释。\n"
                    "  提示: 换 --supercell / --shell-tol, 让每个壳按亚晶格类型分开。")
        for nm, t in cant_plan:
            sv = np.asarray(MV.spin_cant(sig_base, t, azimuth_deg=args.cant_azimuth), dtype=float)
            if sv.shape != (n, 3):
                raise SystemExit(
                    f"mag_vectors.spin_cant 返回 {sv.shape}, 期望 ({n}, 3) —— 每个磁原子一个"
                    "三分量自旋; 拒绝静默降级成共线。")
            if t == 0.0 and not np.array_equal(sv[:, 2], sig_base):
                raise SystemExit(
                    "θ=0 时 spin_cant 的 z 分量与共线 σ 不逐位相等 —— mag_vectors 的契约被破坏"
                    "(spin_cant(...,0) 必须退化为 spin_collinear); 拒绝生成,"
                    " 否则 θ=0 的回归比对失去意义。")
            row = [float(-0.5 * s) for s in
                   classify_bonds(mag, coords_sc, lat_sc, sv, shells, tol=args.shell_tol)]
            if t != 0.0 and row == base_row:
                raise SystemExit(
                    f"cant 序 {nm} (θ={t:g}) 的系数行与基序 {args.cant_base} 逐位相同 —— "
                    "该行不携带任何新信息, 提交它只是浪费机时 (典型原因: 基序是全同向 fm,"
                    " ±θ 对整个自旋系只是一次旋转)。\n"
                    "  提示: --cant-base afm_c (C 心反平行图样, J1 跨亚晶格)。")
            if t < 0:
                tw = cant_order_name(-t)
                if tw in cant_coeff and cant_coeff[tw] != row:
                    raise SystemExit(
                        f"θ={t:g} 与 θ={-t:g} 的系数不逐位相同 —— cos2θ 是偶函数, 能量的"
                        " θ 符号本来就不可辨识; 这里只可能是自旋矢量构造出了非物理的 θ 不对称。")
            th_r = math.radians(t)
            z = [float(x) for x in coeff["fm"]]
            pred = [math.sin(th_r) ** 2 * z[i] + math.cos(th_r) ** 2 * base_row[i]
                    for i in range(3)]
            dev = max(abs(row[i] - pred[i]) for i in range(3))
            cant_names.append(nm)
            cant_theta[nm] = float(t)
            cant_coeff[nm] = row
            # 键名冻结: `base` / `coefficients_per_cell` / 闭式预测 / max_abs_dev(iation) 是
            # 下游求解器 (solve_J_crect_cant.py) 直接读的拼写, 一个都不许改名。
            # ⚠ 口径陷阱: 本仓库历史命名里的 `coefficients_per_cell` 实际是**每超胞** ——
            #   cant_report 里同名键沿用同口径 (与顶层既有键、与求解器 canted_coeff 完全自洽);
            #   **每原胞**的量一律显式叫 `*_per_primitive_cell` (已除以 n_cells)。
            pred_pc = [x / n_cells for x in pred]
            cant_report[nm] = dict(
                order=nm, name=nm, kind="cant", theta_deg=float(t),
                azimuth_deg=float(args.cant_azimuth),
                base=args.cant_base,
                n_cells=int(n_cells), n_mag_cell=2,
                coefficients_per_cell=list(row),
                coeff_per_primitive_cell=[x / n_cells for x in row],
                coeff_fm_per_primitive_cell=[x / n_cells for x in z],
                coeff_base_per_primitive_cell=[x / n_cells for x in base_row],
                prediction_per_cell=list(pred),
                prediction_per_primitive_cell=pred_pc,
                prediction_rule="c_k(theta) = sin^2(theta)*c_k(fm) + cos^2(theta)*c_k(base) "
                                "(每胞; 也等于 sin²θ*(z_k/2) + cos²θ*c_base_k 的每格点口径)",
                closed_form=dict(
                    formula="c(theta) = sin^2(theta)*z + cos^2(theta)*c_base",
                    z_per_cell=list(z), c_base_per_cell=list(base_row),
                    prediction_per_cell=list(pred),
                    prediction_per_primitive_cell=pred_pc,
                    max_abs_deviation=float(dev), max_abs_dev=float(dev)),
                max_abs_dev=float(dev), max_abs_deviation=float(dev),
                coefficients_per_cell_base=list(base_row),
                per_site_bonds={k: float(coeff["fm"][i] / n) for i, k in enumerate(SH)},
                per_shell_intra_bonds={k: cant_intra_shells[k] for k in SH},
                per_shell_inter_bonds={k: cant_inter_shells[k] for k in SH},
                intra_bonds_per_cell={k: cant_intra_shells[k] / n_cells for k in SH},
                inter_bonds_per_cell={k: cant_inter_shells[k] / n_cells for k in SH},
                shell_theta_scaling={k: ("mixed" if (cant_intra_shells[k] and
                                                     cant_inter_shells[k])
                                         else ("theta_independent" if cant_intra_shells[k]
                                               else "cos2theta")) for k in SH},
                theta_dependent_shells=[k for k in SH if cant_inter_shells[k]
                                        and not cant_intra_shells[k]],
                theta_independent_shells=[k for k in SH if cant_intra_shells[k]
                                          and not cant_inter_shells[k]],
                soc_required=bool(MV.needs_soc(nm)), requires_soc=bool(MV.needs_soc(nm)),
                sign_of_theta_unidentifiable=True,
                azimuth_not_determined_by_geometry=True,
                note_theta_sign="cos2θ 是偶函数 => 系数逐位满足 c(+θ)=c(-θ): θ 的**符号**不能由"
                                "总能量确定, 只有 |θ| (等价 cos2θ) 可辨识; cant_p/cant_m 两份"
                                "能量的差别只能来自数值噪声, 不能当物理。",
                note_azimuth="方位角 φ 由构建者选定 (本族 φ=%g°), 几何不决定它; 各向同性"
                             "Heisenberg 下能量只通过 cos2θ 依赖倾斜, 与 φ 无关 —— 要定 φ 必须"
                             "引入 SOC/各向异性 (needs a real SOC DFT run)。" % args.cant_azimuth,
                note="每一行 c(θ)=sin²θ·z+cos²θ·c(base) 都是 fm 行 (z=每格点键数) 与基序行的"
                     "**仿射组合** => 倾斜族不给设计矩阵增加任何新秩: 能拟合同一组 "
                     "{e0,J1,J2,J3}, 但**分不开 DMI 项 D, 也分不开 J1 与 J1′**; θ 的符号不可"
                     "辨识, 方位 φ 不由几何决定。",
            )
        for nm in cant_names:
            d = cant_report[nm]
            print(f"cant {nm}: θ={d['theta_deg']:+g}°  φ={d['azimuth_deg']:g}°  "
                  f"系数/胞 {['%.6f' % x for x in d['coefficients_per_cell']]}  "
                  f"(基序 {d['base']} {['%.6f' % x for x in d['coefficients_per_cell_base']]})  "
                  f"闭式偏差 {d['max_abs_deviation']:.3e}  SOC={d['requires_soc']}")

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}

    # ---------------- 非磁元素占位 + 化学计量校验 (oxyhalide 的 O/X 占位契约) ----------------
    # 顺序刻意: 先查「声明了但胞里没有」(多半是元素符号写错), 再查「胞里有但没声明」
    # (多半是复制了别人的 POSCAR 忘了改 --nonmag)。两种都会让 MAGMOM/占位静默错位。
    _declared = set(args.nonmag)
    _missing = [el for el in args.nonmag if nonmag_counts[el] == 0]
    if _missing:
        raise SystemExit(
            f"--nonmag 里的 {_missing} 在该胞里一个原子都没有; 检查元素符号是否写错。\n"
            f"  胞内实际元素: {sorted(set(atoms_sc))}")
    _undeclared = sorted(set(atoms_sc) - {args.element} - _declared)
    if _undeclared:
        raise SystemExit(
            f"POSCAR 里有未声明的元素 {_undeclared}; 请用 --nonmag 显式列出。\n"
            "  本脚本要求每个非磁元素都被声明, 否则占位与化学式会静默错位。")
    n_cells = int(sc[0] * sc[1] * sc[2])
    _per_cell_float = {el: cnt / n_cells for el, cnt in
                       [(args.element, n)] + list(nonmag_counts.items())}
    _nonint = {el: v for el, v in _per_cell_float.items()
               if abs(v - round(v)) > 1e-9}
    if _nonint:
        raise SystemExit(
            f"化学计量不是整数: 每胞 {_nonint} (超胞 {n_cells} 个胞)。\n"
            "  提示: 原胞选择或 --supercell 倍数不对, 或者 --element 选错了。")
    per_cell = {el: int(round(v)) for el, v in _per_cell_float.items()}
    _formula = "".join(f"{el}{per_cell[el] if per_cell[el] != 1 else ''}"
                       for el in [args.element] + list(args.nonmag))
    print(f"化学计量 (每胞): {_formula}  {per_cell}")
    # MAGMOM 必须与 POSCAR 的**元素分组顺序**一致 (write_poscar 按首次出现顺序分组)。
    # 不能假设磁性元素排在第一位: CrOCl 的 JARVIS 原胞就是 `Cl Cr O` 开头,
    # 早期版本把 Cr 的矩写在最前面 -> MAGMOM 与原子错位, 而 DFT 照样能跑完。
    sp_order, blocks = [], []
    for a in atoms_sc:
        if a not in sp_order:
            sp_order.append(a)
    for sp in sp_order:
        blocks.append((sp, [i for i, a in enumerate(atoms_sc) if a == sp]))
    pos_in_mag = {i: k for k, i in enumerate(mag)}
    for name in ORDERS:
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"{args.element} centered-rectangular {name}")
        vals = []
        for sp, idxs in blocks:
            if sp == args.element:
                vals += [f"{args.moment * sigmas[name][pos_in_mag[i]]:.1f}" for i in idxs]
            else:
                vals += ["0.0"] * len(idxs)
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(vals)}\n")

    # ---------------- cant 序落盘: 每原子**三分量**, 完全展开 (不用 N*val 压缩) ----------------
    # 为什么必须三分量 + 完全展开: VASP 的 MAGMOM 三列就是 (mx,my,mz); 写成 N*0.0 之类的压缩
    #   形式会让"每个原子三个数"的对应关系断掉 (而 LSORBIT 下非共线是**逐原子**三分量)。
    #   非磁原子给三个 0.0 (它们的自旋自由度不参与模型, 但 MAGMOM 的条数必须与 POSCAR 一致)。
    # 幅值沿用既有 --moment 约定: 单位矢量 × moment (自旋矢量是单位长度的方向, 见 mag_vectors)。
    soc_orders, nc_text = [], {}
    if cant_plan:
        sig_cant = sigmas_of(n0, n1, basis)[args.cant_base]
        for nm in cant_names:
            d = os.path.join(args.outdir, nm)
            os.makedirs(d, exist_ok=True)
            write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                         comment=f"{args.element} centered-rectangular {nm}")
            sv = np.asarray(MV.spin_cant(sig_cant, cant_theta[nm],
                                         azimuth_deg=args.cant_azimuth), dtype=float)
            values = MV.magmom_values(sv, args.moment, mag, len(atoms_sc))
            text = MV.render_magmom(values)
            # 落盘前**自己先校验**: 值个数/分量数/有磁矩原子数/|m| 全部必须过 ——
            #   (noncollinear => want_ncomp=3; 非磁原子的三个 0.0 由 spin_all 填)
            ok, ncomp, problems = MV.check_magmom(text, len(atoms_sc), moment=args.moment,
                                                  want_ncomp=3, n_mag=n)
            if not ok:
                raise SystemExit(
                    f"cant 序 {nm} 的 magmom.txt 自校验失败 ({ncomp} 分量): {problems}; "
                    "拒绝落盘 (宁可不生成, 也不要一个 VASP 会读错的 MAGMOM)。")
            nc_text[nm] = text
            with open(os.path.join(d, "magmom.txt"), "w", encoding="utf-8") as f:
                f.write(text)
            if MV.needs_soc(nm):
                soc_orders.append(nm)

    # order_spec: 既有 7 个共线条目**逐字不变** (kind="collinear"), cant 条目追加在后面 ——
    # 设计矩阵只用 SELECTED 里的 4 个共线序, 所以追加 cant 不会动反解结果。
    order_spec = [dict(name=nm, kind="collinear", rule=RULES[nm],
                       coeffs=[float(x) for x in coeff[nm]]) for nm in ORDERS]
    for nm in cant_names:
        order_spec.append(dict(name=nm, kind="cant", theta_deg=cant_theta[nm],
                               azimuth_deg=float(args.cant_azimuth),
                               base=args.cant_base, rule=CANT_RULE,
                               coeffs=[float(x) for x in cant_coeff[nm]],
                               requires_soc=bool(cant_report[nm]["requires_soc"])))
    model_spec = dict(
        schema="model_spec/1", lattice="rect_centered",
        label="C 心矩形多亚晶格 (2 磁性原子/胞)",
        element=args.element, nonmag=list(args.nonmag), spin=args.spin, S2=args.spin ** 2,
        supercell=sc, policy="shift_search", nshift=1,
        noncollinear=bool(cant_names), n_mag=n, n_basis_mag=2,
        shells=[dict(name="J1", d=d_nn, tol=args.shell_tol, pair=None),
                dict(name="J2", d=d_nnn, tol=args.shell_tol, pair=None),
                dict(name="J3", d=d_3rd, tol=args.shell_tol, pair=None)],
    )
    # cant 专属的 model_spec 键**只在真生成了 cant 序时才加** —— 保证不给 --cant-theta 时
    # model_spec (以及整个共线部分) 与旧版逐位相同 (CANT_NEW_MODEL_KEYS)。
    if cant_names:
        model_spec["soc_orders"] = list(soc_orders)
        model_spec["cant_base"] = args.cant_base

    # coefficients_per_cell / _per_mag 在 cant 打开时**追加** cant 行 (requirement: 必须含 cant);
    # 关闭时两者内容与插入顺序都与旧版逐位相同 (dict(coeff) 再 update 的顺序保证)。
    all_coeff = dict(coeff)
    all_coeff.update(cant_coeff)

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        # 每胞化学计量 (P2 契约: O/X 占位与化学式必须自洽, 见 --nonmag 校验)
        "stoichiometry_per_cell": per_cell,
        "formula_per_cell": _formula,
        "n_cells": n_cells,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shell_tol": args.shell_tol,
        "shells_angstrom": {"J1": d_nn, "J2": d_nnn, "J3": d_3rd},
        # 跨口径不变量: 每格点键数。poscar_to_mc3.py 在 2D 投影里工作, 它的壳距离
        # 读到 2.506/3.198/3.860 (面内) 而不是 3.0347/3.198/3.860 (3D 键长), 所以
        # 驱动层不能拿距离去对, 只能拿这个键数三元组去对。
        "per_site_bonds": {"J1": per_site_expected[0], "J2": per_site_expected[1],
                           "J3": per_site_expected[2]},
        "selected_orders": SELECTED,
        "all_orders": ORDERS,
        "noncollinear": bool(cant_names),
        "coefficients_per_cell": {k: list(v) for k, v in all_coeff.items()},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in all_coeff.items()},
        "inversion_matrix_per_cell": Ainv.tolist(),
        "design_det": detA,
        "design_cond": cond,
        "order_spec": order_spec,
        "model_spec": model_spec,
        # ---- P2′ cant 族的附加键 (全部向后兼容; 求解器按这些**确切拼写**读) ----
        # 不给 --cant-theta 时: cant_orders=[] / cant_theta_deg={} / cant_report={} /
        #   soc_orders=[] / cant_azimuth_deg=0.0 / noncollinear=false。
        "cant_orders": list(cant_names),
        "cant_theta_deg": {k: float(v) for k, v in cant_theta.items()},
        "cant_azimuth_deg": float(args.cant_azimuth),
        "cant_report": cant_report,
        # 需要 SOC (=> 下游必须 --lsorbit) 的序: cant 族全部在其中 (needs_soc("cant") 为真)。
        # 顺带把共线序里同样需要 SOC 的也带出来 (今天为空), 免得下游各自判断。
        "soc_orders": list(soc_orders),
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"模型: E = e0 + c1*J1~ + c2*J2~ + c3*J3~  (J~=J*S^2, S={args.spin})")
    print(f"超胞 {tuple(sc)} -> {len(atoms_sc)} 原子 ({n} {args.element}), nonmag {nonmag_counts}")
    print(f"距离壳 (3D 物理键长): J1={d_nn:.4f} A (C 心 (1/2,1/2) 键, 连同胞两磁原子)"
          f"  J2={d_nnn:.4f} A (a)  J3={d_3rd:.4f} A (b)")
    print(f"每格点键数 (跨口径不变量, 驱动层据此校验 MC 模型): "
          f"J1={per_site_expected[0]:.0f}  J2={per_site_expected[1]:.0f}  "
          f"J3={per_site_expected[2]:.0f}")
    print(f"  注: poscar_to_mc3.py 在**面内投影**里搜索, 同一个 J1 键它读作 2.506 A "
          f"(= 3.0347 A 去掉 1.71 A 面外起伏) —— 尺子不同, 键数相同。")
    print(f"相邻壳间距: J2-J1={gaps[0]:.4f} A   J3-J2={gaps[1]:.4f} A   (壳容差 {args.shell_tol})")
    if min(gaps) < 0.3:
        print(f"[warn] 相邻壳间距 {min(gaps):.4f} A < 0.3 A —— 用 tol=0.2 的旧脚本会把这两个"
              "物理上不同的壳静默合并成一个 J; 本脚本已按 --shell-tol 分开。")
    print(f"选定 {SELECTED}, cond(A)={cond:.4g}, det(A)={detA:.6g}")
    for k in ORDERS:
        tag = " *" if k in SELECTED else "  "
        print(f" {tag}{k:10s} c1={coeff[k][0]:8.3f} c2={coeff[k][1]:8.3f} c3={coeff[k][2]:8.3f}"
              f"   (per mag: {coeff[k][0] / n:+.3f}, {coeff[k][1] / n:+.3f}, {coeff[k][2] / n:+.3f})")
    for k in cant_names:
        d = cant_report[k]
        print(f"  cant {k:10s} c1={d['coefficients_per_cell'][0]:8.3f} "
              f"c2={d['coefficients_per_cell'][1]:8.3f} "
              f"c3={d['coefficients_per_cell'][2]:8.3f}   (θ={d['theta_deg']:+g}°, "
              f"φ={d['azimuth_deg']:g}°, 基序 {d['base']}, SOC={d['requires_soc']})")
    if cant_names:
        print(f"cant 族 {cant_names}: 需要 SOC (--lsorbit) 的序 = {soc_orders}")
        print("  诚实边界: 每一行 c(θ)=sin²θ·z+cos²θ·c(base) 都是 fm 行与基序行的仿射组合 -> "
              "倾斜族**不给设计矩阵增加新秩** (能拟合同一组 e0/J1/J2/J3, 但分不开 DMI, "
              "也分不开 J1 与 J1′); θ 的符号不可辨识, 方位 φ 不由几何决定。")

    if args.selftest:
        bad = []

        # 1) 核心结构事实: basis 翻倍版本与它的 (n0,n1)-奇偶母序逐位相同
        for a_nm, b_nm in BASIS_DEGENERATE:
            dg = max(abs(x - y) for x, y in zip(coeff[a_nm], coeff[b_nm]))
            print(f"\n--selftest: {a_nm} vs {b_nm} 简并度 = {dg:.3g} "
                  f"{'PASS' if dg < 1e-9 else 'FAIL'}")
            if dg >= 1e-9:
                bad.append(f"{a_nm} != {b_nm}; C 心矢量连接同一胞内两个 basis, 两者应当简并")

        # 2) fm 与 afm_c 只在 c1 上不同 —— c1 的分辨完全靠这一对
        dc = [abs(coeff["fm"][m] - coeff["afm_c"][m]) for m in range(3)]
        ok = dc[0] > 1e-9 and dc[1] < 1e-9 and dc[2] < 1e-9
        print(f"--selftest: fm 与 afm_c 只差 c1 (|dc|={['%.3g' % x for x in dc]}) "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append("fm 与 afm_c 不再只差 c1")

        # 3) 奇异组合必须被识别出来 (afm_a 与 afm_c_a 相同 => 同选必奇异)
        Asing = np.array([[1.0] + coeff[k] for k in ("fm", "afm_c", "afm_a", "afm_c_a")])
        dsing = abs(float(np.linalg.det(Asing)))
        print(f"--selftest: (fm,afm_c,afm_a,afm_c_a) 行列式 = {dsing:.3g} "
              f"{'PASS' if dsing < 1e-6 else 'FAIL'} (须为 0: 该组合不可选)")
        if dsing >= 1e-6:
            bad.append("退化组合 (fm,afm_c,afm_a,afm_c_a) 竟然非奇异")

        # 4) 七个磁序必须两两不同 (防命名/切片写重复)
        keys = [tuple(sigmas[k].tolist()) for k in ORDERS]
        uniq = len(set(keys))
        print(f"--selftest: 磁序数组互异 {uniq}/{len(ORDERS)} "
              f"{'PASS' if uniq == len(ORDERS) else 'FAIL'}")
        if uniq != len(ORDERS):
            bad.append(f"磁序数组有重复: 只有 {uniq} 个不同 (期望 {len(ORDERS)})")

        # 5) 合成 J -> 反解 (验证 A 可逆 + 系数自洽)
        J_true = [5.0, -0.5, 0.2]
        Jt = [j * args.spin ** 2 for j in J_true]
        E = {k: sum(coeff[k][m] * Jt[m] for m in range(3)) for k in SELECTED}
        Evec = np.array([E[k] for k in SELECTED]) / 1000.0
        sol = Ainv @ Evec
        got = [sol[1 + m] / args.spin ** 2 * 1000.0 for m in range(3)]
        inv_ok = all(abs(g - t) < 1e-6 for g, t in zip(got, J_true))
        print(f"--selftest: 合成 J={J_true} -> 反解 {[f'{x:.6f}' for x in got]}  "
              f"{'PASS' if inv_ok else 'FAIL'}")
        if not inv_ok:
            bad.append("合成 J 反解不闭合")

        # 6) 4x4x1 下的解析理想值 (per mag): 每胞 4 根 C 心键 / 2 根 a 键 / 2 根 b 键
        if tuple(sc) == (4, 4, 1):
            ideal = {"fm": (-2.0, -1.0, -1.0), "afm_c": (+2.0, -1.0, -1.0),
                     "afm_a": (0.0, +1.0, -1.0), "afm_b": (0.0, -1.0, +1.0),
                     "afm_diag": (0.0, +1.0, +1.0), "afm_c_a": (0.0, +1.0, -1.0),
                     "afm_c_diag": (0.0, +1.0, +1.0)}
            err = max(abs(coeff[k][m] / n - ideal[k][m])
                      for k in ideal for m in range(3))
            print(f"--selftest: 4x4x1 解析理想值最大偏差 = {err:.3g} "
                  f"{'PASS' if err < 1e-9 else 'FAIL'}")
            if err >= 1e-9:
                bad.append(f"4x4x1 系数偏离解析理想值 (最大 {err:.3g})")

        # 7) 与 P0.1 规格层交叉校验 (两套独立实现必须逐位一致)
        try:
            ok, msg = _cross_check(
                dict(poscar=os.path.abspath(args.poscar), element=args.element,
                     nonmag=args.nonmag, supercell=sc, shells=list(shells),
                     spin=args.spin, shell_tol=args.shell_tol,
                     basis_prim=int(cr_prim[0])), coeff, n)
            print(f"--selftest: 与 P0.1 规格层交叉校验 {'PASS' if ok else 'FAIL'} ({msg})")
            if not ok:
                bad.append(f"与规格层不一致: {msg}")
        except ImportError as exc:
            print(f"--selftest: 跳过与规格层交叉校验 (无法 import magnetic_orders: {exc})")

        # 8) 条件数闸门
        print(f"--selftest: cond(A)={cond:.4g} < 1e4 {'PASS' if cond < 1e4 else 'FAIL'}")

        # ---------------- 9) P2′ cant 族 ----------------
        # 分两半: 前半在**本进程**里核对物理与闭式 (快); 后半在**子进程**里做双建字节回归
        # 和守卫检查 (真实退出码/真实报错文本 —— 同进程断言不出这两件事)。
        # 没给 --cant-theta 时**不能什么都不查**: 那正是"默认路径逐位不变"的守门人, 所以这里
        # 用样本角在子进程里把 cant 自检整条跑一遍, 并把本进程产物与子进程 no-cant 产物逐字节比。
        if not cant_plan:
            print(f"--selftest: 9) cant 族未请求 (没给 --cant-theta): 共线路径不应受影响; "
                  f"用样本角 '{CANT_SELFTEST_SAMPLE}' 在子进程里真跑一遍 cant 自检 + 双建字节回归")
            d9 = tempfile.mkdtemp(prefix="crect_nocant_regress_")
            base9 = [os.path.abspath(args.poscar), "--element", args.element,
                     "--nonmag", *args.nonmag, "--supercell", *[str(x) for x in sc],
                     "--moment", str(args.moment), "--spin", str(args.spin),
                     "--shell-tol", str(args.shell_tol)]
            d9_no, d9_ca = os.path.join(d9, "no"), os.path.join(d9, "ca")
            rc9, out9 = _run_builder_sub(base9 + ["--outdir", d9_no])
            rc9c, out9c = _run_builder_sub(base9 + ["--outdir", d9_ca, "--cant-theta",
                                                    CANT_SELFTEST_SAMPLE, "--selftest"])
            ok_no = (rc9 == 0)
            ok_ca = (rc9c == 0 and "ALL PASS" in out9c)
            print(f"--selftest:   [{'PASS' if ok_no else 'FAIL'}] no-cant 子进程构建成功 (rc={rc9})")
            if not ok_no:
                bad.append(f"no-cant 子进程构建失败 rc={rc9}\n{out9[-400:]}")
            print(f"--selftest:   [{'PASS' if ok_ca else 'FAIL'}] "
                  f"--cant-theta '{CANT_SELFTEST_SAMPLE}' --selftest 子进程 ALL PASS (rc={rc9c})")
            if not ok_ca:
                bad.append(f"cant 自检子进程失败 rc={rc9c}\n{out9c[-800:]}")
            if ok_no:
                # 本进程 (无 cant) 的整份 JSON 与子进程 no-cant 的必须**整文件逐字节**相同 ——
                # 这是"不给开关时输出与旧版一致"的最强形式。
                j1 = os.path.join(args.outdir, "magnetic_orders.json")
                j2 = os.path.join(d9_no, "magnetic_orders.json")
                same_json = (os.path.isfile(j1) and os.path.isfile(j2)
                             and _sha12(j1) == _sha12(j2))
                print(f"--selftest:   [{'PASS' if same_json else 'FAIL'}] "
                      f"本进程 JSON 与 no-cant 子进程 JSON 整文件逐字节相同 "
                      f"sha12={_sha12(j1) if os.path.isfile(j1) else 'n/a'}")
                if not same_json:
                    bad.append("本进程 magnetic_orders.json 与 no-cant 子进程不同 (默认路径被改动)")
                diff = [nm for nm in ORDERS
                        if _sha12(os.path.join(args.outdir, nm, "magmom.txt"))
                        != _sha12(os.path.join(d9_no, nm, "magmom.txt"))]
                print(f"--selftest:   [{'PASS' if not diff else 'FAIL'}] "
                      f"7 个共线 magmom.txt 与本进程逐字节相同  不同={diff}")
                if diff:
                    bad.append(f"共线 magmom.txt 与 no-cant 子进程不同: {diff}")
            print(f"--selftest:   scratch 目录 {d9}")
        else:
            MV = _mag_vectors()
            SH3 = ("J1", "J2", "J3")

            def ck9(name, cond, extra=""):
                print(f"--selftest:   [{'PASS' if cond else 'FAIL'}] {name}"
                      + (f"  {extra}" if extra else ""))
                if not cond:
                    bad.append("cant 检查失败: " + name + (f"  [{extra}]" if extra else ""))

            # (9.1) θ 与 -θ 系数**逐位**相同 (cos2θ 是偶函数 => 能量的 θ 符号不可辨识)
            for nm in [k for k in cant_names if k.startswith("cant_p")]:
                tw = "cant_m" + nm[len("cant_p"):]
                if tw in cant_coeff:
                    a, b = cant_coeff[nm], cant_coeff[tw]
                    ck9(f"θ 符号不可辨识: {nm} 与 {tw} 系数逐位相同",
                        all(a[i] == b[i] for i in range(3)),
                        f"{['%.9f' % x for x in a]} vs {['%.9f' % x for x in b]}")

            # (9.2) 该族不是死重: θ≠0 的系数行必须与基序行**不同** (否则提交它纯属浪费机时)
            base_row_ = coeff[args.cant_base]
            for nm in cant_names:
                if cant_theta[nm] == 0.0:
                    continue
                # "行不同" 必须带容差看: 同亚晶格键的点积是 sin²θ+cos²θ (不是精确的 1.0),
                # 所以 θ 无关壳的系数只在 ~1e-14 (N_bonds 项浮点累加) 内不变, 不是逐位相同。
                big = [SH3[i] for i in range(3)
                       if abs(cant_coeff[nm][i] - base_row_[i]) > 1e-9]
                dep = cant_report[nm]["theta_dependent_shells"]
                resid = max([abs(cant_coeff[nm][i] - base_row_[i]) for i in range(3)
                             if SH3[i] not in dep] or [0.0])
                ck9(f"{nm} (θ={cant_theta[nm]:g}°) 的系数行不同于基序 {args.cant_base} 行, "
                    f"且差异壳 = 跨亚晶格壳",
                    bool(big) and big == dep,
                    f"差异壳={big} 跨亚晶格壳={dep} 其余壳残差={resid:.2e} (浮点累加, 非逐位)")

            # (9.3) (b) θ=0 必须**精确**退化为共线结果 (spin_cant(...,0) == spin_collinear)
            sig_b = sigmas_of(n0, n1, basis)[args.cant_base]
            row0 = [float(-0.5 * s) for s in classify_bonds(
                mag, coords_sc, lat_sc, MV.spin_cant(sig_b, 0.0, azimuth_deg=args.cant_azimuth),
                shells, tol=args.shell_tol)]
            ck9("θ=0 的系数与共线基序**逐位**相同 (spin_cant(...,0) ≡ spin_collinear)",
                all(row0[i] == base_row_[i] for i in range(3)),
                f"{['%.9f' % x for x in row0]} vs {['%.9f' % x for x in base_row_]}")

            # (9.4) (c) 闭式 θ 标度 c(θ)=sin²θ·z+cos²θ·c_base 的最大偏差必须 < 1e-12
            #         (残差只是 N_bonds 项浮点累加噪声, 会随超胞变大, 见 report 的 max_abs_dev)
            devs = {nm: float(cant_report[nm]["max_abs_dev"]) for nm in cant_names}
            n_ang = len({abs(cant_theta[nm]) for nm in cant_names})
            ck9(f"闭式 θ 标度偏差 < 1e-12 ({n_ang} 个不同 |θ|)",
                max(devs.values()) < 1e-12 and n_ang >= 2,
                f"max={max(devs.values()):.3e} " + str({k: '%.2e' % v for k, v in devs.items()}))

            # (9.5) (d) 壳归属 + 标度: 同亚晶格壳 θ 无关; 跨亚晶格壳 c_k(θ)=+0.5*N_inter*cos2θ
            for nm in cant_names:
                r = cant_report[nm]
                th = cant_theta[nm]
                c2 = math.cos(2.0 * math.radians(th))
                intra_dev = inter_dev = 0.0
                ratio = []
                for i, k in enumerate(SH3):
                    if k in r["theta_independent_shells"]:
                        intra_dev = max(intra_dev, abs(cant_coeff[nm][i] - base_row_[i]))
                    if k in r["theta_dependent_shells"]:
                        inter_dev = max(inter_dev, abs(cant_coeff[nm][i]
                                                       - 0.5 * r["per_shell_inter_bonds"][k] * c2))
                        ratio.append(cant_coeff[nm][i] / base_row_[i])
                ck9(f"{nm}: 同亚晶格壳 {r['theta_independent_shells']} 与 θ 无关 / "
                    f"跨亚晶格壳 {r['theta_dependent_shells']} 按 -cos2θ 标度",
                    intra_dev < 1e-12 and inter_dev < 1e-9
                    and all(abs(x - c2) < 1e-12 for x in ratio),
                    f"θ 无关偏差={intra_dev:.2e}  -cos2θ 偏差={inter_dev:.2e}  "
                    f"c(θ)/c(base)={['%.12f' % x for x in ratio]} vs cos2θ={c2:.12f}")
            # 本晶格 (C 心矩形) 的壳归属必须是 J1 跨 / J2,J3 同 —— 物理上唯一能把 J1 的 θ 标度
            # 与 J2/J3 分开的基序 (afm_c)。换个 POSCAR 就会变, 所以只在数值对上时断言。
            if cant_report[cant_names[0]]["base"] == "afm_c" \
                    and cant_intra_shells == {"J1": 0, "J2": 2 * n, "J3": 2 * n} \
                    and cant_inter_shells == {"J1": 4 * n, "J2": 0, "J3": 0}:
                ck9("本晶格壳归属: J1 (C 心 (1/2,1/2) 键) 跨亚晶格 => θ 相关; J2/J3 同亚晶格 => θ 无关",
                    True, f"有向键条数 intra={list(cant_intra_shells.values())} "
                          f"inter={list(cant_inter_shells.values())} (2*n={2 * n}, 4*n={4 * n})")
            else:
                print(f"--selftest:   [skip] 本晶格的壳归属不是 J1 跨/ J2,J3 同 "
                      f"(intra={list(cant_intra_shells.values())}, "
                      f"inter={list(cant_inter_shells.values())}) —— 只按上面的标度关系检查")

            # (9.6) (e) 三分量 MAGMOM 必须过 want_ncomp=3, 并被 1 分量检查**拒绝**
            for nm in cant_names:
                ok3, nc3, pb3 = MV.check_magmom(nc_text[nm], len(atoms_sc), moment=args.moment,
                                                want_ncomp=3, n_mag=n)
                ok1, nc1, pb1 = MV.check_magmom(nc_text[nm], len(atoms_sc), moment=args.moment,
                                                want_ncomp=1, n_mag=n)
                ck9(f"{nm}: magmom.txt 过 3 分量校验 / 被 1 分量校验拒绝",
                    ok3 and (not ok1) and nc3 == 3 and nc1 == 3, f"3分量={'OK' if ok3 else pb3} "
                    f"1分量={'拒绝: ' + pb1[0] if pb1 else '竟然通过(!!)'}")

            # ---------------- 子进程部分: 双建字节回归 + θ=0 落盘 + 守卫 ----------------
            print(f"--selftest:   (子进程双建与守卫检查, 每个构建约几秒)")
            cant_tmp = tempfile.mkdtemp(prefix="crect_cant_selftest_")
            print(f"--selftest:   scratch 目录: {cant_tmp}")
            # POSCAR 是**位置参数** (不是 --poscar); 子进程 argv 必须照 CLI 原样拼。
            base_argv = [os.path.abspath(args.poscar), "--element", args.element,
                         "--nonmag", *args.nonmag, "--supercell", *[str(x) for x in sc],
                         "--moment", str(args.moment), "--spin", str(args.spin),
                         "--shell-tol", str(args.shell_tol)]
            d_no = os.path.join(cant_tmp, "no_cant")
            d_ca = os.path.join(cant_tmp, "with_cant")
            d_z0 = os.path.join(cant_tmp, "theta0")
            rc_no, out_no = _run_builder_sub(base_argv + ["--outdir", d_no])
            rc_ca, out_ca = _run_builder_sub(base_argv + ["--outdir", d_ca, "--cant-theta",
                                                          args.cant_theta, "--cant-base",
                                                          args.cant_base, "--cant-azimuth",
                                                          str(args.cant_azimuth)])
            rc_z0, out_z0 = _run_builder_sub(base_argv + ["--outdir", d_z0, "--cant-theta", "0",
                                                          "--cant-base", args.cant_base])
            if rc_no or rc_ca or rc_z0:
                ck9("子进程构建 (no-cant / cant / θ=0) 全部成功", False,
                    f"rc={rc_no},{rc_ca},{rc_z0}\n{out_no[-500:]}\n{out_ca[-500:]}\n{out_z0[-500:]}")

            # (9.7) (a) 不给 --cant-theta 的 JSON 与 cant 版逐字段字节相同 (只多 5 个声明的新键;
            #        noncollinear 是唯一语义上翻转的键)
            if not (rc_no or rc_ca):
                j_no = json.load(open(os.path.join(d_no, "magnetic_orders.json"), encoding="utf-8"))
                j_ca = json.load(open(os.path.join(d_ca, "magnetic_orders.json"), encoding="utf-8"))
                extra = sorted(set(j_ca) - set(j_no))
                miss = sorted(set(j_no) - set(j_ca))
                # 需求 5: 5 个新键在**不给** --cant-theta 时也必须以中性值存在 => 两侧键集完全相同,
                # 不是"cant 版多 5 个键"。(早先写过"多 5 个键", 那是比法错, 不是回归坏。)
                ck9("两侧 JSON 键集完全相同 (5 个新键在无 cant 时也以中性值存在)",
                    set(j_ca) == set(j_no) and all(k in j_no for k in CANT_NEW_KEYS),
                    f"新增={extra} 缺失={miss}")

                # 共线部分的比较必须**分片**做: coefficients_per_cell / _per_mag / order_spec /
                # model_spec 在 cant 打开时被追加/扩展 (需求 3、5 要求的), 整键比必然不同 ——
                # 那是比法错。只把"共线切片"拿出来逐字节比 (json.dumps 保留插入序, 所以键序也比了)。
                def _collo(j):
                    # 5 个新键 + 语义翻转的 noncollinear 由后面的专项检查覆盖, 这里先摘掉,
                    # 只留"旧版就存在、且必须逐位不变"的键。
                    m = {k: v for k, v in j.items()
                         if k not in ("coefficients_per_cell", "coefficients_per_mag",
                                      "order_spec", "model_spec")
                         and k not in CANT_NEW_KEYS and k not in CANT_DIFF_KEYS}
                    m["coefficients_per_cell"] = {k: j["coefficients_per_cell"][k]
                                                  for k in ORDERS}
                    m["coefficients_per_mag"] = {k: j["coefficients_per_mag"][k]
                                                 for k in ORDERS}
                    m["order_spec"] = [e for e in j["order_spec"] if e.get("kind") == "collinear"]
                    ms = dict(j["model_spec"])
                    for k in CANT_NEW_MODEL_KEYS + CANT_DIFF_KEYS:
                        ms.pop(k, None)
                    m["model_spec"] = ms
                    return m
                c_no, c_ca = _collo(j_no), _collo(j_ca)
                diff_keys = [k for k in c_no
                             if json.dumps(c_no[k], ensure_ascii=False)
                             != json.dumps(c_ca[k], ensure_ascii=False)]
                ck9("no-cant 版 JSON 的共线切片与 cant 版**字节**相同 (含键序)",
                    not diff_keys, f"不同的键={diff_keys}")
                ck9("cant 版共线行仍排在最前且顺序不变 (下游按插入序读)",
                    list(j_ca["coefficients_per_cell"])[:len(ORDERS)] == ORDERS
                    and list(j_ca["coefficients_per_cell"])[len(ORDERS):] == cant_names)
                ck9("noncollinear: no-cant=False, cant=True",
                    j_no["noncollinear"] is False and j_ca["noncollinear"] is True)
                ck9("no-cant 版中性键: cant_orders=[] / cant_theta_deg={} / cant_report={} / "
                    "soc_orders=[] / cant_azimuth_deg=0.0",
                    j_no["cant_orders"] == [] and j_no["cant_theta_deg"] == {}
                    and j_no["cant_report"] == {} and j_no["soc_orders"] == []
                    and j_no["cant_azimuth_deg"] == 0.0)
                ck9("cant 版 cant_orders/theta/soc_orders/cant_report 自洽",
                    j_ca["cant_orders"] == cant_names
                    and j_ca["cant_theta_deg"] == {k: float(v) for k, v in cant_theta.items()}
                    and j_ca["soc_orders"] == soc_orders
                    and sorted(j_ca["cant_report"]) == sorted(cant_names)
                    and all(j_ca["cant_report"][k]["base"] == args.cant_base
                            for k in cant_names))
                # 共线 magmom.txt 与 POSCAR: 三个目录 (本进程 outdir / no-cant / cant) 逐字节相同
                same_bytes = []
                for nm in ORDERS:
                    b0 = open(os.path.join(d_no, nm, "magmom.txt"), "rb").read()
                    b1 = open(os.path.join(d_ca, nm, "magmom.txt"), "rb").read()
                    b2 = open(os.path.join(args.outdir, nm, "magmom.txt"), "rb").read()
                    if not (b0 == b1 == b2):
                        same_bytes.append(nm)
                ck9("7 个共线 magmom.txt 在 no-cant / cant / 本次运行三个目录里逐字节相同",
                    not same_bytes,
                    f"不同={same_bytes} 例 sha12(no-cant afm_c)={_sha12(os.path.join(d_no, 'afm_c', 'magmom.txt'))} "
                    f"sha12(cant afm_c)={_sha12(os.path.join(d_ca, 'afm_c', 'magmom.txt'))}")
                posc_diff = [nm for nm in ORDERS
                             if open(os.path.join(d_no, nm, "POSCAR"), "rb").read()
                             != open(os.path.join(d_ca, nm, "POSCAR"), "rb").read()]
                ck9("7 个共线 POSCAR 在 no-cant / cant 两版里逐字节相同", not posc_diff, str(posc_diff))

            # (9.8) (g) cant 序的落盘文件存在, 且 POSCAR 只与基序 POSCAR 差**注释行**
            if not rc_ca:
                miss_files = [f"{nm}/{f}" for nm in cant_names for f in ("POSCAR", "magmom.txt")
                              if not os.path.isfile(os.path.join(d_ca, nm, f))]
                ck9("每个 cant 序的 POSCAR/magmom.txt 都落盘", not miss_files, str(miss_files))
                only_comment = [nm for nm in cant_names
                                if not _poscar_same_except_comment(os.path.join(d_ca, nm, "POSCAR"),
                                                                   os.path.join(d_ca, args.cant_base,
                                                                                "POSCAR"))]
                ck9("cant 序 POSCAR 与基序 POSCAR 只差第 1 行注释 (原子/坐标逐字节相同)",
                    not only_comment, str(only_comment))

            # (9.9) (b 落盘) θ=0: 磁盘上的系数行与共线基序逐位相同; MAGMOM 三分量且 x=y≡0,
            #       z=±moment (state 与共线**完全一致**, 只是按非共线契约写成三列)
            if not rc_z0:
                j_z0 = json.load(open(os.path.join(d_z0, "magnetic_orders.json"), encoding="utf-8"))
                row_z0 = [float(x) for x in j_z0["coefficients_per_cell"]["cant_p0"]]
                ck9("θ=0 的磁盘系数行与共线基序行逐位相同 (cant_p0 == afm_c)",
                    all(row_z0[i] == base_row_[i] for i in range(3))
                    and j_z0["cant_report"]["cant_p0"]["base"] == args.cant_base
                    and j_z0["coefficients_per_cell"]["cant_p0"]
                    == j_z0["coefficients_per_cell"][args.cant_base],
                    f"{row_z0} vs {base_row_}")
                txt0 = open(os.path.join(d_z0, "cant_p0", "magmom.txt"), encoding="utf-8").read()
                v0, nc0 = MV.parse_magmom(txt0, n_atoms=len(atoms_sc))
                m0 = v0.reshape(len(atoms_sc), nc0)
                ck9("θ=0 的 MAGMOM 三分量, x=y 恰为 0, z=σ·moment (与共线态逐位一致)",
                    nc0 == 3 and np.all(m0[:, 0] == 0.0) and np.all(m0[:, 1] == 0.0)
                    and np.all(m0[mag, 2] == args.moment * sig_b),
                    f"ncomp={nc0} x/y 全 0={bool(np.all(m0[:, :2] == 0.0))} "
                    f"z 与 σ·moment 逐位相同={bool(np.all(m0[mag, 2] == args.moment * sig_b))}")

            # (9.10) (f) 全部守卫必须**非零退出**并给出预期中文提示
            g_dir = os.path.join(cant_tmp, "guards")
            guard_cases = [
                ("θ = 90°", ["--cant-theta", "90"], "|θ| < 90"),
                ("θ = -90°", ["--cant-theta", "-90"], "|θ| < 90"),
                ("θ = inf", ["--cant-theta", "inf"], "不是有限值"),
                ("θ = nan", ["--cant-theta", "nan"], "不是有限值"),
                ("θ 非数", ["--cant-theta", "abc"], "不是数"),
                ("θ 空值", ["--cant-theta", ""], "是空值"),
                ("方位角 = 360 (越界)", ["--cant-theta", "8", "--cant-azimuth", "360"], "[0, 360)"),
                ("方位角 = -1 (越界)", ["--cant-theta", "8", "--cant-azimuth", "-1"], "[0, 360)"),
                ("方位角 = nan", ["--cant-theta", "8", "--cant-azimuth", "nan"], "有限"),
                ("cant-base = fm (整体旋转)", ["--cant-theta", "8", "--cant-base", "fm"], "整体"),
                ("cant-base 名字非法", ["--cant-theta", "8", "--cant-base", "afm_z"],
                 "不是允许的 cant 基序"),
                ("两个角度格式化成同一个序名", ["--cant-theta", "8 8.0000000000001"], "序名冲突"),
            ]
            for label, extra, want in guard_cases:
                rc, out = _run_builder_sub(base_argv + ["--outdir", g_dir] + extra)
                ck9(f"守卫非零退出: {label}", rc != 0 and want in out,
                    f"rc={rc} 关键提示 {'命中' if want in out else '缺失: ' + want!r}")
            # 磁性原子数 != 2 (拿本 POSCAR 删到只剩 1 个磁性原子, 写进 scratch)
            lat0, atoms0, coords0, _m0 = parse_poscar(os.path.abspath(args.poscar))
            first_mag = next(k for k, a in enumerate(atoms0) if a == args.element)
            ks = [k for k, a in enumerate(atoms0) if a != args.element or k == first_mag]
            p1 = os.path.join(cant_tmp, "one_mag_POSCAR")
            write_poscar(p1, lat0, [atoms0[k] for k in ks], coords0[ks],
                         comment="selftest: 1 magnetic atom per cell")
            argv1 = [p1] + base_argv[1:]
            rc, out = _run_builder_sub(argv1 + ["--outdir", g_dir, "--cant-theta", "8"])
            ck9("守卫非零退出: 磁性原子数 != 2", rc != 0 and "2 个磁性原子/胞" in out,
                f"rc={rc} 关键提示 {'命中' if '2 个磁性原子/胞' in out else '缺失'}")
            # 亚晶格守卫 (纯函数直调: 单取值 σ / σ 图样与声明不符)
            try:
                cant_sublattice_guard(np.ones(n), args.cant_base, basis, n0, n1)
                msg_u = ""
            except SystemExit as exc:
                msg_u = str(exc)
            ck9("亚晶格守卫: 单取值 σ (= 整体旋转) 被拒", "恰好 2 个自旋亚晶格" in msg_u)
            stripe = sign_pow(n0)
            if np.array_equal(stripe, basis):
                print("--selftest:   [skip] 本晶格的 (-1)^n0 条纹图样恰好等于 basis, 无法用它测错配")
            else:
                try:
                    cant_sublattice_guard(stripe, args.cant_base, basis, n0, n1)
                    msg_w = ""
                except SystemExit as exc:
                    msg_w = str(exc)
                ck9("亚晶格守卫: σ 图样与声明不符 (条纹 σ 冒充 basis) 被拒", "与声明不符" in msg_w)
            # 序名与已有共线序重名 (注入一个假的重名, 用完立刻还原)
            ORDERS.append(cant_order_name(cant_theta[cant_names[0]]))
            try:
                try:
                    cant_order_plan([abs(cant_theta[cant_names[0]])], args.cant_base)
                    msg_d = ""
                except SystemExit as exc:
                    msg_d = str(exc)
            finally:
                ORDERS.pop()
            ck9("cant 序名与已有共线序重名被拒 (非零退出级硬报错)", "重名" in msg_d)

            print(f"--selftest:   cant 检查完成, scratch 目录 {cant_tmp} "
                  f"({'保留(有 FAIL)' if bad else '全部通过 -> 删除'})")
            if not bad:
                shutil.rmtree(cant_tmp, ignore_errors=True)

        print(f"\n--selftest: {'ALL PASS' if not bad else 'FAIL'}")
        for b in bad:
            print(f"  !! {b}")
        if bad:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
