#!/usr/bin/env python
"""test_magnetic_orders.py — 磁序规格层 (P0.1) 的契约测试

覆盖:
  1. 注册表覆盖 9 条既有链路 (8 条 P0–P3 + P4.1 的 kagome), 且每条都带 fixture;
  2. **逐位回归**: 规格层算出的 c_k 必须与 9 条旧链路落盘的 magnetic_orders.json
     完全一致 (容差 1e-9) —— 这是"统一层没改变物理"的唯一证据;
  3. 不变式: 有向计数 ×(-0.5) 与 无向 i<j 两种键枚举恒等 (超胞够大时);
     小超胞下最小镜像会别名坍缩, 正是矩形/蜂窝必须用平移搜索的原因;
  4. 设计矩阵可解性: 条件数有限且不病态; 退化矩阵必须被 check_design 拒绝
     (绝不静默出 J)。

运行:
  python tests/test_magnetic_orders.py
  pytest tests/test_magnetic_orders.py -q
"""
import dataclasses
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
sys.path.insert(0, SCRIPTS)

import magnetic_orders as mo  # noqa: E402

ALL_CHAINS = ("tri_mni2", "rect_crsbr", "sq_k2nif4", "sq_la2cuo4",
              "honeycomb_mnps3", "cluster_nb3i8", "fgt_fe3gete2", "crect_crocl",
              "kagome_fe3sn")

# 每条链路喂进反解的设计矩阵行 (与 build_*.py 的默认选中构型一致; P6 起全共线)
DESIGN_ROWS = {
    # 三角 (P6): 默认不再生成 120°, 由 fm + stripy_a + stripy_ab + uud 给出 rank=4
    # (cond≈3.7e3)。注意 stripy_a 与 stripy_b 的系数列完全相同, 两者同时入选必 rank=3,
    # builder 对重复列有硬闸门。
    "tri_mni2": ("fm", "stripy_a", "stripy_ab", "uud"),
    "rect_crsbr": ("fm", "afm_a", "afm_b", "afm_diag"),
    "sq_k2nif4": ("fm", "stripe_a", "p4_a", "mix_ab"),
    "sq_la2cuo4": ("fm", "stripe_a", "p4_a", "mix_ab"),
    "honeycomb_mnps3": ("fm", "neel", "zigzag", "stripy"),
    "fgt_fe3gete2": ("fm", "afm_ab", "afm_b_vert", "afm_a", "afm_b_plane"),
    # afm_c (胞内两个 basis 反平行) 是唯一翻转 basis 乘积的选定序: C 心键 J1 连的是
    # 同一胞内两个 Cr, 所以 fm/afm_c 只差 c1 (±64)。注意它**不是数学必需**
    # (纯奇偶序组 (fm,afm_a,afm_b,afm_diag) det=-262144 cond=59.17 也能解);
    # 真正的奇异来自"序 + 它的 basis 翻倍版"同时入选 (afm_c_a == afm_a, det 精确为 0)。
    "crect_crocl": ("fm", "afm_c", "afm_a", "afm_b"),
    # kagome (P6): q=0 一族 (fm/q0_ud/q0_duu) 满足 c1 == c2 且 c3 相同 ⇒ 只用 q=0 序时
    # 子块 rank=2 (cond≈2.3e34), J1/J2 不可分离。必须含 √3×√3 调制: sq3_updntri 提供
    # c3 ≠ q0 的信息, sq3_alt 提供 c1 ≠ c2 的信息, 合起来 rank=4 (cond≈1.0e2)。
    # 默认四态全共线; q0_120/sq3_120 只在 --include-noncollinear 时生成。
    "kagome_fe3sn": ("fm", "q0_ud", "sq3_updntri", "sq3_alt"),
}
COND_MAX = 1e4


def test_registry_covers_all_chains():
    assert set(mo.REGISTRY) == set(ALL_CHAINS)
    for key in ALL_CHAINS:
        spec = mo.REGISTRY[key]
        assert spec.fixture, f"{key} 没挂 fixture, 回归会静默跳过"
        assert os.path.exists(mo.resolve(spec.fixture)), f"{key} fixture 不存在: {spec.fixture}"
        assert os.path.exists(mo.resolve(spec.structure.poscar)), f"{key} POSCAR 不存在"


def test_regression_matches_legacy_json():
    bad = {}
    for key in ALL_CHAINS:
        rows, _ctx, verdict = mo.regress(mo.REGISTRY[key], tol=1e-9)
        assert rows, f"{key} 回归没有任何可比对的构型"
        if verdict != "PASS":
            bad[key] = [(n, note) for n, _g, _w, note in rows if note != "ok"]
    assert not bad, f"与旧链路 JSON 不一致: {bad}"


def test_cluster_geometry_is_really_clustered():
    """簇链路必须真的识别出 25 个三聚体, 否则系数会碰巧对上而物理是错的。"""
    _res, ctx = mo.compute(mo.REGISTRY["cluster_nb3i8"])
    assert ctx.n == 75
    assert ctx.n_groups == 25
    assert {len(g) for g in ctx.groups} == {3}


def test_fgt_needs_pair_filter():
    """J3(AA) 与 J4(BB) 距离相同 (3.99 A), 只靠距离无法区分 —— 必须按亚晶格对过滤。"""
    spec = mo.REGISTRY["fgt_fe3gete2"]
    d3 = [s.d for s in spec.shells if s.name == "J3"][0]
    d4 = [s.d for s in spec.shells if s.name == "J4"][0]
    assert d3 == d4
    assert {s.pair for s in spec.shells if s.name in ("J3", "J4")} == {"AA", "BB"}
    # 去掉 pair 过滤后, J3/J4 会因为同距离互相遮蔽 -> fm 的系数必然错
    nopair = tuple(dataclasses.replace(s, pair=None) for s in spec.shells)
    res = mo.compute(dataclasses.replace(spec, shells=nopair))[0]
    want = np.array([-54.0, -9.0, -27.0, -54.0])
    assert not np.allclose(res["fm"].coeffs, want), "pair 过滤没起作用"


def test_directed_and_undirected_are_identical_when_supercell_is_large():
    """docstring 里那条恒等式: -0.5·Σ_有向 ≡ -Σ_{i<j} (超胞够大时不别名)。"""
    spec = mo.REGISTRY["rect_crsbr"]
    und = mo.compute(dataclasses.replace(spec, policy="min_image"))[0]
    shf = mo.compute(dataclasses.replace(spec, policy="shift_search"))[0]
    for name in DESIGN_ROWS["rect_crsbr"]:
        assert np.allclose(und[name].coeffs, shf[name].coeffs, atol=1e-12), name


def test_design_matrices_are_well_conditioned():
    for key, names in DESIGN_ROWS.items():
        res, _ctx = mo.compute(mo.REGISTRY[key])
        coeffs = {n: res[n].coeffs for n in names}
        A = mo.design_matrix(names, coeffs)
        cd = mo.check_design(A, COND_MAX, key)
        assert np.isfinite(cd) and cd <= COND_MAX


def test_check_design_refuses_degenerate_matrix():
    """两行相同的设计矩阵必须被拒绝 —— 不允许静默出一个假的 J。"""
    A = np.array([[1.0, -12.0, -24.0], [1.0, -12.0, -24.0], [1.0, 12.0, -24.0]])
    try:
        mo.check_design(A, COND_MAX, "degenerate")
    except ValueError:
        return
    raise AssertionError("退化设计矩阵没有被 check_design 拒绝")


def test_lattice_presets_are_declared_not_inferred():
    """规格必须显式声明壳距离 (而非运行时自动聚类), 否则无法保证与旧路径同壳。"""
    for key in ALL_CHAINS:
        spec = mo.REGISTRY[key]
        assert spec.shells or spec.cluster, f"{key} 既没有 shells 也没有 cluster 模型"
        for s in spec.shells:
            assert s.d > 0
            assert 0 < s.tol < 1


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    bad = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except AssertionError as exc:
            bad += 1
            print(f"  FAIL  {fn.__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            bad += 1
            print(f"  ERROR {fn.__name__}: {type(exc).__name__}: {exc}")
    print(f"\n{len(fns) - bad}/{len(fns)} 通过")
    return 1 if bad else 0


if __name__ == "__main__":
    print("=" * 72)
    print("磁序规格层 (P0.1) 契约测试")
    print("=" * 72)
    raise SystemExit(_run_all())
