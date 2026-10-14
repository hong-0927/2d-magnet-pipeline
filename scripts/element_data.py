"""元素数据: f 块(稀土/锕系)识别, 供 DFT+U 组装脚本决定 LDAUL。

VASP 的 LDAUL 指定 DFT+U 作用于哪个轨道量子数:
  2 = d 轨道 (3d/4d/5d 过渡金属),  3 = f 轨道 (4f/5f 稀土/锕系)。
流水线原本面向 3d 磁体, 组装脚本硬编码 LDAUL=2。引入 4f 磁体(如 GdI2 的
Gd^{2+} 4f^7)后, U 必须加到 f 轨道(LDAUL=3), 否则会把 U 错误地加到 5d 巡游
电子上。此模块集中维护 f 块元素表, 供各组装脚本按元素自动选择 d/f 轨道。
"""

# 4f 镧系 + 5f 锕系。这些元素的磁性源自 f 电子, U 应加在 f 轨道 (LDAUL=3)。
# 其余磁性元素(3d/4d/5d 过渡金属)默认 U 加在 d 轨道 (LDAUL=2)。
F_BLOCK = {
    # 镧系 (4f)
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    # 锕系 (5f)
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
}


def ldaul_for(element):
    """返回磁性元素 element 的 DFT+U 轨道量子数 LDAUL (字符串): f 块=3, 其余=2。"""
    return "3" if element in F_BLOCK else "2"
