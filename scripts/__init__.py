"""scripts 包: 2d_magnet 流水线的可执行脚本集合。

本目录既是「按脚本文件路径调用」(driver 用 subprocess 跑 `scripts/<name>.py`) 的脚本库，
也是「可安装的 Python 包」(非可编辑 `pip install .` 会把整个目录连 `literature_J.json`
一起装进 site-packages)。run_2d_magnet.py 用 `import scripts` + `scripts.__file__` 定位
本目录, 因此三种运行方式都能找到脚本: 源码树直跑 / 可编辑安装 / 非可编辑安装。
"""
