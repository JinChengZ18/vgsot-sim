# 实验原始数据：物理 ID ↔ 论文标签 映射

## 权威约定（与 06/07/08 拟合脚本一致）

```
physical ID (in filename)   论文/分析标签
─────────────────────────   ──────────────
device 4                    Device A          ← 主基准器件
device 2                    Device B          ← 交叉验证器件
```

数据采集时器件以物理 ID 命名（device 2、device 4）；在论文与分析脚本（`06_psw_t_fitting/fit_from_loops.py`、`07_process_variability/*.py`、`08_sampling_effect/*.py`）中以 Device A、Device B 标签引用。

## 目录结构

本目录下两个子文件夹的命名同时编码物理 ID 与论文标签，便于在任何上下文（数据查找、脚本引用、文档比对）下直接读出对应关系：

```
folder                     物理器件    论文标签
────────────────────────   ─────────   ─────────
raw/device2_paperB/        device 2    Device B
raw/device4_paperA/        device 4    Device A
```

## 数据文件命名约定

```
device{N}pulse width_{T} ns  V_SOT {V} mV Hx 200 Oe Psw= {P}.txt   ← 单点 P_sw 测量
device{N}pulse width_{T} ns 200 Oe.txt                              ← 完整 R-V 滞回回线
```

其中：
- `N` ∈ {2, 4}：物理器件 ID（与所在文件夹的 `device{N}_paper{X}/` 编码一致）
- `T`：脉冲宽度，单位 ns（0.750、1.000、2.000、5.000）
- `V`：写入电压，单位 mV，正负号代表脉冲极性
- `P`：该电压下 100 次重复的成功翻转比例（Wilson 95% 区间在拟合脚本中重新计算）
- `Hx 200 Oe`：固定面内偏置场

## 数据-文件名-代码-结果-文档的一致性

完整调用链以三个等价标识符贯穿：

| 层级 | 标识符 | 示例 |
|---|---|---|
| 数据文件 | `device{N}` 前缀 | `device4pulse width_0.750 ns 200 Oe.txt` |
| 目录 | `device{N}_paper{X}` | `raw/device4_paperA/` |
| 分析脚本 | `Device {X}` | `DEVICE_FILES["A"]` in `fit_from_loops.py` |
| 论文与结果 | `Device {X}` | "Device A、P→AP、$$t_w$$ = 0.75 ns" |

由此 `device 4 ↔ paperA ↔ Device A` 三处标识符在任何文件路径、源码或论文文字中都互为同义，无需额外查表。

## 引用与一致性

需要使用本目录数据的脚本应在头部 docstring 引用本 MAPPING.md，以保证 Device A/B 标签在整条数据管线（06 → 07 → 08）中始终对应 device 4 / device 2。
