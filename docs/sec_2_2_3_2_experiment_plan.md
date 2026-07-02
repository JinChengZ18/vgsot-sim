# §2.2.3.2「保模长的隐式中点法与Cayley变换」实证实验计划与结果

> 由对抗性多智能体审计生成（4 个敌意审稿视角 → 归并 13 条可证伪命题 → 逐条设计实验 → 红队攻击 → 修订），随后逐组实现、运行并经独立验证者复现。
> 本文档 v2：原计划文档在并行会话的 docs 重组中丢失，本版重建并补入实测结果。
> 状态图例：✅ 已完成且独立复现 | 🔄 运行中 | ⏸ 待全量。

---

## 0. 核心症结（贯穿全节的结构性根因）

正文 §2.2.3.2 用**隐式中点法**的方程
`(m_{n+1}-m_n)/Δt = w_n × (m_{n+1}+m_n)/2`
立论，并据此援引隐式中点的全部性质：二阶精度[28]、强1.0/弱2.0收敛[31]、**无需Itô-Stratonovich修正即收敛于正确Boltzmann分布**[13,32]。

但已发布的 `switching_vector()`（[dynamic_switching_vector.py:200-207](../src/vgsot_sim/dynamic_switching_vector.py)）只在**当前态 m_n** 处对 ω 求值一次，再施加单步闭式 Cayley 旋转——**没有任何 fixed-point / 中点迭代逻辑**。即正文第521行「w_n 用当前步状态 m_n 评估」与同段的中点方程及第557行的收敛定理**自相矛盾**：运行的是「显式-ω 几何 Cayley 步」，而非隐式中点法。

**实验证实的后果**：Cayley 的*几何*性质（|m|=1 无条件保持，SC1）为真且与求值点无关（实验 C ✅）；但*精度*性质已被证伪——发布算法全局**一阶**（p=1.013），而用同一批原语搭出的真·中点恢复**二阶**（p=2.001）（实验 A ✅）。统计性质（SC5 Boltzmann）严谨 pilot 尚不能在 3σ 水平裁决，需 60 ns 全量（实验 B 🔄）。

**由此引出的关键抉择（等实验数据齐后定，用户已确认）**：
- **路线①（改文本）**：如实把该节描述为「显式-ω 几何积分器」——保模长稳定、全局一阶；删除/弱化二阶、强弱阶、Boltzmann-无修正等论断。
- **路线②（改代码）**：把内核升级为**真·中点**（对 ω 做 2–3 次不动点迭代到 (m_n+m_{n+1})/2，实验 A 已给出参考实现），恢复二阶（已实证）与可能的 Boltzmann 性质，再重跑全部图表。

---

## 1. 命题清单与最终判定

| ID | 主题 | 优先级 | 嫌疑判定 | 实验 | 实测判定 |
|----|------|:---:|------|:---:|------|
| C1 | 显式-ω vs 隐式中点（根因） | 5 | 被实现证伪 | A | ✅ **证实**：p_pub=1.013 vs p_mid=2.001 |
| C3 | SC2 确定性二阶 O(Δt²) | 5 | 被实现证伪 | A | ✅ **证伪**（对发布代码），高置信 |
| C8 | 全节零数值物证 | 4 | likely_true | A/B | ✅ 已补：`tests/test_integrator_order.py` 5 用例通过 |
| C2 | SC5 Boltzmann 无需修正 | 5 | 被实现证伪(探针) | B | ✅（初步，cayley 段 60ns 已收）**不被证伪**：T_eff/T = 0.978/0.987/1.003/1.007（dt=0.05–0.4ps，σ≈0.015–0.021，~90τ_int 充分平衡）；Richardson dt→0 截距 0.979±0.014（-1.5σ）与 0.990±0.015（-0.7σ，两估计器），dt 斜率 +1.2σ——均距 3σ 判伪线甚远 → 发布 Cayley 稳态在 1–2% 内即 Boltzmann。探针 1.8–2.0× 与 pilot 0.83 均为伪影（前者估计器偏置、后者欠平衡）。euler/midpoint 对照段补跑中（首跑被 11:22 系统重启杀死于 euler 段，已加逐 cell 检查点防复发） |
| C10 | 归一化注入稳态偏置 | 2 | depends | B | 🔄 随 B 全量一并裁决 |
| C9 | SC1 无条件保模长 | 2 | likely_true | C | ✅ **成立**：物理网格 ‖m‖-1≤2.2e-16；裸轨迹 0.75ns max 1.9e-15；但「无条件/任意Δt」须软化（广告网格 \|ω\|→1e14 处矩阵式 12 处奇异，闭式仍精确） |
| C6 | SC3 强收敛阶 1.0 | 3 | depends(探针~0.5) | E | ✅ **成立**：headline p=0.972 CI[0.945,1.000] R²=0.998（Ms=0.003 档，噪声主导比 2112）；三档 1.020/0.972/1.096；GBM 门控 EM=0.5✓/Milstein=0.992✓；自参考对照 p=1.211 演示独立参考规避的膨胀。保留：最强噪声档 p_finest3=0.858（轻微细端退化，值得一句限定） |
| C7 | SC4 弱收敛阶 2.0 | 3 | depends | E | ✅ **成立**：p=1.97–2.04 双参考一致 R²≈1.000（轴对称点）；破缺对称探针（H_x=6.7% 阱强）池化 p=1.825 存活——非对称性伪影（细端受参考差异地板 3e-4 限制、粗端大转角饱和，干净窗口斜率 2.2–3.6）。**体制限定**：此为噪声主导体制；生产 P_sw 在脉冲期驱动主导，实验 A 的确定性一阶误差才是大步长统计精度的控制项 |
| C5 | 大步长宽容度有支撑 | 4 | likely_false | B/E | 🔄 随 B/E 裁决 |
| C11 | SC6 球坐标 1/sinθ 奇点 | 2 | likely_true | D | ✅ **成立**：guard-OFF 斜率 -1.0004；但 (a) 极点保护实际 θ≤1e-9 才触发（非 1e-8）；(b)「无极点」须改为「内核场组装规避 1/sinθ 除法（笛卡尔 dm/dt 处处有限，span 1.087）」 |
| C12 | SC7 显式 Euler 模长爆炸 | 2 | likely_true | D | ✅ **成立**：\|m\|-1 精确循解析律 √(1+(Δt\|ω⊥\|)²)-1（dt=1e-12→+0.41421356=√2-1；1e-11→+9.0499=√101-1）；但真实 \|ω\|≈3.0e10 下发散阈在 dt·\|ω\|≳1 即 ~3.3e-11，增长律应以 \|ω⊥\| 表述 |
| C4 | 图表用哪个积分器 | 4 | 被实现证伪 | F | ✅ **证实**：fig2.10 全部 cayley 调用（0 次 euler）、`run_piecewise_direct_excitation` 默认='cayley' → 第743行「默认 Euler-球坐标」陈旧、第757行为真 |
| C13 | θ_SH=0.04 还是 0.066 | 2 | 被实现证伪 | F | ✅ **证实**：运行时 θ_SH=0.066，出图脚本无覆盖 → 图2.10标题/第743/76-77行溯源陈旧。200 种子系综（NEGATIVE I_SOT，P→AP）：θ_SH=0.066 → I_50=1289.9 µA CI[1276,1305]（距实验 1152.3 µA 差 12%，为较优匹配）；θ_SH=0.04 → I_50=2076.5 µA（偏 80%，排除）。**新发现**：(a) §2.2.4.6「匹配至 1%」在当前代码路径下不复现（对应仓库挂起项 #12 重标定，粗估 θ_SH≈0.074 可回到 1152 µA，须用与原标定一致的判据重扫）；(b) 高过驱动 P_sw 饱和于 ~0.78 而非 1（1400–2000 µA 平台，回跳/过冲迹象） |

**B pilot 的意外新发现**：Euler-球坐标积分器的平衡有效温度 **T_eff/T≈0.44（冷 6.3σ，dt-稳健）**——比 Cayley 的问题更严重且方向相反。鉴于 fig2.11 SER MC 默认走 euler 路径（F 证实），此发现直接关系热统计图表的可信度，须在 B 全量中确认并写入正文修订。

---

## 2. 六组实验设计与产物

（设计细节、红队修订与判据同 v1；此处保留执行要点与产物索引。）

### A — 确定性收敛阶 ✅
脚本 `scripts/02_integrator/order_of_accuracy.py` + `tests/test_integrator_order.py`（5 用例通过）。
NON=0，I_SOT=-1.5mA，T=0.5ns，dt∈{8..0.125}ps，参考=中点@T/131072。
**结果**：p_pub=1.013(R²=1.000)、p_mid=2.001(R²=1.000)、p_euler=1.019、常ω对照 p=2.0000（对解析 Rodrigues）→ 一阶损失源于 ω(m(t)) 的**时间变化**（左端点冻结），与阻尼无关。验证者用独立 scipy DOP853 参考复现全部斜率至 4 位有效数字；有/无归一化逐位相同（排除归一化掩蔽）；轨迹 min sinθ=0.296（排除极点混杂）。
产物：`result/sec_2_2_3_2/A/order_of_accuracy.{png,json}`。

### B — Boltzmann 有效温度 🔄（全量运行中）
脚本 `scripts/11_boltzmann_teff/teff_dtsweep.py`。零驱动纯单轴势阱（h_ex 置零 + ENE=0，Δ=48.515 解析）；T_eff 用**精确有界基矩反演 + WLS 直方图斜率**（谐估计器构造偏置 0.989 仅列表不使用）；固定物理采样时长；τ_int + 长度平台稳态证明；dt→0 Richardson + block-bootstrap。
**Pilot（3ns×12轨，独立验证者逐位复现）**：cayley dt→0 截距 0.829±0.133（-1.3σ，不显著）；「单调过热趋势」被验证者证明是种子涨落；τ_int~180-330ps → 3ns 仅 ~17τ_int 欠平衡，**诚实判定 inconclusive**。euler_spherical 锚点 0.444±0.089（**冷 6.3σ**，dt-稳健——新发现）。
**全量**：60ns × 4dt{0.05,0.1,0.2,0.4}ps × 48轨 × 3积分器，n_boot=2000，~10.5h，已后台启动（`--full`，log: `result/sec_2_2_3_2/B/full_run.log`）。

### C — 保模长与归一化装饰性 ✅
脚本 `scripts/02_integrator/test_norm_preservation.py`。
**结果**：物理网格(63格×8方向) max‖m‖-1=2.22e-16；广告网格闭式仍 8.88e-16 但矩阵式在 \|ω\|∈{1e11..1e14} 12 处奇异（闭式 vs 矩阵最大失配 1.08e3）→「无条件/任意Δt」措辞须限定为闭式实现；**载荷测试**：裸笛卡尔 0.75ns（无归一化、无球坐标往返、m_z 从 +1 摆到 ~0.02 过赤道）max‖m‖-1=1.89e-15 ≪ 1e-12 → 归一化纯装饰。Part C（管线 ON/OFF 1.02e-12）仅为鲁棒性对照（驱动层球坐标重建本身重置模长）。
产物：`result/sec_2_2_3_2/C/`。

### D — Euler 爆炸 + 球坐标极点 ✅（部分成立，纠两处夸大）
脚本 `scripts/09_simulation_figures/c12_norm_stability_sweep.py`、`scripts/11_pole_singularity/probe_pole_singularity.py`。
**结果**：见上表 C11/C12 行。验证者逐位复现（确定性无 RNG）；「发散onset 5e-11」为扫描离散化措辞，真实过界 3.29e-11。
产物：`result/sec_2_2_3_2/D/`。

### E — 随机强/弱收敛阶 ✅（两项均"成立"，与探针预期相反）
脚本 `scripts/11_strong_order/strong_order_brownian.py`（标量，GATE-0 GBM 门控 + GATE-1 确定性锚 + GATE-2 噪声主导比≥50 + 独立中点参考 + 自参考退化对照）、`scripts/10_weak_order/weak_order_audit.py`（按轨迹向量化，M=150k，双参考交叉验证，预平衡短 horizon）。
首轮工作流中执行智能体因长段无输出被停滞判死（6 次）；脚本完整，直接后台重跑完成。向量化实现与标量公共 API 的一致性已抽查（ω⊥ 相对差 5.3e-15，单步差 2.5e-16）。
**结果**：见上表 C6/C7 行。净结论——援引的随机阶（强 1.0/弱 2.0）对发布的显式-ω 格式**经验成立**（实测体制内）；被证伪的是确定性二阶（实验 A）与「隐式中点」标签本身（C1）。修订方向：改正格式命名与确定性阶数，随机阶表述的依据从「援引隐式中点定理」改为「本工作实测验证」（可作方法学小贡献点）。
产物：`result/sec_2_2_3_2/E/strong_order.{png,json}`、`weak_order.{png,json}`、`weak_order_symbreak.{log,json}`。

### F — 图表与参数溯源审计 ✅（含 200 种子系综，全部完成）
脚本 `scripts/09_simulation_figures/integrator_audit.py`、`audit_theta_sh.py` + `tests/test_figure_provenance.py`。
**结果**：见上表 C4/C13 行（判定"claim_false_confirmed"高置信）。系综复核确认 0.066 为运行实况且是较优匹配，但距实验阈值仍有 12% 系统性偏差——§2.2.4.6「1% 匹配」声称属 FL-SOT 修正前旧代码路径，须随 #12 重标定后改写；高过驱动 P_sw 平台 ~0.78 可作图 2.11 讨论素材。
产物：`result/sec_2_2_3_2/F/`（`theta_sh_provenance.json`、`theta_sh_ensemble.png`、`integrator_provenance.json`、`full_ensemble.log`）。

---

## 3. 对正文的修订清单（待 B/E 收口后定稿）

1. **第521行 + 方程矛盾**：按抉择路线①或②统一「隐式中点」表述与实现。
2. **SC2（第557行二阶）**：路线① → 改述全局一阶并引实验 A 图；路线② → 升级内核后保留并引 A 的中点曲线。
3. **SC5（Boltzmann 无修正）**：以 B 全量结果为准改写；无论方向，需补充 Euler-球坐标 0.44× 过冷的发现及其对 fig2.11 的影响评估。
4. **SC1（第555行）**：保留结论，软化「排除了手动归一化的需要」（代码每步仍归一化，属防御性 O(eps) 清理）与「无论 Δt 取何值」（限定为闭式实现；矩阵式在极端 \|ω\|Δt 下奇异）。
5. **SC6「无极点」**：改为「内核场组装规避 1/sinθ 除法（笛卡尔 dm/dt 处处有限）」；补注极点保护实际触发边界 θ≤1e-9。
6. **SC7 稳定性阈**：增长律以 \|ω⊥\| 表述；真实工况发散阈 dt·\|ω\|≳1（~3.3e-11 s）。
7. **第743行**：「复现章节图表默认 Euler-球坐标」→ 改为实况（fig2.10=cayley、fig2.11 SER=euler 且两者混用需说明）；θ_SH 溯源统一为 0.066（FL-SOT 修正后 Cayley 标定）。
7b. **§2.2.4.6 第751行「匹配至 1% 以内」**：当前代码路径下系综 I_50=1290 µA vs 实验 1152 µA（12%），该声称属修正前旧路径——随 #12 重标定 θ_SH（粗估 ~0.074，判据须与原标定一致）后更新数字，或如实改述为「校准至 ~12% 并列明重标定为待办」。
8. **SC3/SC4（强/弱阶）**：按 E 结果改写；若强阶 ~0.5、弱阶 ~1 证实，「大步长宽容度核心优势」须重新论证或删除。
9. **试错-修正记录素材**（按论文规范第10条）：探针 1.8–2.0× 过热被严谨估计器推翻、Euler 意外过冷、执行智能体停滞与重启等，可作脚注级真实迭代素材。

## 4. 复现命令

```bash
cd vgsot-sim
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/02_integrator/order_of_accuracy.py      # A
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/02_integrator/test_norm_preservation.py # C
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/09_simulation_figures/c12_norm_stability_sweep.py  # D1
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/11_pole_singularity/probe_pole_singularity.py      # D2
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/09_simulation_figures/integrator_audit.py          # F-A
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/09_simulation_figures/audit_theta_sh.py            # F-B (加 --full --n-seeds 200 出系综)
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/11_boltzmann_teff/teff_dtsweep.py                  # B pilot (加 --full 过夜)
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/11_strong_order/strong_order_brownian.py           # E-strong
PYTHONIOENCODING=utf-8 PYTHONPATH=src python scripts/10_weak_order/weak_order_audit.py                  # E-weak
PYTHONPATH=src python -m pytest tests/test_integrator_order.py tests/test_figure_provenance.py -q       # 回归
```
