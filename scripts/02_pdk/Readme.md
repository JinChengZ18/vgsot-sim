# 40 nm SOT‑MRAM PDK 行为模型文件结构与信息内容说明

## 一、概述

在先进非易失存储器的工艺设计与电路验证流程中，工艺设计套件（Process Design Kit, PDK）承担着连接器件物理特性与电路级仿真的关键作用。对于 SOT‑MRAM 这类多物理场耦合器件而言，PDK 不仅需要描述静态电学参数，还必须刻画写入动力学、磁化翻转判据以及工艺波动对器件行为的影响。

本文所讨论的 PDK 由两个核心文件构成：`hiksot_40_spectre_v01.scs` 与 `sotmodel.va`。前者是面向 Spectre 仿真环境的模型组织与统计配置文件，后者则是以 Verilog‑A 形式实现的器件级行为模型内核。二者共同构成了完整的 SOT‑MRAM 行为建模体系，其分工与信息层级具有明确的工程逻辑。为保证模型可复现性，本文在原有公式描述基础上进一步补充了各式中参数的具体取值及其在源文件中的定义位置，便于研究人员在模型校验、二次开发或跨平台移植时直接核对。

---

## 二、Spectre 模型封装文件的功能与信息内容

### 2.1 文件定位与作用

`hiksot_40_spectre_v01.scs` 是该 PDK 在 Spectre 仿真环境中的入口文件，其主要职责并非描述器件物理行为本身，而是定义模型的调用方式、工艺角划分以及统计仿真规则。该文件通过 Spectre 的 `section` 机制，将同一器件模型在不同工艺假设下进行组织，使得电路设计人员能够在不修改底层模型代码的情况下，灵活选择仿真条件。

从工程角度看，该文件承担的是模型管理层的角色，而非物理建模层。尤其是在统计建模方面，`.scs` 文件给出了 Monte‑Carlo 与 mismatch 的分布类型与标准差，从而决定了参数扰动的强度与含义；而扰动因子如何作用到器件方程中，则由 Verilog‑A 内核实现。

---

### 2.2 工艺角与统计模型的组织方式

文件中定义了多个 `section`，包括 `sot_tt`、`sot_ttg`、`sot_ttl` 以及 `sot_mc`。这些 section 对应不同的工艺与统计假设，其核心差异体现在对器件参数波动的处理方式上。

在典型工艺角（TT）下，模型通过 `statistics { mismatch { ... } }` 结构引入器件级失配，其中写入通道电阻、MTJ 并联态电阻以及 TMR 的波动均被建模为高斯分布。其数学形式可概括为：

$$
X = X_0 \cdot (1 + \delta_X),
\quad \delta_X \sim \mathcal{N}(0, \sigma_X),
$$

其中 $X$ 表示某一器件参数，$X_0$ 表示其标称值，$\delta_X$ 是零均值随机扰动，$\sigma_X$ 为标准差。需要指出的是，`.scs` 文件中给出的标准差是“扰动因子”的标准差；在该 PDK 的实现中，扰动因子本身以乘法形式作用于标称参数（见 `sotmodel.va` 中 `Rsot = Rsot0*Rsot_mis` 等语句），因此当 $R_{\mathrm{sot\_mis}}$ 服从 $\mathcal{N}(1, 0.07)$ 时，$R_{\mathrm{sot}}$ 的相对标准差为 7%。

在 `sot_tt` 与 `sot_ttg` 两个 section 中，标准差取值分别为：

$$
\sigma(R_{\mathrm{sot\_mis}})=0.07,\quad
\sigma(R_{p\_\mathrm{mis}})=0.07,\quad
\sigma(\mathrm{TMR}_{\mathrm{mis}})=0.04,
$$

在 `sot_ttl` section 中，局部失配假设更弱，其标准差为：

$$
\sigma(R_{\mathrm{sot\_mis}})=0.07,\quad
\sigma(R_{p\_\mathrm{mis}})=0.04,\quad
\sigma(\mathrm{TMR}_{\mathrm{mis}})=0.02.
$$

在 `sot_mc` section 中，统计建模提升到 process 级别，通过 `statistics { process { ... } }` 结构实现 Monte‑Carlo 仿真，其标准差与 `sot_tt` 保持一致：

$$
\sigma(R_{\mathrm{sot\_mis}})=0.07,\quad
\sigma(R_{p\_\mathrm{mis}})=0.07,\quad
\sigma(\mathrm{TMR}_{\mathrm{mis}})=0.04.
$$

这些取值均直接由 `hiksot_40_spectre_v01.scs` 文件中的 `vary ... dist=gauss std=...` 语句给出。

---

### 2.3 子电路封装与参数接口

在 `sot_mc` section 中，文件定义了子电路 `sot_mtj_ckt(p q n)`，该子电路是电路设计人员在原理图或网表中直接实例化的对象。其参数接口包括：

$$
\{R_{\mathrm{sot0}}, R_{p0}, \mathrm{TMR}_0, I_{\mathrm{op0}}, H_x, \mathrm{IniState}\},
$$

并在 `.scs` 文件中给出了默认值：

$$
R_{\mathrm{sot0}}=800~\Omega,\quad
R_{p0}=11~\mathrm{k}\Omega,\quad
\mathrm{TMR}_0=1.2,\quad
I_{\mathrm{op0}}=800~\mu\mathrm{A},\quad
H_x=20,\quad
\mathrm{IniState}=1.
$$

此外，子电路还显式暴露了统计扰动因子参数：

$$
R_{\mathrm{sot\_mis}}=1,\quad
R_{p\_\mathrm{mis}}=1,\quad
\mathrm{TMR}_{\mathrm{mis}}=1,
$$

并在 `statistics { process {...} }` 中定义其抽样分布。在工程使用中，设计者通常不直接手动设置这些扰动因子，而是由 Spectre 的 Monte‑Carlo 引擎根据分布自动生成样本。

---

## 三、Verilog‑A 行为模型文件的功能与信息内容

### 3.1 文件定位与建模层级

`sotmodel.va` 是整个 PDK 的核心文件，其以 Verilog‑A 语言实现了 SOT‑MRAM 单元的完整行为模型。该文件直接描述了器件在电学端口上的电流‑电压关系、磁化状态的演化规则以及写入与读出过程中的物理约束。

与 `.scs` 文件不同，该文件属于物理与行为建模层，其内容决定了仿真结果的物理合理性与工程可用性。更具体地说，`.scs` 文件决定“模型如何被调用以及参数如何被随机化”，而 `.va` 文件决定“给定端口激励与参数样本时，器件如何响应”。

---

### 3.2 端口定义与基本电学关系

模型定义了三个电学端口 $p$、$q$ 与 $n$，分别对应读出端、写入端与公共参考端。写入通道电流由下式给出：

$$
I_{\mathrm{sot}} = \frac{V(q,n)}{R_{\mathrm{sot}}},
$$

其中 $R_{\mathrm{sot}}$ 在 `@(initial_step)` 中被初始化为：

$$
R_{\mathrm{sot}} = R_{\mathrm{sot0}} \cdot R_{\mathrm{sot\_mis}}.
$$

在 `sotmodel.va` 中，$R_{\mathrm{sot0}}$ 与 $R_{\mathrm{sot\_mis}}$ 的默认值分别为：

$$
R_{\mathrm{sot0}} = 800~\Omega,\quad R_{\mathrm{sot\_mis}} = 1.
$$

读出电流在模型中以 `Istt` 表示，并通过一个等效电阻网络计算，其形式为：

$$
I_{\mathrm{stt}} = \frac{2V(p,q) + 2V(p,n)}{4R_{\mathrm{mtj}} + R_{\mathrm{sot}}}.
$$

该表达式体现了三端口结构中写入通道与 MTJ 通道之间的电流分配关系。在模型输出端口电流时，Verilog‑A 采用如下“电流注入”约束：

$$
I(p) \leftarrow I_{\mathrm{stt}},
\quad
I(q) \leftarrow \frac{V(q,n)}{R_{\mathrm{sot}}}-\frac{1}{2}I_{\mathrm{stt}},
\quad
I(n) \leftarrow -\frac{V(q,n)}{R_{\mathrm{sot}}}-\frac{1}{2}I_{\mathrm{stt}}.
$$

这组方程保证了端口电流满足基尔霍夫电流定律 (KCL)，同时将读出支路电流对写入支路进行对称分配，从而形成可收敛的三端口等效网络。

---

### 3.3 MTJ 电阻态与 TMR 建模

模型中 MTJ 的等效电阻由磁化状态决定：

$$
R_{\mathrm{mtj}} =
\begin{cases}
R_p, & \text{P 状态}, \\
R_{ap} = R_p (1 + \mathrm{TMR}), & \text{AP 状态}.
\end{cases}
$$

其中 $R_p$ 在 `@(initial_step)` 中初始化为：

$$
R_p = R_{p0} \cdot R_{p\_\mathrm{mis}},
$$

并且 `sotmodel.va` 中默认参数为：

$$
R_{p0} = 11~\mathrm{k}\Omega,\quad R_{p\_\mathrm{mis}} = 1.
$$

初始时刻模型采用常数 TMR 初始化反并联态电阻：

$$
R_{ap}^{(0)} = R_p\cdot(\mathrm{TMR}_0\cdot \mathrm{TMR}_{\mathrm{mis}} + 1),
$$

其中默认值为：

$$
\mathrm{TMR}_0 = 1.2,\quad \mathrm{TMR}_{\mathrm{mis}} = 1.
$$

需要强调的是，上式仅为初始值；在进入连续仿真后，模型将 TMR 设为 MTJ 电压 $V_{\mathrm{mtj}}$ 的函数。模型先计算 MTJ 电压：

$$
V_{\mathrm{mtj}} = I_{\mathrm{stt}}\cdot R_{\mathrm{mtj}},
$$

随后计算电压依赖 TMR：

$$
\mathrm{TMR}(V_{\mathrm{mtj}}) =
\left(\frac{\mathrm{TMR}_0}{1.2}\right)
\left(
\frac{1}{a V_{\mathrm{mtj}}^2 + b |V_{\mathrm{mtj}}| + c} - 1
\right)
\cdot \mathrm{TMR}_{\mathrm{mis}}.
$$

其中拟合系数在 `@(initial_step)` 中由模型直接给定为：

$$
a = 0.1729,\quad
b = 0.1315,\quad
c = 0.4475.
$$

由此更新反并联态电阻：

$$
R_{ap} = (1+\mathrm{TMR})\cdot R_p.
$$

模型还对工作区间施加了显式约束：当

$$
|V_{\mathrm{mtj}}| > 0.5~\mathrm{V}
$$

时输出 warning，表明该偏压下模型未被定义或未保证有效性。这一约束对于读出电压选择与写入驱动网络的安全性评估具有直接意义。

---

### 3.4 写入动力学与脉宽依赖判据

写入翻转的核心逻辑体现在对写入电流幅值与持续时间的联合判定上。模型首先定义了一个脉宽依赖函数 `Ipw(t_pw,k1,a1,b1,c1)`：

$$
I_{\mathrm{pw}}(t) = \frac{k_1\left(a_1 e^{-b_1 t} + c_1\right)}{700}.
$$

该函数的参数在 `@(initial_step)` 中显式给出：

$$
k_1=0.915,\quad
a_1=1096,\quad
b_1=0.3745,\quad
c_1=596.5.
$$

模型同时定义了用于写入判据的脉宽尺度参数（单位均为 ns，因模型中通过 `$abstime*1E9` 将时间转换为 ns）：

$$
t_{\min} = pw_{\min}=0.2,\quad
pw_{\max}=20,\quad
pw_{\mathrm{std}}=5.
$$

其中 $pw_{\min}$ 表征最小有效脉宽，$pw_{\max}$ 用于定义“进入计时窗口”的最小电流阈值，$pw_{\mathrm{std}}$ 用于 DC 分析下的等效脉宽。

写入阈值电流通过标称操作电流 $I_{\mathrm{op}}$ 与 $I_{\mathrm{pw}}(t)$ 的乘积给出：

$$
I_{\mathrm{th}}(t) = I_{\mathrm{op}}\cdot I_{\mathrm{pw}}(t),
$$

其中 $I_{\mathrm{op}}$ 在 `@(initial_step)` 中被赋值为：

$$
I_{\mathrm{op}} = I_{\mathrm{op0}},\quad I_{\mathrm{op0}}=800~\mu\mathrm{A}.
$$

在瞬态分析中，为了避免对所有微小电流扰动进行计时，模型先以最大脉宽计算一个最小触发电流：

$$
I_{\min} = I_{\mathrm{op}} \cdot I_{\mathrm{pw}}(pw_{\max}).
$$

当写入电流 $I_{\mathrm{sot}}$ 在正确方向上超过 $I_{\min}$ 时，模型触发计时器并开始累计有效脉宽 $pw(t)$。随后模型以实时累计的 $pw(t)$ 计算动态阈值电流 $I_{\mathrm{th}}(pw)$，并在满足

$$
|I_{\mathrm{sot}}| \ge I_{\mathrm{th}}(pw),
\quad pw \ge pw_{\min}
$$

时执行磁化翻转。该实现方式在数值上等价于以经验形式将“翻转难度”表征为脉宽的函数，从而避免显式求解磁化动力学微分方程。

模型在翻转边界附近还引入了一个以 $10~\mu\mathrm{A}$ 为尺度的电阻线性过渡段，用以抑制电阻瞬时跳变引发的收敛困难。以 AP→P 翻转为例，当翻转刚发生且

$$
|I_{\mathrm{sot}}| \le I_{\mathrm{th}}(pw)+10~\mu\mathrm{A}
$$

时，模型将 $R_{\mathrm{mtj}}$ 在 $R_{ap}$ 与 $R_p$ 之间线性插值，其形式在代码中体现为：

$$
R_{\mathrm{mtj}} = R_{ap} - \frac{R_{ap}-R_p}{10~\mu\mathrm{A}}\left(|I_{\mathrm{sot}}|-I_{\mathrm{th}}(pw)\right).
$$

类似的线性插值也用于 P→AP 翻转过程。该处理并不改变翻转判据，而主要服务于数值稳定性，因此在进行 Python 解析化简时可根据仿真目标选择保留或忽略。

---

### 3.5 磁场方向与翻转极性

模型通过外加磁场 $H_x$ 的符号确定写入电流的有效方向。其逻辑关系为：

$$
\mathrm{sign}_{H_x} =
\begin{cases}
+1, & H_x>0,\\
-1, & H_x\le 0,
\end{cases}
$$

并据此定义两类翻转的有效写入极性：

$$
\mathrm{sign}_{\mathrm{AP\rightarrow P}} = -\mathrm{sign}_{H_x},
\quad
\mathrm{sign}_{\mathrm{P\rightarrow AP}} = -\mathrm{sign}_{\mathrm{AP\rightarrow P}}.
$$

在 `sotmodel.va` 中，$H_x$ 的默认值为：

$$
H_x = 20.
$$

因此，在默认配置下 $H_x>0$，从而 $\mathrm{sign}_{\mathrm{AP\rightarrow P}}=-1$，$\mathrm{sign}_{\mathrm{P\rightarrow AP}}=+1$。这意味着当器件初始为 AP 状态时，写入电流需要在模型定义的“负方向”达到阈值才会触发 AP→P 翻转，而 P→AP 则对应相反电流方向。该方向性在 `.scs` 文件层面并未显式呈现，而完全由 `.va` 文件中的符号规则与所给 $H_x$ 决定。







