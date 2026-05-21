# 2.1三端SOT-sMTJ器件模型

sMTJ的器件层模型需要在静态电学特性、自由层磁化动力学、写入概率以及行为级抽象四个层次之间保持自洽，才能够支撑后续阵列级与系统级的建模工作。本节首先建立三端SOT-sMTJ的T型等效电路并阐明读写路径解耦的物理意义，随后由Néel–Brown热激活模型给出受驱动调制的翻转概率，综合SOT与VCMA两类驱动机制得到联合写入概率表达式，最终将复杂的指数嵌套形式压缩为工作区内便于工程使用的Sigmoid近似。[1] [2]

## 2.1.1三端SOT-sMTJ器件结构与基础电学模型

三端SOT-sMTJ由垂直磁各向异性（PMA）的MTJ堆叠（参考层/隧穿势垒/自由层）置于重金属（HM）沟道之上构成。MTJ顶端为读出端$$T_1$$，沟道两端分别为$$T_2$$与$$T_3$$。该拓扑使读出与写入路径在物理上分离：读出操作仅在$$T_1$$与沟道接地端之间施加低偏压，依靠TMR效应产生与磁化态相关的隧穿电流；写入操作则将大电流沿$$T_2$$至$$T_3$$流过HM沟道，由自旋霍尔效应产生的横向自旋积累驱动自由层翻转。读写解耦避免了STT-MTJ读写共用同一隧穿路径所引起的可靠性退化问题[5]，是SOT-sMTJ作为可编程随机源的结构基础。

将MTJ等效为电阻$$R_{MTJ}$$、HM沟道沿长度方向均分为两段电阻$$R_{SOT}/2$$，三端器件可建模为T型电阻网络，其中沟道电阻由几何与材料参数给定

$$
R_{SOT} = \frac{\rho_{SOT} L_{SOT}}{W_{SOT} T_{SOT}},
$$

$$\rho_{SOT}$$、$$L_{SOT}$$、$$W_{SOT}$$与$$T_{SOT}$$分别为HM沟道的电阻率、长度、宽度与厚度。

![三端SOT-sMTJ器件结构与T型等效电路网络](figures/fig_01_t_circuit.png)

**图2.1：三端SOT-sMTJ器件结构与T型等效电路网络。** (a)器件物理堆叠结构，包括位于顶端的MTJ（参考层/MgO势垒/自由层）与底部的HM沟道；MTJ顶电极对应端口$$T_1$$，HM沟道两端分别引出端口$$T_2$$与$$T_3$$。(b)T型电阻网络，$$R_{MTJ}$$表示自由层与参考层之间的隧穿电阻，沟道被等分为两段$$R_{SOT}/2$$分别接$$T_2$$与$$T_3$$，三段电阻在中间节点N处汇合。读出操作通过$$T_1$$施加小偏压获取$$R_{MTJ}$$状态，写入操作通过$$T_2$$、$$T_3$$沿沟道施加大电流产生SOT，两条电流路径在物理上解耦。

设端口$$T_1$$、$$T_2$$、$$T_3$$电压分别为$$V_1$$、$$V_2$$、$$V_3$$，中间节点电压为$$V_N$$。节点N处的KCL给出

$$
\frac{V_1-V_N}{R_{MTJ}} + \frac{V_2-V_N}{R_{SOT}/2} + \frac{V_3-V_N}{R_{SOT}/2} = 0,
$$

解得$$V_N = [R_{SOT}V_1 + 2R_{MTJ}(V_2+V_3)]/(R_{SOT}+4R_{MTJ})$$。代入欧姆关系即可得到流过MTJ的电流

$$
I_{MTJ} = \frac{V_1-V_N}{R_{MTJ}} = \frac{4V_1 - 2V_2 - 2V_3}{4R_{MTJ}+R_{SOT}},
$$

以及流入端口$$T_2$$、$$T_3$$的电流

$$
I_2 = \frac{V_2-V_3}{R_{SOT}} - \frac{1}{2}I_{MTJ}, \qquad I_3 = -\frac{V_2-V_3}{R_{SOT}} - \frac{1}{2}I_{MTJ}.
$$

后两式中的$$\pm\tfrac{1}{2}I_{MTJ}$$项分别对应$$I_{MTJ}$$经由沟道向$$T_2$$与$$T_3$$对称分流的部分。这两段分流在沟道左右两侧方向相反，所产生的横向自旋积累在自由层下方相互抵消，对SHE-SOT效应没有净贡献。驱动自由层翻转的有效SOT电流由$$T_2$$至$$T_3$$的贯穿分量给出

$$
I_{SOT} = \frac{V_2-V_3}{R_{SOT}}.
$$

该表达式从电路层面验证了三端结构的读写解耦：施加于$$T_1$$的MTJ偏压$$V_{MTJ} = V_1 - V_N$$即使较大也不会改变$$I_{SOT}$$的方向与幅值，从而允许将$$V_{MTJ}$$独立用作VCMA能垒调制变量（详见2.1.3节）。

MTJ的隧穿电阻由Brinkman模型描述，该模型基于WKB近似处理梯形势垒的电压依赖隧穿[3]。在低偏压下，平行态电阻可表为势垒厚度与势垒高度的函数

$$
R_P = \frac{t_{ox}}{F \cdot A_{MTJ} \cdot \phi_{ox}}\exp\!\left(\frac{2 t_{ox}}{\hbar}\sqrt{2 m_e e\, \phi_{ox}}\right),
$$

其中$$t_{ox}$$与$$\phi_{ox}$$分别为隧穿势垒层厚度与势垒高度，$$F$$为根据电阻—面积乘积$$RA$$标定的拟合因子，$$A_{MTJ}$$为MTJ截面积，$$e$$为元电荷，$$m_e$$为电子质量，$$\hbar$$为约化普朗克常数。指数因子主导了$$R_P$$对势垒厚度的强敏感性，$$t_{ox}$$每变化$$0.1\,\mathrm{nm}$$量级即可使$$R_P$$变化数倍，构成MTJ阻值工艺难以精细控制的根本来源之一。

反平行态电阻$$R_{AP}$$由TMR比定义。在Jullière二流模型下，理论TMR可由两侧铁磁电极的有效自旋极化率$$P$$给出

$$
\mathrm{TMR}_0 = \frac{R_{AP}-R_P}{R_P} = \frac{2P^2}{1-P^2}.
$$

实际TMR随MTJ两端电压$$V_{MTJ}$$增大而衰减，主要源自界面处的非弹性磁子辅助隧穿过程，即磁子吸收/发射通道随偏压打开后破坏了自旋守恒[4]。Moodera等的实验观测显示，TMR随$$V_{MTJ}$$近似呈Lorentzian形衰减，常用的工程参数化形式为

$$
\mathrm{TMR}(V_{MTJ}) = \frac{\mathrm{TMR}_0}{1 + V_{MTJ}^2/V_h^2},
$$

其中特征电压$$V_h$$定义为TMR下降至理论值一半时所对应的偏压，CoFeB/MgO/CoFeB体系中典型值约为$$0.3\text{--}0.7\,\mathrm{V}$$。当自由层与参考层磁化方向夹角为$$\theta$$时，MTJ电阻满足如下含偏压依赖的角度公式

$$
R_{MTJ}(V_{MTJ},\theta) = R_P\,\frac{1+V_{MTJ}^2/V_h^2 + \mathrm{TMR}_0}{1+V_{MTJ}^2/V_h^2 + 0.5\,\mathrm{TMR}_0(1+\cos\theta)}.
$$

在$$\theta = 0$$（平行）时分子分母相同，$$R_{MTJ} = R_P$$；在$$\theta = \pi$$（反平行）且$$V_{MTJ} \ll V_h$$时，$$R_{MTJ} \to R_P(1+\mathrm{TMR}_0) = R_{AP}$$。该公式同时刻画了TMR随偏压的压缩与随磁化夹角的连续变化，是磁动力学求解过程中读出电流计算的基础。

读路径电压$$V_{MTJ}$$与写路径电流$$I_{SOT}$$在端口拓扑层面已实现解耦，但二者在磁化层内部仍通过共享自由层的能垒调制相互耦合：$$V_{MTJ}$$经VCMA压低能垒，$$I_{SOT}$$经SHE提供翻转驱动力，两者协同决定动态翻转概率。这一耦合关系将在下两节中具体展开。

## 2.1.2自由层磁化动力学与热激活翻转概率

为给出MTJ在不同驱动条件下的翻转概率，本节采用宏自旋近似下的紧凑模型展开推导，主要包含Néel–Brown热激活模型与受驱动调制的能垒缩减表达式。

对具有PMA的MTJ自由层，热稳定因子定义为$$\Delta = E_b/(k_B T)$$，度量两稳定磁化态之间的能垒高度，其中$$E_b = K_{\mathrm{eff}}V \approx \tfrac{1}{2}\mu_0 H_k M_s V$$，$$V$$为自由层体积，$$K_{\mathrm{eff}}$$为有效各向异性能密度，$$H_k$$为等效各向异性场，$$M_s$$为饱和磁化强度。在亚临界写入条件下，sMTJ的翻转受热涨落主导，由Néel–Brown热激活模型描述[1,2]：磁化在能垒约束下发生热激活逃逸，平均停留时间满足Arrhenius关系$$\tau = \tau_0 \exp(\Delta)$$，其中$$\tau_0$$为亚纳秒至纳秒量级的尝试时间，典型取值约为$$1\,\mathrm{ns}$$。在持续时间为$$t_w$$的写脉冲下以泊松过程近似，翻转概率为

$$
P_{\mathrm{sw}}(t_w) = 1 - \exp\!\left(-\frac{t_w}{\tau_0 e^{\Delta}}\right).
$$

外场、自旋转移矩或SOT驱动存在时，常将其等效为能垒压低$$\Delta_{\mathrm{eff}} = \Delta(1-I/I_{c0})^n$$（$$I < I_{c0}$$），其中$$I_{c0}$$为零温临界驱动，指数$$n$$由器件物理工作区间决定：亚临界热激活区（驱动远低于$$I_{c0}$$、热涨落主导翻转）能垒随驱动呈线性压低，对应$$n = 1$$；接近零温弹道极限（驱动接近$$I_{c0}$$、热涨落可忽略）能垒压低呈抛物线行为，对应$$n = 2$$。概率计算工作区位于亚临界热激活区，故下文取$$n = 1$$。综合得带驱动的翻转概率

$$
P_{\mathrm{sw}}(t_w, I) = 1 - \exp\!\left[-\frac{t_w}{\tau_0}\exp\!\left(-\Delta\!\left(1-\frac{I}{I_{c0}}\right)^n\right)\right].
$$

定义二值随机变量$$m \in \{0,1\}$$表示器件翻转事件，则$$m \sim \operatorname{Bernoulli}(P_{\mathrm{sw}})$$，sMTJ自然成为以驱动量为参数的可编程伯努利源。该模型在概率计算工作区与器件实测的吻合度将在2.3节通过专门的实验对比给出，本节不再赘述。

## 2.1.3 SOT与VCMA联合驱动模型

SOT-MRAM的写入电流流经HM沟道，由自旋霍尔效应或Rashba–Edelstein效应产生横向自旋积累，对自由层施加阻尼型等效场

$$
H_{\mathrm{DL}} = \frac{\hbar \theta_{\mathrm{SH}} I_{SOT}}{2 e \mu_0 M_s t_f A_c},
$$

其中$$\theta_{\mathrm{SH}}$$为有效自旋霍尔角，$$I_{SOT}$$即2.1.1节由$$(V_2-V_3)/R_{SOT}$$给出的沟道贯穿电流，$$t_f$$为自由层厚度，$$A_c$$为沟道截面积。零温近似下临界SOT电流标度为

$$
I_{c0}^{\mathrm{SOT}} \propto \frac{2e}{\hbar}\,\frac{\alpha M_s t_f A_c}{\theta_{\mathrm{SH}}}\,H_k,
$$

降低$$M_s$$、$$t_f$$、阻尼系数$$\alpha$$或提高$$\theta_{\mathrm{SH}}$$均可减小写入电流[6]。

VCMA机制通过MTJ两端的偏压$$V_{MTJ}$$改变自由层/氧化层界面的各向异性，界面各向异性能面密度$$K_i$$随电场线性变化[7]

$$
K_i(V_{MTJ}) = K_i(0) - \xi\,\frac{V_{MTJ}}{t_{ox}},
$$

其中$$\xi$$为VCMA系数（量纲$$\mathrm{fJ/(V\cdot m)}$$），$$t_{ox}$$已由2.1.1节给出。此处$$V_{MTJ}$$正是2.1.1节T型等效电路中$$T_1$$与中点N之间的电压降，电学模型与磁动力学模型由此通过共享变量自然衔接。将界面项折算到体各向异性后，能垒变为

$$
E_b(V_{MTJ}) = \left[K_{\mathrm{eff}}(0) - \frac{\xi}{t_f}\,\frac{V_{MTJ}}{t_{ox}}\right]V,
$$

代入Néel–Brown模型可得VCMA辅助下的翻转概率

$$
P_{\mathrm{sw}}(t_w, V_{MTJ}) = 1 - \exp\!\left(-\frac{t_w}{\tau_0}\exp\!\left[-\Delta_0 + \beta_{\mathrm{VCMA}} V_{MTJ}\right]\right),
$$

其中$$\Delta_0 = E_b(0)/(k_B T)$$，$$\beta_{\mathrm{VCMA}} = \xi V/(k_B T \cdot t_f t_{ox})$$。VCMA将翻转概率曲线沿$$V_{MTJ}$$方向左移，使较低写入能量下即可获得可控的随机翻转。

VCMA与SOT联合作用时，自然的建模顺序是VCMA先降低各向异性与能垒（$$H_k \to H_k(V_{MTJ})$$），SOT再提供翻转驱动力，由此得到偏压依赖的SOT临界电流$$I_{c0}^{\mathrm{VCMA}}(V_{MTJ})$$及联合翻转概率

$$
P_{\mathrm{sw}}(t_w, I_{SOT}, V_{MTJ}) = 1 - \exp\!\left[-\frac{t_w}{\tau_0}\exp\!\left(-(\Delta_0 - \beta_{\mathrm{VCMA}}V_{MTJ})\!\left(1-\frac{I_{SOT}}{I_{c0}^{\mathrm{VCMA}}(V_{MTJ})}\right)^n\right)\right].
$$

由于$$H_k(V_{MTJ})$$随VCMA偏压下降，$$I_{c0}^{\mathrm{VCMA}}$$亦相应降低，VCMA辅助SOT得以减小写电流与写能耗[8]。该联合模型给出了三端sMTJ作为可编程伯努利源的完整描述：外部电路只需指定脉宽$$t_w$$、SOT电流$$I_{SOT}$$与MTJ偏压$$V_{MTJ}$$这一三元组，器件即输出对应概率的随机比特。

## 2.1.4工作区Sigmoid近似

联合翻转概率$$P_{\mathrm{sw}}$$是驱动量的指数嵌套函数，严格数学形式属于Gumbel型累积分布而非标准Sigmoid（Logistic）函数。但在临界过渡区（$$P_{\mathrm{sw}} \approx 0.05$$至$$0.95$$）高阶项影响甚小，可作Sigmoid近似。设等效驱动变量$$u$$综合脉宽、SOT电流与VCMA偏压的影响，将有效能垒在工作点$$u_0$$附近一阶展开$$\Delta E_{\mathrm{eff}}(u) \approx \Delta E_{\mathrm{eff}}(u_0) - \kappa(u-u_0)$$，其中$$\kappa = -\partial \Delta E_{\mathrm{eff}}/\partial u|_{u_0}$$为能垒压低效率，则

$$
P_{\mathrm{sw}}(u) = 1 - \exp[-A\exp(B(u-u_0))], \qquad B = \frac{\kappa}{k_B T}.
$$

定义中值翻转点$$u_{\mathrm{th}}$$满足$$P_{\mathrm{sw}}(u_{\mathrm{th}}) = 0.5$$，将$$P_{\mathrm{sw}}$$转换到logit域$$L(u) = \ln[P_{\mathrm{sw}}/(1-P_{\mathrm{sw}})]$$并在$$u_{\mathrm{th}}$$处一阶泰勒展开（零次项$$L(u_{\mathrm{th}}) = 0$$），得$$L(u) \approx 2B\ln 2 \cdot (u-u_{\mathrm{th}})$$。记$$\beta_s = 2\kappa\ln 2/(k_B T)$$，还原概率形式即得过渡区Sigmoid近似

$$
P_{\mathrm{sw}}(u) \approx \sigma(\beta_s(u-u_{\mathrm{th}})) = \frac{1}{1 + e^{-\beta_s(u-u_{\mathrm{th}})}}.
$$

参数$$u_{\mathrm{th}}$$为等效中值翻转点，随脉宽增加、温度升高或VCMA/SOT效率增大而左移；$$\beta_s$$控制曲线陡峭程度，热稳定因子越高、器件离散性越小则$$\beta_s$$越大，工艺波动导致的曲线展宽等效为$$\beta_s$$下降。从微观磁化动力学角度看，自旋翻转由随机Landau-Lifshitz-Gilbert（sLLG）方程主导，对应Fokker-Planck框架下概率团在双稳态势阱中的演化，严格求解的计算代价在SPICE仿真或阵列级评估中难以承受。Sigmoid近似将复杂的微观随机物理过程压缩为$$(u_{\mathrm{th}}, \beta_s)$$这两个兼具物理意义与工程可测性的拟合参数，输出天然位于$$[0]$$而无需额外裁剪，便于在SPICE紧凑模型与系统级仿真中直接调用。改写为$$P_{\mathrm{sw}}(u) \approx \sigma(\beta_s u + b)$$形式时，截距$$b = -\beta_s u_{\mathrm{th}}$$；增益$$\beta_s$$由热稳定因子$$\Delta$$、VCMA系数$$\xi$$、SHE效率$$\theta_{\mathrm{SH}}$$等底层物理参数共同决定，参数$$(\beta_s, u_{\mathrm{th}})$$因此可由工艺标定直接给出。

![面向概率计算的sMTJ行为级模型分层架构](figures/fig_02_behavioral_layers.png)

**图2.2：面向概率计算的sMTJ行为级模型分层架构。** 自上而下四个层次。(a)计算层：从外部驱动量$$(t_w, I_{SOT}, V_{MTJ})$$到等效驱动$$u$$的映射，以及基于sMTJ随机性的伯努利采样，输出由$$P_{\mathrm{sw}}(u)$$参数化的随机比特流，作为可编程随机源供后续电路与系统级使用。(b)行为层：临界过渡区内$$P_{\mathrm{sw}}$$的Sigmoid近似，由等效驱动量$$u$$、阈值$$u_{\mathrm{th}}$$与斜率$$\beta_s$$定义，并展示其与底层物理参数$$\Delta$$、$$\xi$$、$$\theta_{\mathrm{SH}}$$的依赖关系。(c)物理模型层：包含SOT与热力矩的LLG方程、综合VCMA能垒压低与SOT驱动的调制能垒$$\Delta E_b$$，以及用于推导统一翻转概率$$P_{\mathrm{sw}}$$的Néel–Brown模型。(d)器件层：具有PMA的sMTJ结构与三类物理驱动机制，包括施加于MgO势垒两端的VCMA偏压$$V_{MTJ}$$、流经HM沟道的SOT电流$$I_{SOT}$$以及热涨落。

## 参考文献

[1] L. Néel. Théorie du traînage magnétique des ferromagnétiques en grains fins avec application aux terres cuites. *Annales de Géophysique*, 5: 99–136, 1949.

[2] W. F. Brown. Thermal Fluctuations of a Single-Domain Particle. *Physical Review*, 130(5): 1677–1686, 1963. [https://doi.org/10.1103/PhysRev.130.1677](https://doi.org/10.1103/PhysRev.130.1677)

[3] W. F. Brinkman, R. C. Dynes, J. M. Rowell. Tunneling Conductance of Asymmetrical Barriers. *Journal of Applied Physics*, 41(5): 1915–1921, 1970. [https://doi.org/10.1063/1.1659141](https://doi.org/10.1063/1.1659141)

[4] J. S. Moodera, L. R. Kinder, T. M. Wong, R. Meservey. Large Magnetoresistance at Room Temperature in Ferromagnetic Thin Film Tunnel Junctions. *Physical Review Letters*, 74(16): 3273–3276, 1995. [https://doi.org/10.1103/PhysRevLett.74.3273](https://doi.org/10.1103/PhysRevLett.74.3273)

[5] L. Liu, O. J. Lee, T. J. Gudmundsen, D. C. Ralph, R. A. Buhrman. Current-Induced Switching of Perpendicularly Magnetized Magnetic Layers Using Spin Torque from the Spin Hall Effect. *Physical Review Letters*, 109(9): 096602, 2012. [https://doi.org/10.1103/PhysRevLett.109.096602](https://doi.org/10.1103/PhysRevLett.109.096602)

[6] V. Krizakova, M. Perumkunnil, S. Couet, P. Gambardella, K. Garello. Spin-orbit torque switching of magnetic tunnel junctions for memory applications. *Journal of Magnetism and Magnetic Materials*, 562: 169692, 2022. [https://doi.org/10.1016/j.jmmm.2022.169692](https://doi.org/10.1016/j.jmmm.2022.169692)

[7] T. Nozaki, J. Okabayashi, S. Tamaru, *et al.* Understanding voltage-controlled magnetic anisotropy effect at Co/oxide interface. *Scientific Reports*, 13: 10640, 2023. [https://doi.org/10.1038/s41598-023-37422-4](https://doi.org/10.1038/s41598-023-37422-4)

[8] Y. Lv, B. Dixit, J.-P. Wang. Modulation of switching dynamics in magnetic tunnel junctions for low-error-rate computational random-access memory. *AIP Advances*, 16(2): 025134, 2026. [https://doi.org/10.1063/9.0001026](https://doi.org/10.1063/9.0001026)

[2]:
