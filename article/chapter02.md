# 第2章 三端SOT-sMTJ器件建模、仿真与实验验证

第一章所提出的全自旋三位一体架构与时域展开范式以sMTJ作为唯一的硬件单元，把存储、随机源与乘加三类功能统一在同一MRAM工艺之下。该架构能否落地，首先取决于sMTJ自身能否提供一个可参数化、可工艺标定，并且与算法语义直接对接的Bernoulli采样接口；否则上层的伊辛求解与PBNN推断都只能停留在概念可行而工程不成立的层面。本章因此把第三、四章共同依赖的器件层物理基础前置完成，将三端SOT-sMTJ从器件端口、磁化动力学、统计翻转行为与实验标定四个层面接成同一套可调用模型，使后续两类任务在器件层共享单一、可信、可复现的物理基础。

为此，2.1节从T型电路与SOT-VCMA联合驱动出发建立概率翻转的器件级抽象，明确读写解耦与Sigmoid参数化的物理来源；2.2节以含热噪声、自热反馈和温度依赖材料参数的sLLG求解流程说明该抽象的微观依据，并以开源仿真平台vgsot-sim给出可调用的工程实现；2.3节用300 mm工艺平台上的单器件实测与PDK失配数据校准Néel-Brown与Sigmoid参数，并提出D2D-C2C双层正交分解框架以分离热激活物理与分布形态偏差；2.4节归纳后续伊辛与PBNN两章可直接继承的Sigmoid接口、五参数行为模型与工艺裕度结论。本章工作的定位是为整条研究主线提供唯一的物理入口——后两章对sMTJ的任何调用均不再重复底层物理推导，只需经由本章所定义的接口接入。

## 2.1三端SOT-sMTJ器件模型

sMTJ的器件层模型需要在静态电学特性、自由层磁化动力学、写入概率以及行为级抽象四个层次之间保持自洽，才能够支撑后续阵列级与系统级的建模工作。本节首先建立三端SOT-sMTJ的T型等效电路并阐明读写路径解耦的物理意义，随后由Néel–Brown热激活模型给出受驱动调制的翻转概率，综合SOT与VCMA两类驱动机制得到联合写入概率表达式，最终将复杂的指数嵌套形式压缩为工作区内便于工程使用的Sigmoid近似。[^ref-brown-thermal]

### 2.1.1三端SOT-sMTJ器件结构与基础电学模型

三端SOT-sMTJ由垂直磁各向异性 (PMA) 的MTJ堆叠 (参考层/隧穿势垒/自由层) 置于重金属 (HM) 沟道之上构成。MTJ顶端为读出端$$T_1$$，沟道两端分别为$$T_2$$与$$T_3$$。该拓扑使读出与写入路径在物理上分离：读出操作仅在$$T_1$$与沟道接地端之间施加低偏压，依靠TMR效应产生与磁化态相关的隧穿电流；写入操作则将大电流沿$$T_2$$至$$T_3$$流过HM沟道，由自旋霍尔效应产生的横向自旋积累驱动自由层翻转。读写解耦避免了STT-MTJ读写共用同一隧穿路径所引起的可靠性退化问题[^ref-liu-sot-prl]，是SOT-sMTJ作为可编程随机源的结构基础。

将MTJ等效为电阻$$R_{MTJ}$$、HM沟道沿长度方向均分为两段电阻$$R_{SOT}/2$$，三端器件可建模为T型电阻网络，其中沟道电阻由几何与材料参数给定

$$
R_{SOT} = \frac{\rho_{SOT} L_{SOT}}{W_{SOT} T_{SOT}},
$$

$$\rho_{SOT}$$、$$L_{SOT}$$、$$W_{SOT}$$与$$T_{SOT}$$分别为HM沟道的电阻率、长度、宽度与厚度。

![三端SOT-sMTJ器件结构与T型等效电路网络](figs/Chapter02_local_01.png)

**图2.1：三端SOT-sMTJ器件结构与T型等效电路网络。** (a)器件物理堆叠结构，包括位于顶端的MTJ (参考层/MgO势垒/自由层) 与底部的HM沟道；MTJ顶电极对应端口$$T_1$$，HM沟道两端分别引出端口$$T_2$$与$$T_3$$。(b)T型电阻网络，$$R_{MTJ}$$表示自由层与参考层之间的隧穿电阻，沟道被等分为两段$$R_{SOT}/2$$分别接$$T_2$$与$$T_3$$，三段电阻在中间节点N处汇合。读出操作通过$$T_1$$施加小偏压获取$$R_{MTJ}$$状态，写入操作通过$$T_2$$、$$T_3$$沿沟道施加大电流产生SOT，两条电流路径在物理上解耦。

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

该表达式从电路层面验证了三端结构的读写解耦：施加于$$T_1$$的MTJ偏压$$V_{MTJ} = V_1 - V_N$$即使较大也不会改变$$I_{SOT}$$的方向与幅值，从而允许将$$V_{MTJ}$$独立用作VCMA能垒调制变量 (详见2.1.3节)。

MTJ的隧穿电阻由Brinkman模型描述，该模型基于WKB近似处理梯形势垒的电压依赖隧穿[^ref-brinkman-bdr]。在低偏压下，平行态电阻可表为势垒厚度与势垒高度的函数

$$
R_P = \frac{t_{ox}}{F \cdot A_{MTJ} \cdot \phi_{ox}}\exp\!\left(\frac{2 t_{ox}}{\hbar}\sqrt{2 m_e e\, \phi_{ox}}\right),
$$

其中$$t_{ox}$$与$$\phi_{ox}$$分别为隧穿势垒层厚度与势垒高度，$$F$$为根据电阻—面积乘积$$RA$$标定的拟合因子，$$A_{MTJ}$$为MTJ截面积，$$e$$为元电荷，$$m_e$$为电子质量，$$\hbar$$为约化普朗克常数。指数因子主导了$$R_P$$对势垒厚度的强敏感性，$$t_{ox}$$每变化$$0.1\,\mathrm{nm}$$量级即可使$$R_P$$变化数倍，构成MTJ阻值工艺难以精细控制的根本来源之一。

反平行态电阻$$R_{AP}$$由TMR比定义。在Jullière二流模型下，理论TMR可由两侧铁磁电极的有效自旋极化率$$P$$给出

$$
\mathrm{TMR}_0 = \frac{R_{AP}-R_P}{R_P} = \frac{2P^2}{1-P^2}.
$$

实际TMR随MTJ两端电压$$V_{MTJ}$$增大而衰减，主要源自界面处的非弹性磁子辅助隧穿过程，即磁子吸收/发射通道随偏压打开后破坏了自旋守恒[^ref-moodera-tmr]。Moodera等的实验观测显示，TMR随$$V_{MTJ}$$近似呈Lorentzian形衰减，常用的工程参数化形式为

$$
\mathrm{TMR}(V_{MTJ}) = \frac{\mathrm{TMR}_0}{1 + V_{MTJ}^2/V_h^2},
$$

其中特征电压$$V_h$$定义为TMR下降至理论值一半时所对应的偏压，CoFeB/MgO/CoFeB体系中典型值约为$$0.3\text{--}0.7\,\mathrm{V}$$。当自由层与参考层磁化方向夹角为$$\theta$$时，MTJ电阻满足如下含偏压依赖的角度公式

$$
R_{MTJ}(V_{MTJ},\theta) = R_P\,\frac{1+V_{MTJ}^2/V_h^2 + \mathrm{TMR}_0}{1+V_{MTJ}^2/V_h^2 + 0.5\,\mathrm{TMR}_0(1+\cos\theta)}.
$$

在$$\theta = 0$$ (平行) 时分子分母相同，$$R_{MTJ} = R_P$$；在$$\theta = \pi$$ (反平行) 且$$V_{MTJ} \ll V_h$$时，$$R_{MTJ} \to R_P(1+\mathrm{TMR}_0) = R_{AP}$$。该公式同时刻画了TMR随偏压的压缩与随磁化夹角的连续变化，是磁动力学求解过程中读出电流计算的基础。

读路径电压$$V_{MTJ}$$与写路径电流$$I_{SOT}$$在端口拓扑层面已实现解耦，但二者在磁化层内部仍通过共享自由层的能垒调制相互耦合：$$V_{MTJ}$$经VCMA压低能垒，$$I_{SOT}$$经SHE提供翻转驱动力，两者协同决定动态翻转概率。这一耦合关系将在下两节中具体展开。

### 2.1.2自由层磁化动力学与热激活翻转概率

为给出MTJ在不同驱动条件下的翻转概率，本节采用宏自旋近似下的紧凑模型展开推导，主要包含Néel–Brown热激活模型与受驱动调制的能垒缩减表达式。

对具有PMA的MTJ自由层，热稳定因子定义为$$\Delta = E_b/(k_B T)$$，度量两稳定磁化态之间的能垒高度，其中$$E_b = K_{\mathrm{eff}}V \approx \tfrac{1}{2}\mu_0 H_k M_s V$$，$$V$$为自由层体积，$$K_{\mathrm{eff}}$$为有效各向异性能密度，$$H_k$$为等效各向异性场，$$M_s$$为饱和磁化强度。在亚临界写入条件下，sMTJ的翻转受热涨落主导，由Néel–Brown热激活模型描述[^ref-brown-thermal]：磁化在能垒约束下发生热激活逃逸，平均停留时间满足Arrhenius关系$$\tau = \tau_0 \exp(\Delta)$$，其中$$\tau_0$$为亚纳秒至纳秒量级的尝试时间，典型取值约为$$1\,\mathrm{ns}$$。在持续时间为$$t_w$$的写脉冲下以泊松过程近似，翻转概率为

$$
P_{\mathrm{sw}}(t_w) = 1 - \exp\!\left(-\frac{t_w}{\tau_0 e^{\Delta}}\right).
$$

外场、自旋转移矩或SOT驱动存在时，常将其等效为能垒压低$$\Delta_{\mathrm{eff}} = \Delta(1-I/I_{c0})^n$$ ($$I < I_{c0}$$)，其中$$I_{c0}$$为零温临界驱动，指数$$n$$由器件物理工作区间决定：亚临界热激活区 (驱动远低于$$I_{c0}$$、热涨落主导翻转) 能垒随驱动呈线性压低，对应$$n = 1$$；接近零温弹道极限 (驱动接近$$I_{c0}$$、热涨落可忽略) 能垒压低呈抛物线行为，对应$$n = 2$$。概率计算工作区位于亚临界热激活区，故下文取$$n = 1$$。综合得带驱动的翻转概率

$$
P_{\mathrm{sw}}(t_w, I) = 1 - \exp\!\left[-\frac{t_w}{\tau_0}\exp\!\left(-\Delta\!\left(1-\frac{I}{I_{c0}}\right)^n\right)\right].
$$

定义二值随机变量$$m \in \{0,1\}$$表示器件翻转事件，则$$m \sim \operatorname{Bernoulli}(P_{\mathrm{sw}})$$，sMTJ自然成为以驱动量为参数的可编程伯努利源。该模型在概率计算工作区与器件实测的吻合度将在2.3节通过专门的实验对比给出，本节不再赘述。

### 2.1.3 SOT与VCMA联合驱动模型

SOT-MRAM的写入电流流经HM沟道，由自旋霍尔效应或Rashba–Edelstein效应产生横向自旋积累，对自由层施加阻尼型等效场

$$
H_{\mathrm{DL}} = \frac{\hbar \theta_{\mathrm{SH}} I_{SOT}}{2 e \mu_0 M_s t_f A_c},
$$

其中$$\theta_{\mathrm{SH}}$$为有效自旋霍尔角，$$I_{SOT}$$即2.1.1节由$$(V_2-V_3)/R_{SOT}$$给出的沟道贯穿电流，$$t_f$$为自由层厚度，$$A_c$$为沟道截面积。零温近似下临界SOT电流标度为

$$
I_{c0}^{\mathrm{SOT}} \propto \frac{2e}{\hbar}\,\frac{\alpha M_s t_f A_c}{\theta_{\mathrm{SH}}}\,H_k,
$$

降低$$M_s$$、$$t_f$$、阻尼系数$$\alpha$$或提高$$\theta_{\mathrm{SH}}$$均可减小写入电流[^ref-krizakova-sot-review]。

VCMA机制通过MTJ两端的偏压$$V_{MTJ}$$改变自由层/氧化层界面的各向异性，界面各向异性能面密度$$K_i$$随电场线性变化[^ref-nozaki-vcma]

$$
K_i(V_{MTJ}) = K_i(0) - \xi\,\frac{V_{MTJ}}{t_{ox}},
$$

其中$$\xi$$为VCMA系数 (量纲$$\mathrm{fJ/(V\cdot m)}$$)，$$t_{ox}$$已由2.1.1节给出。此处$$V_{MTJ}$$正是2.1.1节T型等效电路中$$T_1$$与中点N之间的电压降，电学模型与磁动力学模型由此通过共享变量自然衔接。将界面项折算到体各向异性后，能垒变为

$$
E_b(V_{MTJ}) = \left[K_{\mathrm{eff}}(0) - \frac{\xi}{t_f}\,\frac{V_{MTJ}}{t_{ox}}\right]V,
$$

代入Néel–Brown模型可得VCMA辅助下的翻转概率

$$
P_{\mathrm{sw}}(t_w, V_{MTJ}) = 1 - \exp\!\left(-\frac{t_w}{\tau_0}\exp\!\left[-\Delta_0 + \beta_{\mathrm{VCMA}} V_{MTJ}\right]\right),
$$

其中$$\Delta_0 = E_b(0)/(k_B T)$$，$$\beta_{\mathrm{VCMA}} = \xi V/(k_B T \cdot t_f t_{ox})$$。VCMA将翻转概率曲线沿$$V_{MTJ}$$方向左移，使较低写入能量下即可获得可控的随机翻转。

VCMA与SOT联合作用时，自然的建模顺序是VCMA先降低各向异性与能垒 ($$H_k \to H_k(V_{MTJ})$$)，SOT再提供翻转驱动力，由此得到偏压依赖的SOT临界电流$$I_{c0}^{\mathrm{VCMA}}(V_{MTJ})$$及联合翻转概率

$$
P_{\mathrm{sw}}(t_w, I_{SOT}, V_{MTJ}) = 1 - \exp\!\left[-\frac{t_w}{\tau_0}\exp\!\left(-(\Delta_0 - \beta_{\mathrm{VCMA}}V_{MTJ})\!\left(1-\frac{I_{SOT}}{I_{c0}^{\mathrm{VCMA}}(V_{MTJ})}\right)^n\right)\right].
$$

由于$$H_k(V_{MTJ})$$随VCMA偏压下降，$$I_{c0}^{\mathrm{VCMA}}$$亦相应降低，VCMA辅助SOT得以减小写电流与写能耗[^ref-zhang-vgsot]。该联合模型给出了三端sMTJ作为可编程伯努利源的完整描述：外部电路只需指定脉宽$$t_w$$、SOT电流$$I_{SOT}$$与MTJ偏压$$V_{MTJ}$$这一三元组，器件即输出对应概率的随机比特。

### 2.1.4工作区Sigmoid近似

联合翻转概率$$P_{\mathrm{sw}}$$是驱动量的指数嵌套函数，严格数学形式属于Gumbel型累积分布而非标准Sigmoid (Logistic) 函数。但在临界过渡区 ($$P_{\mathrm{sw}} \approx 0.05$$至$$0.95$$) 高阶项影响甚小，可作Sigmoid近似。设等效驱动变量$$u$$综合脉宽、SOT电流与VCMA偏压的影响，将有效能垒在工作点$$u_0$$附近一阶展开$$\Delta E_{\mathrm{eff}}(u) \approx \Delta E_{\mathrm{eff}}(u_0) - \kappa(u-u_0)$$，其中$$\kappa = -\partial \Delta E_{\mathrm{eff}}/\partial u|_{u_0}$$为能垒压低效率，则

$$
P_{\mathrm{sw}}(u) = 1 - \exp[-A\exp(B(u-u_0))], \qquad B = \frac{\kappa}{k_B T}.
$$

定义中值翻转点$$u_{\mathrm{th}}$$满足$$P_{\mathrm{sw}}(u_{\mathrm{th}}) = 0.5$$，将$$P_{\mathrm{sw}}$$转换到logit域$$L(u) = \ln[P_{\mathrm{sw}}/(1-P_{\mathrm{sw}})]$$并在$$u_{\mathrm{th}}$$处一阶泰勒展开 (零次项$$L(u_{\mathrm{th}}) = 0$$)，得$$L(u) \approx 2B\ln 2 \cdot (u-u_{\mathrm{th}})$$。记$$\beta_s = 2\kappa\ln 2/(k_B T)$$，还原概率形式即得过渡区Sigmoid近似

$$
P_{\mathrm{sw}}(u) \approx \sigma(\beta_s(u-u_{\mathrm{th}})) = \frac{1}{1 + e^{-\beta_s(u-u_{\mathrm{th}})}}.
$$

参数$$u_{\mathrm{th}}$$为等效中值翻转点，随脉宽增加、温度升高或VCMA/SOT效率增大而左移；$$\beta_s$$控制曲线陡峭程度，热稳定因子越高、器件离散性越小则$$\beta_s$$越大，工艺波动导致的曲线展宽等效为$$\beta_s$$下降。从微观磁化动力学角度看，自旋翻转由随机Landau-Lifshitz-Gilbert (sLLG) 方程主导，对应Fokker-Planck框架下概率团在双稳态势阱中的演化，严格求解的计算代价在SPICE仿真或阵列级评估中难以承受。Sigmoid近似将复杂的微观随机物理过程压缩为$$(u_{\mathrm{th}}, \beta_s)$$这两个兼具物理意义与工程可测性的拟合参数，输出天然位于$$[0,1]$$而无需额外裁剪，便于在SPICE紧凑模型与系统级仿真中直接调用。改写为$$P_{\mathrm{sw}}(u) \approx \sigma(\beta_s u + b)$$形式时，截距$$b = -\beta_s u_{\mathrm{th}}$$；增益$$\beta_s$$由热稳定因子$$\Delta$$、VCMA系数$$\xi$$、SHE效率$$\theta_{\mathrm{SH}}$$等底层物理参数共同决定，参数$$(\beta_s, u_{\mathrm{th}})$$因此可由工艺标定直接给出。

![面向概率计算的sMTJ行为级模型分层架构](figs/Chapter02_local_02.png)

**图2.2：面向概率计算的sMTJ行为级模型分层架构。** 自上而下四个层次。(a)计算层：从外部驱动量$$(t_w, I_{SOT}, V_{MTJ})$$到等效驱动$$u$$的映射，以及基于sMTJ随机性的伯努利采样，输出由$$P_{\mathrm{sw}}(u)$$参数化的随机比特流，作为可编程随机源供后续电路与系统级使用。(b)行为层：临界过渡区内$$P_{\mathrm{sw}}$$的Sigmoid近似，由等效驱动量$$u$$、阈值$$u_{\mathrm{th}}$$与斜率$$\beta_s$$定义，并展示其与底层物理参数$$\Delta$$、$$\xi$$、$$\theta_{\mathrm{SH}}$$的依赖关系。(c)物理模型层：包含SOT与热力矩的LLG方程、综合VCMA能垒压低与SOT驱动的调制能垒$$\Delta E_b$$，以及用于推导统一翻转概率$$P_{\mathrm{sw}}$$的Néel–Brown模型。(d)器件层：具有PMA的sMTJ结构与三类物理驱动机制，包括施加于MgO势垒两端的VCMA偏压$$V_{MTJ}$$、流经HM沟道的SOT电流$$I_{SOT}$$以及热涨落。

## 2.2 磁学仿真方法与平台

在基于MRAM的概率计算体系中，器件层面的随机磁化动力学决定了计算行为的本质。与传统确定性存储器不同，sMTJ在纳米尺度下的磁化翻转过程不可避免地受到热涨落影响，从而表现出显著的cycle-to-cycle (C2C) 随机性。这种C2C随机性在概率计算框架中被主动利用，以用于实现硬件级随机采样与概率推断。前文已经从统计角度给出了基于Néel–Brown模型的sMTJ翻转概率行为模型，本节则从磁化动力学出发，通过sLLG方程与磁场热噪声项建模推导出概率行为的物理起源，从而实现从微观动力学到宏观概率模型的统一描述。

磁化动力学的建模通常可以分为宏自旋模型与微磁模型两种层级。宏自旋模型假设自由层磁化在空间上均匀，仅随时间演化，适用于尺寸较小、接近单畴的器件结构；微磁模型则进一步引入空间离散化与交换相互作用，能够描述复杂磁畴结构与非均匀动力学行为。当器件物理直径$$D_{\mathrm{phys}}$$小于单畴临界尺寸$$l_{\mathrm{ex}} = \sqrt{2A_{\mathrm{ex}}/(\mu_0 M_s^2)}$$时，自由层可被合理地视为空间均匀的单一磁矩[^ref-kittel-domain]。考虑到本文所研究的80 nm级MTJ器件满足该条件，同时研究重点在于概率翻转统计特性与电路级耦合效率，因此采用宏自旋模型作为基础，并在其上引入热噪声与VCMA调制项，从而在保证物理准确性的同时兼顾计算效率。

全节安排如下：2.2.1节建立包含STT、SOT与VCMA项的扩展LLG方程，并推导有效磁场的各物理分量；2.2.2节在此基础上补充温度效应与器件非理想性，包括自热建模、材料参数温度依赖、TMR非线性以及退磁形状效应；2.2.3节介绍适用于随机微分方程的保模长数值积分算法与Monte Carlo统计方法；2.2.4节介绍与上述理论框架配套的开源仿真平台vgsot-sim；2.2.5节对基于该平台的数值仿真流程与可输出观测量给出形式化说明。

---

### 2.2.1 磁化动力学基础模型

#### 2.2.1.1 Landau–Lifshitz–Gilbert方程与自旋力矩项

自由层磁化矢量$$\mathbf{m}$$ (已归一化，$$|\mathbf{m}|=1$$) 的时间演化由扩展的Landau–Lifshitz–Gilbert (LLG) 方程描述。该方程最初由Landau和Lifshitz于1935年在铁磁共振理论框架下提出[^ref-landau-lifshitz]，Gilbert于1955年引入耗散项并将其改写为与实验更为吻合的阻尼形式[^ref-gilbert-damping]，后经Slonczewski[^ref-slonczewski-stt]与Berger[^ref-berger-stt]在1996年分别独立地将自旋转移力矩 (STT) 引入；本文所关注的SOT-MTJ体系进一步纳入了由重金属层自旋霍尔效应产生的自旋轨道力矩 (SOT) 项[^ref-miron-sot][^ref-liu-spin-hall]。其完整隐式形式为

$$
\frac{\partial \mathbf{m}}{\partial t}
= -\gamma \mathbf{m} \times \mathbf{H}_{\mathrm{eff}}
+ \alpha \mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t}
+ \boldsymbol{\tau}_{\mathrm{STT}}
+ \boldsymbol{\tau}_{\mathrm{SOT}}
$$

其中$$\gamma$$为旋磁比，$$\alpha$$为Gilbert阻尼系数，$$\mathbf{H}_{\mathrm{eff}}$$为有效磁场 (详见2.2.1.2节)，$$\boldsymbol{\tau}_{\mathrm{STT}}$$、$$\boldsymbol{\tau}_{\mathrm{SOT}}$$分别为两类自旋力矩项。本小节先给出二者的物理表达，再给出显式数值求解形式。

**STT项。** 自旋转移力矩源自穿过MTJ势垒的隧穿电流中的极化电子，经过自由层时把自旋角动量转移给局域磁矩。以自由层归一化磁化矢量$$\mathbf{m}$$与参考层固定磁化方向$$\hat{\mathbf{m}}_p$$表示，Slonczewski与Berger独立推导得到阻尼型分量
$$
\boldsymbol{\tau}_{\mathrm{STT,DL}}
= -\gamma\,a_J\,\mathbf{m}\times(\mathbf{m}\times\hat{\mathbf{m}}_p),
\qquad
a_J \equiv \frac{\hbar P J_{\mathrm{STT}}}{2 e \mu_0 M_s t_f}
$$

其中$$J_{\mathrm{STT}}=I_{\mathrm{MTJ}}/A_{\mathrm{MTJ}}$$为隧穿电流密度，$$P$$为自旋极化率 (CoFeB/MgO界面取约0.58[^ref-ikeda-pma])。$$a_J$$具有磁场量纲 (A/m)，等价于将STT折算到与PMA同形式的等效阻尼型场$$H_{\mathrm{DL}}^{\mathrm{STT}}=a_J$$。同时存在共轭的场型分量

$$
\boldsymbol{\tau}_{\mathrm{STT,FL}} = +\gamma\,b_J\,\mathbf{m}\times\hat{\mathbf{m}}_p,
\qquad
b_J = \beta_{\mathrm{STT}}^{\mathrm{FL/DL}}\,a_J,
$$

工程上以比例系数$$\beta_{\mathrm{STT}}^{\mathrm{FL/DL}}\in[0,0.3]$$参数化；对垂直MTJ这一比值偏小，本工作中按0处理。注意$$\mathbf{m}\times\hat{\mathbf{m}}_p$$在$$\mathbf{m}$$靠近极轴 ($$\hat{\mathbf{m}}_p=\pm\hat{\mathbf{z}}$$) 时随$$\sin\theta$$线性消失，因此纯STT写入需要依赖热涨落把磁化从极轴附近偏离才能启动翻转。

**SOT项。** 自旋轨道力矩由重金属层 (如Ta、$$\beta$$-W) 中的自旋霍尔效应产生：横向电流通过自旋霍尔角$$\theta_{\mathrm{SH}}$$转换为纵向自旋流，在重金属/铁磁层界面累积出极化方向$$\boldsymbol{\sigma}=\hat{\mathbf{z}}\times\hat{\mathbf{j}}$$的自旋极化，进而对自由层磁矩施加力矩。该力矩按物理机制分解为阻尼型 (DL) 与场型 (FL) 两个分量[^ref-manchon-sot]：
$$
\boldsymbol{\tau}_{\mathrm{SOT}}
= -\gamma H_{\mathrm{DL}}^{\mathrm{SOT}}\,\mathbf{m}\times(\mathbf{m}\times\boldsymbol{\sigma})
+ \gamma H_{\mathrm{FL}}^{\mathrm{SOT}}\,\mathbf{m}\times\boldsymbol{\sigma}
$$

其中DL分量的等效场幅值由自旋霍尔角$$\theta_{\mathrm{SH}}$$、沟道电流密度$$J_{\mathrm{SOT}}=I_{\mathrm{SOT}}/A_{\mathrm{SOT}}$$以及自由层厚度$$t_f$$共同决定：

$$
H_{\mathrm{DL}}^{\mathrm{SOT}}
= \frac{\hbar \theta_{\mathrm{SH}} J_{\mathrm{SOT}}}{2 e \mu_0 M_s t_f},
\qquad
H_{\mathrm{FL}}^{\mathrm{SOT}} = \beta_{\mathrm{SOT}}^{\mathrm{FL/DL}}\,H_{\mathrm{DL}}^{\mathrm{SOT}}.
$$

DL分量与Gilbert阻尼等价但方向相反，超过临界电流密度后可驱动磁化发生确定性翻转或进入概率性翻转窗口。与STT不同，$$\mathbf{m}\times(\mathbf{m}\times\boldsymbol{\sigma})$$在$$\boldsymbol{\sigma}$$与极轴正交 (如$$\boldsymbol{\sigma}\perp\hat{\mathbf{z}}$$) 时即便$$\mathbf{m}\parallel\hat{\mathbf{z}}$$也保持有限，故SOT写入无需热涨落即可在极轴附近启动，亚纳秒级翻转因此成为可能。FL分量等效为附加偏置场，对翻转阈值有修正作用；对称界面结构中其幅值通常远小于DL分量，但在某些重金属-铁磁体系中FL/DL比值可达$$0.5\sim1$$，仍须显式纳入。

**STT与SOT的对照。** 二者在器件依赖性与工作机制上的核心差异由下表归纳。

**表2.1** STT与SOT两种自旋力矩驱动机制的对照表

| 物理量 | STT (Slonczewski) | SOT (Spin Hall) |
|---|---|---|
| 驱动电流路径 | $$I_{\mathrm{MTJ}}=V_{\mathrm{MTJ}}/R_{\mathrm{MTJ}}$$，纵贯隧道势垒 | $$I_{\mathrm{SOT}}=(V_2-V_3)/R_{\mathrm{SOT}}$$，沿重金属沟道横向 |
| 阻尼型等效场 | $$H_{\mathrm{DL}}^{\mathrm{STT}}=\hbar P J_{\mathrm{STT}}/(2e\mu_0 M_s t_f)$$ | $$H_{\mathrm{DL}}^{\mathrm{SOT}}=\hbar\theta_{\mathrm{SH}}J_{\mathrm{SOT}}/(2e\mu_0 M_s t_f)$$ |
| 效率因子 | 自旋极化率$$P\sim 0.58$$ | 自旋霍尔角$$\theta_{\mathrm{SH}}$$ ($$\beta$$-W体系约0.04–0.4，与厚度强相关) |
| 自旋极化方向 | $$\hat{\mathbf{m}}_p$$ (参考层磁化方向) | $$\boldsymbol{\sigma}=\hat{\mathbf{z}}\times\hat{\mathbf{j}}$$ (由电流方向决定) |
| 极轴处力矩 | $$\propto\sin\theta\to 0$$ (启动需热涨落) | 始终有限 (无需启动阈值) |
| 翻转能耗主因 | 隧穿耗散$$\propto V_{\mathrm{MTJ}}^2/R_{\mathrm{MTJ}}$$ | 沟道焦耳热$$\propto V_{\mathrm{SOT}}^2/R_{\mathrm{SOT}}$$ |
| 写入–读取通道 | 共用 (写电流流过MTJ) | 解耦 (写流沟道，读流势垒) |

本文以SOT为单器件主驱动机制：SOT力矩在极轴附近不消失，写入延迟可压至亚纳秒级；读写通道解耦避免了STT写入时势垒长期承受高偏压所引发的可靠性退化；三端结构还允许VCMA偏压与SOT电流独立调度，从而支持2.1.3节给出的SOT-VCMA联合驱动模型。STT项在本框架中保留为可选支路，用于STT-only写入能耗基准与STT+SOT联合校核。

**显式数值形式。** 在数值求解中，直接对隐式LLG积分会导致迭代步骤复杂、计算代价高昂。利用矢量恒等式$$\mathbf{m}\times(\mathbf{m}\times\mathbf{H})=(\mathbf{m}\cdot\mathbf{H})\mathbf{m}-\mathbf{H}$$以及$$|\mathbf{m}|=1$$的约束，可将其改写为不含$$\partial\mathbf{m}/\partial t$$隐式项的显式Landau–Lifshitz–Slonczewski (LLS) 形式：

$$
\frac{d\mathbf{m}}{dt}
= -\frac{\gamma}{1+\alpha^2}
\left[
\mathbf{m}\times\mathbf{H}_{\mathrm{eff}}
+ \alpha\,\mathbf{m}\times(\mathbf{m}\times\mathbf{H}_{\mathrm{eff}})
\right]
+ \mathbf{T}_{\mathrm{ext}}
$$

其中$$\mathbf{T}_{\mathrm{ext}} = \boldsymbol{\tau}_{\mathrm{STT}}+\boldsymbol{\tau}_{\mathrm{SOT}}$$为经$$(1+\alpha^2)$$预因子修正后的全部外加力矩之和，该因子由Gilbert阻尼项的消去过程引入。该显式形式将进动项与阻尼项分离，避免了隐式方程在大$$\alpha$$极限下的数值不稳定问题。在实际仿真中，时间步长须满足$$\Delta t\ll(\gamma H_{\mathrm{eff}})^{-1}$$以保证积分稳定性；具体的保模长几何积分算法在2.2.3节详述。

---

#### 2.2.1.2 有效磁场建模

有效磁场$$\mathbf{H}_{\mathrm{eff}}$$由多个物理机制叠加而成：

$$
\mathbf{H}_{\mathrm{eff}}
= \mathbf{H}_{\mathrm{PMA}} + \mathbf{H}_{\mathrm{VCMA}}
+ \mathbf{H}_{\mathrm{D}} + \mathbf{H}_{\mathrm{EX}} + \mathbf{H}_{\mathrm{TH}}
$$

以下分别给出各分量的物理起源与数学表达。

**垂直磁各向异性场 (PMA)。** 在CoFeB/MgO界面体系中，由于Fe-O键的轨道杂化，界面处积累的垂直各向异性能可在较薄的自由层中克服形状各向异性，从而使易磁化轴沿垂直于薄膜平面的$$z$$方向稳定。其等效场为
$$
\mathbf{H}_{\mathrm{PMA}}
= \frac{2K_i}{\mu_0 M_s t_f}\,m_z\,\hat{z}
$$

其中$$K_i$$为单位面积界面各向异性能密度，$$M_s$$为饱和磁化强度，$$t_f$$为自由层厚度。该场沿$$z$$轴方向对磁化产生恢复力，是PMA-MTJ中能垒的主要来源。

**VCMA等效场。** 电压调控磁各向异性 (VCMA) 效应是指在MgO势垒两端施加电压时，界面电场改变Fe-O键的轨道占据，从而调制$$K_i$$的大小[^ref-nozaki-vcma-feb]。在线性响应范围内，VCMA等效磁场可写为

$$
\mathbf{H}_{\mathrm{VCMA}}
= -\frac{2\beta_{\mathrm{VCMA}}\,V_{\mathrm{MTJ}}}{\mu_0 M_s\,t_{\mathrm{ox}}\,t_f}\,m_z\,\hat{z}
$$

其中$$\beta_{\mathrm{VCMA}}$$为VCMA系数 (单位fJ·V$$^{-1}$$·m$$^{-1}$$)，$$V_{\mathrm{MTJ}}$$为MTJ两端电压，$$t_{\mathrm{ox}}$$为MgO势垒厚度。$$\mathbf{H}_{\mathrm{VCMA}}$$与$$\mathbf{H}_{\mathrm{PMA}}$$同向叠加，即正向电压降低等效各向异性场、削减能垒，这一特性是VGSOT写入方案中实现低功耗辅助翻转的物理基础。

**退磁场。** 有限尺寸的磁性薄层在磁化过程中会产生与磁化方向相反的退磁场，其表达式为
$$
\mathbf{H}_{\mathrm{D}} = -M_s\,\mathbf{N}\cdot\mathbf{m}
$$

其中$$\mathbf{N}$$为退磁张量。对于圆柱形MTJ自由层 ($$t_f \ll D_{\mathrm{phys}}$$)，在薄圆盘极限下退磁因子可近似为

$$
N_x = N_y \approx \frac{\pi t_f}{4D_{\mathrm{phys}}},\qquad
N_z = 1 - 2N_x
$$

需要指出的是，上述退磁因子表达式是对$$t_f/D_{\mathrm{phys}} \ll 1$$极限下扁椭球体的线性化近似，适用于快速定性估算。当器件尺寸缩减或自由层厚度增大、宽厚比不满足薄圆盘近似时，应采用精确的椭球体解析公式，具体形式将在2.2.2.4节给出，实际仿真中亦采用该精确式。此外，该近似表达式中的直径取物理直径$$D_{\mathrm{phys}}$$，而电学有效直径$$D_{\mathrm{elec}}$$因边缘刻蚀效应通常较$$D_{\mathrm{phys}}$$小约5–10 nm，两者的区分在2.2.2.4节一并讨论。

**交换偏置场。** 为实现无外加磁场条件下的确定性场无关SOT翻转，VGSOT结构中引入了合成反铁磁 (SAF) 层或直接的反铁磁钉扎层以提供面内交换偏置场：

$$
\mathbf{H}_{\mathrm{EX}} = H_{\mathrm{EX}}\,\hat{y}
$$

该偏置场打破了SOT翻转中$$\pm z$$方向的等效对称性，使驱动电流的极性与翻转方向之间形成确定性的一一对应关系[^note-dev-hex]。

**热噪声场。** 在有限温度下，自由层磁矩与晶格声子系统之间的热涨落通过涨落-耗散定理 (fluctuation-dissipation theorem) 耦合进入LLG方程[^ref-callen-welton]。Brown于1963年在单畴粒子框架下严格证明，与Gilbert阻尼$$\alpha$$共轭的热随机场$$\mathbf{H}_{\mathrm{TH}}$$必须满足白噪声统计，其二阶相关函数为[^ref-brown-thermal]
$$
\langle H_{\mathrm{TH},i}(t)\,H_{\mathrm{TH},j}(t')\rangle
= \delta_{ij}\,\delta(t-t')\,
\frac{2k_BT\alpha}{\mu_0 M_s \gamma V}
$$

其中$$k_B$$为Boltzmann常数，$$T$$为器件瞬态温度，$$V = \pi D_{\mathrm{phys}}^2 t_f/4$$为自由层体积，$$i,j \in \{x,y,z\}$$。上述二阶相关函数表明热噪声场在各分量间互不相关 ($$\delta_{ij}$$)，且在时间上为白噪声 ($$\delta(t-t')$$)，其幅度正比于$$\sqrt{\alpha T/(M_s\gamma V)}$$，即较大的阻尼、较高的温度或较小的磁矩体积均会增强热涨落强度。

在数值离散化实现中，为与有限时间步$$\Delta t$$相容，对连续白噪声功率谱密度进行时域积分并在每步独立重新采样，得到各分量的等效高斯随机场幅值：

$$
\mathbf{H}_{\mathrm{TH}}
= \boldsymbol{\xi}\sqrt{\frac{2k_BT\alpha}{\mu_0 M_s \gamma V \Delta t}}
$$

其中$$\boldsymbol{\xi} = (\xi_x, \xi_y, \xi_z)^{\mathrm{T}}$$，各分量为相互独立的标准正态随机变量，即$$\xi_i \sim \mathcal{N}(0,1)$$。离散采样形式等价于以时间步$$\Delta t$$对原始相关函数中的$$\delta(t-t')$$进行积分后得到的离散版本，García-Palacios和Lázaro于1998年对此给出了详细推导，并验证了该采样方式在Stratonovich随机积分意义下的自洽性[^ref-garcia-palacios-sllg]。

值得指出的是，在前期一些公开的VGSOT-MTJ紧凑模型的Python移植实现[^ref-zhang-vgsot]中，存在一个对统计学结论影响显著的实现细节差异：$$\boldsymbol{\xi}$$被生成为三维高斯样本后再显式归一化为单位向量$$\boldsymbol{\xi}\leftarrow\boldsymbol{\xi}/|\boldsymbol{\xi}|$$，相当于将$$|\mathbf{H}_{\mathrm{TH}}|^2$$从$$\chi^2_3$$分布锁定为常数，违反了FDT对各分量方差的硬性约束 (正确实现下$$\mathbb{E}[|\mathbf{H}_{\mathrm{TH}}|^2]=3\sigma_{\mathrm{th}}^2$$，而归一化后$$|\mathbf{H}_{\mathrm{TH}}|^2\equiv\sigma_{\mathrm{th}}^2$$)，等效将注入噪声功率压缩了三倍。其物理后果是临界电压附近的$$P_{\mathrm{sw}}(V)$$曲线在仿真中显著陡于实测，从而严重低估Sigmoid斜率的D2D展宽因子$$\eta_c$$。本工作配套仿真器采用三分量独立$$\mathcal{N}(0,1)$$采样以确保与FDT严格自洽，2.3.5节中给出的$$\eta_c$$与$$\mathcal{F}(\mathrm{CV})$$标定数值即基于修正后的实现。需要注意的是，上述采样式中的温度$$T$$在引入自热效应后将成为随时间步演化的动态状态变量，而非固定常数，具体耦合更新机制在2.2.2节与2.2.3节给出。

至此，上述各有效场分量与热噪声采样表达式共同构成了驱动LLS方程的完整有效场模型。各分量的物理参数取值依据2.2.2节各小节末尾列出的器件参数表 (2.2.2.1节末表2.2、2.2.2.2节末表2.3、2.2.2.3节末表2.4)，而离散热噪声采样中温度$$T$$与材料参数$$M_s(T)$$、$$K_i(T)$$之间的耦合反馈关系，则是2.2.2节温度效应建模的核心内容。

---

---

### 2.2.2 温度效应与器件非理想性

前述磁化动力学模型主要刻画了自由层在有效磁场、自旋力矩与热噪声共同作用下的随机演化规律。然而，对于面向工程实现的VGSOT-MTJ或SOT-MTJ紧凑模型而言，仅有理想化的磁化动力学方程仍然不足以准确描述实际器件行为。其原因在于，写入过程本身会引起明显的焦耳热积累，自由层材料参数随温度变化而漂移，隧穿输运特性同时受到偏压与温度共同调制，而器件有限尺寸又会通过退磁场和形状各向异性进一步影响等效能垒与翻转阈值。因此，在概率翻转建模中，温度效应与器件非理想性是决定概率曲线展宽、阈值漂移以及读写不对称性的关键来源。下面在统一符号体系下，对这些因素进行集中建模。

需要特别说明本节各物理量与2.2.1节LLG方程之间的耦合关系。在每一个离散时间步内，热扩散方程先以当前电学状态更新器件温度$$T$$；随后，$$T$$的变化通过2.2.2.2节的温度依赖关系同步更新材料参数$$M_s(T)$$、$$K_i(T)$$和$$\eta(T)$$；更新后的材料参数重新计算有效场各分量$$\mathbf{H}_{\mathrm{PMA}}$$、$$\mathbf{H}_{\mathrm{D}}$$及热噪声幅度；最终以新的有效场驱动Cayley变换推进磁化矢量$$\mathbf{m}$$的演化。这条$$T\to$$材料参数$$\to\mathbf{H}_{\mathrm{eff}}\to\mathbf{m}$$的逐步反馈链，是本文紧凑模型区别于零温LLG仿真的核心特征，也是2.2.3节数值求解框架必须在每步内完整执行的物理约束。图2.3以示意形式汇总了本节所涉及的四条非理想通道，即自热源与热扩散模型、温度依赖的磁性参数漂移、TMR与电输运非线性、以及退磁与形状效应，并标明了各通道对概率翻转曲线阈值漂移与斜率展宽的最终影响路径。

![sMTJ温度效应和非理想性建模总览](figs/Chapter02_local_03.png)

**图2.3** sMTJ温度效应与非理想性建模。自热源与热扩散模型为器件温度提供瞬态演化通道；温度依赖材料参数漂移经由$$M_s(T)$$、$$K_i(T)$$、$$\eta(T)$$调制能垒与驱动力矩；TMR与电输运非线性确定读出电阻对磁化方向、偏压与温度的连续映射；退磁与形状效应经由退磁因子改变有效各向异性。四条通道汇聚至概率翻转$$P_{\mathrm{sw}}$$的宏观响应，表现为阈值漂移与斜率变缓两类非理想特征。

---

#### 2.2.2.1 自热效应

在写入脉冲作用下，器件内部电流耗散会转化为热源，使MTJ实际温度偏离环境温度。其基本控制方程由热扩散方程给出[^ref-li-zhang-thermal]：

$$
C_v \frac{dT}{dt} = \lambda \nabla^2 T + Q
$$

其中$$C_v$$为单位体积热容，$$\lambda$$为热导率，$$Q$$为体热源项。对于电流驱动的磁性器件，焦耳热一般写为

$$
Q = \rho j^2
$$

其中$$\rho$$为等效电阻率，$$j$$为电流密度。对于传统STT-MTJ，热源仅来自穿过MTJ结柱的隧穿电流，可将其折算为单位体积功率密度：

$$
Q_{\mathrm{MTJ}}
=
\frac{V_{\mathrm{MTJ}}^2}{R_{\mathrm{MTJ}} A_{\mathrm{MTJ}}}
\cdot
\frac{1}{t_{\mathrm{MTJ}}}
$$

将上述MTJ热源表达式代入热扩散方程，并将结柱沿厚度方向简化为一维等效热网络，可得MTJ温度演化表达式[^ref-li-zhang-thermal]：

$$
C_v t_{\mathrm{MTJ}} \frac{dT}{dt}
=
\frac{V_{\mathrm{MTJ}}^2}{R_{\mathrm{MTJ}}}
\cdot
\frac{4}{\pi D_{\mathrm{phys}}^2}
-
\frac{\lambda_{\mathrm{MgO}}}{t_{\mathrm{MgO}}}(T-T_0)
$$

其中$$t_{\mathrm{MTJ}}$$为MTJ柱总厚度，$$D_{\mathrm{phys}}$$为物理直径，$$t_{\mathrm{MgO}}$$为MgO势垒厚度，$$T_0$$为环境温度，$$\lambda_{\mathrm{MgO}}$$为MgO薄膜热导率。等号右侧第一项表示焦耳热注入，第二项表示器件通过MgO层向外散热。此处以MgO作为主散热路径是基于其热导率显著低于CoFeB和上下金属电极的事实，MgO纳米薄膜的热导率约为4 W/(m·K)，远低于金属层的50–100 W/(m·K)量级，因此MgO热阻构成整个热路的主要瓶颈，一维近似的误差对70–90 nm量级MTJ的温升估算在可接受范围内[^ref-kim-mtj-spice]。

对于SOT-MTJ或VGSOT-MTJ，自热问题更为复杂，原因是热源不再只来自隧穿电流，重金属沟道中的SOT写入电流同样会产生显著焦耳热。SOT沟道的体热源密度为

$$
Q_{\mathrm{SOT}}
=
\frac{V_{\mathrm{SOT}}^2}{R_{\mathrm{SOT}} A_{\mathrm{SOT}}}
\cdot
\frac{1}{t_{\mathrm{SOT}}}
$$

其中$$A_{\mathrm{SOT}} = W_{\mathrm{SOT}} \times L_{\mathrm{SOT}}$$为SOT沟道的俯视面积，$$t_{\mathrm{SOT}}$$为重金属层厚度。将SOT沟道热源的贡献等效折算到MTJ柱所覆盖的面积内，总的温度演化方程修正为

$$
C_v t_{\mathrm{MTJ}} \frac{dT}{dt}
=
\frac{V_{\mathrm{MTJ}}^2}{R_{\mathrm{MTJ}}}
\cdot
\frac{4}{\pi D_{\mathrm{phys}}^2}
+
\frac{V_{\mathrm{SOT}}^2}{R_{\mathrm{SOT}} A_{\mathrm{SOT}}}
\cdot
\frac{t_{\mathrm{MTJ}}}{L_{\mathrm{SOT}}}
-
\frac{\lambda_{\mathrm{MgO}}}{t_{\mathrm{MgO}}}(T-T_0)
$$

需要指出，上述总温度演化方程对散热路径沿用了纯STT情形下以MgO为主的单路径假设。对于SOT-MTJ结构，MTJ与重金属沟道之间的底部界面同样提供热传导通道；然而考虑到该界面处存在TaN或Ta种子层，其等效界面热阻与MgO层量级相近，因此以MgO单路径表征整体热阻仍为合理的工程近似，更精确的处理可引入双路径热网络，对此不在本文讨论范围。

该模型说明，SOT写入过程中器件温升往往高于纯STT情形，原因是沟道电流密度通常更高，且额外热功率直接注入MTJ下方。从仿真实现角度，上述温度演化方程能够以较低计算代价近似取代三维热有限元分析，直接嵌入行为级平台；而其动态性，即$$T$$随每个时间步内$$V_{\mathrm{MTJ}}$$和$$V_{\mathrm{SOT}}$$的变化而更新，正是将自热效应纳入概率翻转建模不可回避的物理前提。

为定量考察上述热模型在本文基准器件上的行为，将器件参数代入温度演化方程，得到如图2.4所示的瞬态温升仿真结果。器件热时间常数$$\tau_{\mathrm{th}} = C_v t_{\mathrm{MTJ}} t_{\mathrm{MgO}}/\lambda_{\mathrm{MgO}}$$在本参数下约为17.5 ps，显著小于典型纳秒级写入脉冲宽度，这意味着在整个写入窗口内器件已充分进入热稳态，脉冲中后段的行为可由稳态温升刻画。图2.4(a)对比了三种写入方式在相同0.8 V驱动下的稳态温升：纯STT写入由于隧穿电流受隧道电阻限制，温升仅约6 K；纯SOT写入由于沟道电阻较低，稳态温升达到28 K；STT与SOT同时施加时温升约为34 K，近似等于两种独立热源稳态温升之和，表明在本文采用的等效一维热路模型中两条焦耳热通道满足线性叠加关系。图2.4(b)进一步给出纯SOT模式下$$V_{\mathrm{SOT}}$$从0.4 V扫至1.2 V的瞬态温升曲线，稳态温升依次为7 K、16 K、28 K、44 K与63 K，其随驱动电压的平方律标度关系与沟道焦耳功率$$P_{\mathrm{SOT}}=V_{\mathrm{SOT}}^2/R_{\mathrm{SOT}}$$的理论预期吻合。该标度律意味着高压驱动下自热与SOT力矩同时被放大，自由层所感受到的热噪声幅度按$$\sqrt{T}$$比例增长，因而在临界电流附近相同幅值电流所对应的翻转概率分布会随驱动强度而系统性展宽。

![瞬态热演化仿真结果](figs/Chapter02_local_04.png)

**图2.4** 基于本文参数的自热瞬态仿真。(a)三种写入模式在0.8 V驱动下的器件温度演化，热时间常数$$\tau_{\mathrm{th}}\approx17.5$$ ps，STT与SOT热源稳态温升近似线性叠加。(b)纯SOT模式下$$V_{\mathrm{SOT}}$$从0.4 V扫至1.2 V的温度演化，稳态温升随驱动电压平方增长，印证了焦耳功率标度关系。

上述自热模型涉及的SOT沟道几何与电学参数、热容与MgO热导率、以及VCMA系数集中列于表2.2。沟道几何与电学参数用于计算沟道焦耳热源项及2.2.1.1节给出的等效DL场幅值；热容$$C_v$$与MgO热导率$$\lambda_{\mathrm{MgO}}$$直接决定温度演化方程中的温升响应速率，取自经典MTJ热传导仿真基准；VCMA系数$$\beta_{\mathrm{VCMA}}$$依据CoFeB/MgO界面的典型实验测量范围选取。

**表2.2** SOT沟道输运、热学与VCMA文献参数。

| **参数符号** | **物理意义** | **设定值** | **说明** |
|---|---|---|---|
| $$P$$ | 自旋极化率 | $$0.58$$ | MTJ隧穿自旋力矩驱动效率 |
| $$\beta_{\mathrm{VCMA}}$$ | VCMA系数 | $$60\,\mathrm{fJ/(V\cdot m)}$$ | 电压对界面各向异性能的调制强度 |
| $$L_{\mathrm{SOT}}$$ | 沟道长度 | $$240\,\mathrm{nm}$$ | SOT重金属写入沟道物理长度 |
| $$W_{\mathrm{SOT}}$$ | 沟道宽度 | $$200\,\mathrm{nm}$$ | SOT重金属写入沟道物理宽度 |
| $$T_{\mathrm{SOT}}$$ | 沟道厚度 | $$4.3\,\mathrm{nm}$$ | 决定SOT电流密度分布 |
| $$\rho_{\mathrm{SOT}}$$ | 沟道电阻率 | $$2.78\times10^{-6}\,\Omega\cdot\mathrm{m}$$ | Ta/W重金属层典型电阻率 |
| $$C_v$$ | 等效体热容 | $$2.5\times10^6\,\mathrm{J/(m^3\cdot K)}$$ | 暂态温升响应$$dT/dt$$计算依据 |
| $$\lambda_{\mathrm{MgO}}$$ | MgO热导率 | $$4\,\mathrm{W/(m\cdot K)}$$ | 纳米薄膜值显著低于块体，是自热主因 |
| $$\phi$$ | MgO有效势垒高度 | $$0.4\,\mathrm{eV}$$ | 决定基础平行态电阻$$R_P$$的能垒 |

---

#### 2.2.2.2 温度依赖磁性参数

温度升高不仅通过热噪声改变翻转统计，还会直接改变自由层的内禀材料参数。为保持模型闭环，须将饱和磁化强度、界面各向异性能密度和自旋极化率均写成温度的函数，以使每一时间步内有效场的计算反映真实的材料状态。

**饱和磁化强度$$M_s(T)$$。** 铁磁体的自发磁化随温度升高而衰减，在远低于居里温度$$T_C$$的工作区间内，由低能自旋波 (magnon) 激发主导的衰减遵从Bloch $$T^{3/2}$$定律[^ref-bloch-law]。对于CoFeB超薄膜，以室温 (RT) 实测值$$M_s(\mathrm{RT})$$作为归一化基准，可写出适合紧凑模型参数提取的表达式：
$$
M_s(T)
=
M_s(\mathrm{RT})
\frac{1-(T/T_C)^{3/2}}
{1-(\mathrm{RT}/T_C)^{3/2}}
$$

该归一化形式避免了对$$M_s(0\,\mathrm{K})$$的直接测量依赖，可直接与振动样品磁强计 (VSM) 或超导量子干涉仪 (SQUID) 的室温测量结果对接。随着$$T$$升高，$$M_s(T)$$单调下降，从而使退磁场、各向异性场以及零温临界驱动电流均发生漂移[^ref-dieny-pma-review]。

**自旋极化率$$\eta(T)$$。** 在Julliere隧穿模型框架下[^ref-julliere-tmr]，MTJ的自旋极化率$$\eta$$正比于费米能级处的自旋劈裂密度之差，而在平均场近似中，该量跟随磁化强度一同衰减。因此，有效自旋极化率遵从与$$M_s(T)$$相同的温度依赖形式：

$$
\eta(T)
=
\eta(\mathrm{RT})
\frac{1-(T/T_C)^{3/2}}
{1-(\mathrm{RT}/T_C)^{3/2}}
$$

温度上升时自旋极化率下降，导致STT及SOT的等效驱动力均减弱。对于本文以SOT为主驱动的VGSOT架构，该项可被合并入等效自旋霍尔注入效率的温度修正，其物理含义是：在较高温度下，相同电流密度所能提供的有效翻转力矩有所削减。

**界面各向异性能密度$$K_i(T)$$。** Callen–Callen理论[^ref-callen-callen]预测，对于量子数$$l=2$$的单轴各向异性，其温度衰减指数在三维铁磁体中为10/3；然而在CoFeB/MgO超薄界面体系中，二维费米面的修正以及界面态对各向异性能的主导贡献使实验观测到的指数偏离该理论值，通常落在2至3之间[^ref-dieny-pma-review]。结合对CoFeB/MgO薄膜的温变铁磁共振 (FMR) 测量数据的拟合，取指数为2.18，$$K_i(T)$$的表达式为：

$$
K_i(T)
=
\frac{1}{2} \mu_0 t_{\mathrm{FL}} M_s(\mathrm{RT})
\left[
M_s(\mathrm{RT}) + H_k^{\mathrm{eff}}(\mathrm{RT})
\right]
\left(
\frac{1-(T/T_C)^{3/2}}
{1-(\mathrm{RT}/T_C)^{3/2}}
\right)^{2.18}
$$

其中$$t_{\mathrm{FL}}$$为自由层厚度，$$H_k^{\mathrm{eff}}(\mathrm{RT})$$为室温下由FMR测得的有效各向异性场。本式中$$M_s$$与$$H_k^{\mathrm{eff}}$$均以SI单位$$\mathrm{A/m}$$表示，$$\mu_0$$作为整体因子置于括号外，可保证括号内同量纲相加，避免因$$\mu_0 M_s$$ (单位$$\mathrm{T}$$) 与$$H_k^{\mathrm{eff}}$$ (单位$$\mathrm{A/m}$$) 混用而造成的量纲不自洽。上述$$K_i(T)$$表达式中的外层因子$$\frac{1}{2}\mu_0 t_{\mathrm{FL}}M_s(\mathrm{RT})[M_s(\mathrm{RT})+H_k^{\mathrm{eff}}(\mathrm{RT})]$$即为室温界面各向异性能密度的校准值 (单位$$\mathrm{J/m^2}$$)，随后的幂律项描述其温度衰减。相应的各向异性场$$z$$分量为

$$
H_{\mathrm{ani},z}
=
\frac{2K_i(T)}{t_{\mathrm{FL}} M_s(T)} m_z
$$

$$H_{\mathrm{ani},x}=H_{\mathrm{ani},y}=0$$。由于$$K_i(T)$$的衰减指数2.18大于$$M_s(T)$$的1.5，各向异性场$$H_{\mathrm{ani},z}$$随温度的下降幅度超过$$M_s$$本身，等效能垒的减小不能简单等同于$$M_s$$的衰减。

在250 K至500 K区间对$$M_s(T)$$、$$K_i(T)$$、$$\eta(T)$$进行数值求值，可得图2.5(a)–(c)所示的参数漂移曲线：$$M_s$$从250 K时的650 kA/m下降至500 K时的505 kA/m，相对衰减约22% ($$\eta$$同步) ；$$K_i$$由0.349 mJ/m²降至0.202 mJ/m²，相对衰减约42%，显著高于$$M_s$$的衰减幅度。

![温度依赖材料参数与TMR偏压响应](figs/Chapter02_local_05.png)

**图2.5** 基于本文参数集的温度依赖材料参数与电输运非线性仿真。(a)饱和磁化$$M_s(T)$$遵从Bloch自旋波标度；(b)界面各向异性$$K_i(T)$$以修正Callen–Callen指数2.18衰减，衰减速率显著高于$$M_s$$；(c)自旋极化率$$\eta(T)$$与$$M_s$$同步；(d)TMR随$$V_{\mathrm{MTJ}}$$的衰减对比：实线为本文采纳的三参数二次-有理形式，虚线为Zhang等人采用的单参数Lorentzian形式 ($$V_h=0.5$$ V)，两者在约1.2 V附近交叉。RT(300 K)为(a)–(c)的室温基准点。

温度依赖关系所归一化的室温基准值$$M_s(\mathrm{RT})$$、$$K_i(0)$$、$$\eta(\mathrm{RT})$$，以及沿用一致的磁学参数$$\alpha$$、$$\gamma$$、$$T_C$$、$$A_{\mathrm{ex}}$$与物理尺寸$$D_{\mathrm{phys}}$$、$$D_{\mathrm{elec}}$$、$$t_f$$、$$t_{\mathrm{ox}}$$、$$t_{\mathrm{MTJ}}$$共同列于表2.3。参数量级以CoFeB/MgO垂直各向异性MTJ的公开材料范围为约束[^ref-dieny-pma-review]，交换偏置场的设计参考VCMA-SOT联合翻转机制的场无关翻转模型[^ref-li-jiang-vcma-sot]；居里温度$$T_C = 1100\,\mathrm{K}$$作为Bloch定律与Callen-Callen定律的温度标度基准，并经本节幂律拟合验证。

**表2.3** 80 nm级SOT-MTJ核心磁学与物理尺寸参数。

| **参数符号** | **物理意义** | **设定值** | **说明** |
|---|---|---|---|
| $$D_{\mathrm{phys}}$$ | MTJ物理直径 | $$80\,\mathrm{nm}$$ | 基础阵列物理节点尺寸 |
| $$D_{\mathrm{elec}}$$ | MTJ电学有效直径 | $$\approx 65\,\mathrm{nm}$$ | 计入边缘刻蚀损伤带 ($$\delta_{\mathrm{edge}}\!\approx\!7.5\,\mathrm{nm}$$)，与表2.5实测$$R_P$$自洽 |
| $$t_f$$ | 自由层厚度 | $$1.1\,\mathrm{nm}$$ | 影响PMA能与热稳定体积 |
| $$t_{\mathrm{ox}}$$ | 隧道势垒厚度 | $$1.4\,\mathrm{nm}$$ | 决定隧穿电阻量级 |
| $$t_{\mathrm{MTJ}}$$ | MTJ柱总厚度 | $$\approx 20\,\mathrm{nm}$$ | 用于热扩散方程的等效一维热路 (与$$\tau_{\mathrm{th}}=C_v t_{\mathrm{MTJ}}t_{\mathrm{MgO}}/\lambda_{\mathrm{MgO}}\!\approx\!17.5\,\mathrm{ps}$$对应) |
| $$M_s$$ | 室温饱和磁化强度 | $$6.25\times10^5\,\mathrm{A/m}$$ | 典型超薄CoFeB实验测量校准值 |
| $$K_i(0)$$ | 界面PMA强度 | $$3.2\times10^{-4}\,\mathrm{J/m^2}$$ | 零偏压本征界面各向异性能密度 |
| $$\alpha$$ | Gilbert阻尼系数 | $$0.05$$ | 界面改性CoFeB薄膜典型经验值 |
| $$\gamma$$ | 旋磁比 | $$2.2127\times10^5\,\mathrm{m/(A\cdot s)}$$ | 磁矩进动频率基础物理常数 |
| $$T_C$$ | 居里温度 | $$1100\,\mathrm{K}$$ | 磁学参数温度衰减速率的拟合基准 |
| $$A_{\mathrm{ex}}(\mathrm{RT})$$ | 室温交换刚度常数 | $$4\,\mathrm{pJ/m}$$ | 畴壁形成能与热稳定性评估依据 |

---

#### 2.2.2.3 TMR与电输运非线性

除磁性参数外，电输运特性本身也具有明显的偏压和温度依赖性。对MTJ而言，这主要体现在TMR的非线性衰减以及平行态电阻$$R_P$$对势垒参数的指数敏感性上。由于读出电路最终观测到的是电阻或电流而非直接的磁化矢量，输运非线性会进一步影响写后判决边界、参考电路设计以及阵列级误码率。

在统一模型中，若自由层$$z$$向磁化分量记为$$m_z$$，参考层磁化方向取$$-z$$，则MTJ阻值可由两极限态电导的线性内插得到：

$$
G_{\mathrm{MTJ}}(m_z) = \frac{1-m_z}{2}\,G_P + \frac{1+m_z}{2}\,G_{AP}
$$

以$$k=\mathrm{TMR}/(\mathrm{TMR}+2)$$为内插系数，上式可等价地改写为常用的有理形式：

$$
R_{\mathrm{MTJ}}(m_z)
=
R_P\,\frac{1+k}{1-k\,m_z}
=
\frac{R_P\,R_{AP}}{\dfrac{1+m_z}{2}R_P + \dfrac{1-m_z}{2}R_{AP}}
$$

直接代入极限值可验证该表达式的物理自洽性：$$m_z=-1$$给出$$R_{\mathrm{MTJ}}=R_P$$ (平行态)，$$m_z=+1$$给出$$R_{\mathrm{MTJ}}=R_P(1+\mathrm{TMR})=R_{AP}$$ (反平行态)。磁化状态与电阻之间的连续映射在随机翻转轨迹尚未完全收敛时尤为重要，因为中间磁化态对应的瞬时电阻既不等于$$R_P$$也不等于$$R_{AP}$$。

TMR受偏压影响而衰减，这一现象源于较大偏压下界面处感应态散射增强、magnon激发以及隧穿相干性减弱[^ref-akerman-tmr]。在紧凑建模文献中针对TMR偏压依赖形成了两类经验形式。其一是单参数Lorentzian形式：

$$
\mathrm{TMR}_{\mathrm{Lor}}(V_{\mathrm{MTJ}})
=
\frac{\mathrm{TMR}_0}{1+V_{\mathrm{MTJ}}^2/V_h^2}
$$

其唯一拟合参数$$V_h$$定义为TMR衰减到零偏压值一半时对应的偏压，典型取值$$V_h=0.5\,\mathrm{V}$$。该式由偏压下磁激发与非弹性散射所引起的TMR衰减唯象推导而来，其中$$V^2$$项对应偏压引起的磁激发功率线性项的积分。其二是本文采纳的、由Hikstor SOT-MRAM工艺PDK提取的三参数二次-有理形式：

$$
\mathrm{TMR}_{\mathrm{PDK}}(V_{\mathrm{MTJ}})
=
\frac{\mathrm{TMR}_0}{k_{\mathrm{TMR}}}
\left[
\frac{1}{a_{\mathrm{TMR}} V_{\mathrm{MTJ}}^2 + b_{\mathrm{TMR}} |V_{\mathrm{MTJ}}| + c_{\mathrm{TMR}}} - 1
\right]
$$

两种模型的定量对比示于图2.5(d)：Lorentzian曲线遵从严格的$$V^2$$对称性，在$$V_h=0.5$$ V处即快速衰减到$$\mathrm{TMR}_0/2$$附近，此后以长拖尾方式渐近趋于零；PDK曲线因含$$|V|$$线性项，在$$0.25$$至$$1.0$$ V的中等偏压区间保留较高TMR (较Lorentzian高约8–17个百分点)，但在1.2 V以上因二次项主导迅速跌至零，两者曲线约在1.2 V附近交叉。两种模型的取舍反映了参数经济性与测量吻合度的权衡：Lorentzian形式仅含单个物理可解释的参数$$V_h$$，对于缺乏详细测量数据的早期器件建模非常方便；PDK形式的三参数$$|V|$$线性项可吸收势垒不对称、反铁磁钉扎界面贡献等实际非理想因素，在具有完整制程数据的工程化模型中拟合精度更高，但参数物理意义不如$$V_h$$直接。本文选用PDK形式的主要原因是其具有明确的实验测量基础，且其在关键写入偏压窗口 ($$0.5$$ V至$$1.0$$ V) 内的精度对翻转能效与读出裕量预估至关重要；在缺少PDK数据的外延器件仿真中可替换为Lorentzian形式，以保持同一仿真接口下的模型可迁移性。

平行态电阻$$R_P$$依赖于MgO势垒的厚度与高度，其物理图像来自Brinkman–Dynes–Rowell (BDR) 隧穿模型[^ref-brinkman-bdr]。在WKB近似下，对抛物线形势垒的零偏压隧穿电导$$G_P\propto\sqrt{\phi_{\mathrm{ox}}}\exp(-2t_{\mathrm{ox}}\sqrt{2m_e e\phi_{\mathrm{ox}}}/\hbar)$$，反演得平行态电阻为

$$
R_P
=
\frac{t_{\mathrm{ox}}}{F\,A_{\mathrm{MTJ}}\sqrt{\phi_{\mathrm{ox}}}}
\exp\!\left(
\frac{2\,t_{\mathrm{ox}}\sqrt{2m_e e\phi_{\mathrm{ox}}}}{\hbar}
\right)
$$

其中$$t_{\mathrm{ox}}$$为势垒厚度，$$\phi_{\mathrm{ox}}$$为MgO有效势垒高度，$$A_{\mathrm{MTJ}}$$为电学有效面积 (按2.2.2.4节取$$D_{\mathrm{elec}}$$对应面积)，$$m_e$$为电子质量，$$e$$为元电荷，$$\hbar$$为约化普朗克常数，$$F$$为由R·A乘积一致性约束所定标的常数 (其量纲为$$1/(\Omega\!\cdot\!\mathrm{m}\!\cdot\!\sqrt{\mathrm{eV}})$$，并非无量纲量；数值上保证了$$R_P\cdot A_{\mathrm{MTJ}}$$回归到实验测得的R·A值)。BDR模型严格成立于$$eV\ll\phi_{\mathrm{ox}}$$的弱偏压极限；在较大偏压下，该模型作为势垒参数到阻值的定性映射仍具工程适用性，但应理解为等效参数化而非严格推导。上述表达式揭示了$$R_P$$对$$t_{\mathrm{ox}}$$与$$\phi_{\mathrm{ox}}$$的指数敏感性，这意味着工艺波动中的势垒厚度涨落会被指数放大为阻值分布，进而通过前述$$R_{\mathrm{MTJ}}(m_z)$$映射关系影响整个阵列的$$R_P/R_{AP}$$离散性。

在概率计算场景下，TMR与电输运非线性的影响体现在以下几个层面。读出端的电阻窗口会随偏压与温度实时变化，从而影响读出参考电压与感放裕量。由于写入过程中自热升温，写后立刻读取与热平衡后读取对应不同的瞬时TMR，导致动态读出误差。在阵列级环境中，阻值分布非线性会在多次统计采样中引入额外的均值偏置。因此，行为级模型中有必要将$$R_{\mathrm{MTJ}}(m_z,T,V)$$而非固定的$$R_P$$/$$R_{AP}$$作为读出接口变量。

上述TMR偏压衰减拟合式中的二次-有理系数直接取自Hikstor SOT-MRAM工艺PDK Verilog-A模型的参数提取结果，可精确重现实测偏压下TMR的非线性衰减行为；TMR$$_0$$、$$R\!\cdot\!A$$、$$\theta_{\mathrm{SH}}$$三项由2.3.2节实验阈值$$(R_P,R_{AP},V_{\mathrm{th}})$$联合反推得到。表2.4集中列出本文采纳的TMR偏压衰减PDK拟合系数与端口级标定参数。

**表2.4** TMR偏压衰减PDK拟合系数与端口级标定参数。

| **参数符号** | **物理意义** | **设定值** | **说明** |
|---|---|---|---|
| $$\mathrm{TMR}_0$$ | 零偏压TMR比值 | $$1.00$$ | 标定至2.3.2节 滞回回线幅度$$R_{AP}/R_P\approx 2$$ |
| $$R\!\cdot\!A$$ (标定值) | 电阻面积积 | $$16.6\,\Omega\!\cdot\!\mu\mathrm{m}^2$$ | 与$$D_{\mathrm{elec}}=65\,\mathrm{nm}$$配合给出$$R_P\approx 5\,\mathrm{k}\Omega$$ |
| $$\theta_{\mathrm{SH}}$$ (标定值) | 有效自旋霍尔角 | $$0.04$$ | 使$$V_{\mathrm{th}}^{\mathrm{sim}}(0.75\,\mathrm{ns})\approx V_{\mathrm{th}}^{\mathrm{exp}}=894\,\mathrm{mV}$$ |
| $$k_{\mathrm{TMR}}$$ | TMR衰减归一化项 | $$1.2346$$ | PDK偏压依赖方程归一化系数 |
| $$a_{\mathrm{TMR}}$$ | 二次项系数 | $$0.1729$$ | $$V_{\mathrm{mtj}}^2$$对TMR的衰减权重 |
| $$b_{\mathrm{TMR}}$$ | 一次项系数 | $$0.1315$$ | $$|V_{\mathrm{mtj}}|$$对TMR的衰减权重 |
| $$c_{\mathrm{TMR}}$$ | 常数项偏移 | $$0.4475$$ | 拟合方程基准常数 |

---

#### 2.2.2.4 退磁与形状效应

有限尺寸MTJ的形状效应主要通过退磁因子进入有效场表达式。虽然宏自旋模型将自由层视为单一磁矩，但器件几何尺寸仍会通过退磁场影响垂直各向异性补偿关系，从而改变能垒高度和热稳定性。2.2.1.2节给出了适用于$$t_f \ll D_{\mathrm{phys}}$$薄圆盘极限的退磁因子线性近似，本节则给出应用于实际仿真的扁椭球体精确解析表达式，后者是计算后续有效各向异性能$$K_U^{\mathrm{eff}}(T)$$的正确基础。

对于圆形MTJ柱，自由层近似为扁椭球体，其横向退磁因子$$N_x = N_y$$满足解析表达式[^ref-stoner-wohlfarth]：

$$
N_x
=
\frac{q^2}{2(1-q^2)}
\left[
\frac{1}{q\sqrt{1-q^2}}\arccos(q)-1
\right]
$$

其中宽厚比参数定义为

$$
q = \frac{t_{\mathrm{FL}}}{D_{\mathrm{elec}}}
$$

$$D_{\mathrm{elec}}$$为电学有效直径，因边缘刻蚀工艺 (通常为离子铣刻) 导致自由层边缘受损区域导电性降低，$$D_{\mathrm{elec}}$$通常比物理直径$$D_{\mathrm{phys}}$$小约5–10 nm。由圆柱对称性，有

$$
N_x=N_y,\qquad N_z=1-2N_x
$$

退磁场各分量则为

$$
H_{\mathrm{demag},x} = -\mu_0 M_s(T) N_x m_x
$$

$$
H_{\mathrm{demag},y} = -\mu_0 M_s(T) N_x m_y
$$

$$
H_{\mathrm{demag},z} = -\mu_0 M_s(T) (1-2N_x) m_z
$$

形状效应通过$$N_z - N_x$$直接进入有效各向异性能密度，结合温度依赖参数，净各向异性能为

$$
K_U^{\mathrm{eff}}(T)
=
K_i(T)
-
\frac{1}{2}\mu_0 M_s^2(T)(N_z-N_x)
$$

由此可见，随着$$T$$升高，$$M_s(T)$$和$$K_i(T)$$均下降，但两者以不同的幂律速率衰减，从而使$$K_U^{\mathrm{eff}}(T)$$的变化不可简化为单一参数的线性插值，而必须逐步在耦合循环中更新。圆形器件$$N_x = N_y$$的对称性还确保了横向退磁场在统计意义上各向同性，不引入额外的翻转方向偏好，这也是圆形MTJ在概率器件建模中优于椭圆或不规则形状的重要理由。

---

#### 2.2.2.5 自热-材料耦合的仿真验证

为定量考察前述四条非理想通道 (自热、温度依赖材料参数、退磁场与TMR非线性) 在LLG时间步内的耦合反馈，本工作在仿真器中接通完整的温度到材料参数再到有效场的反馈链：每个时间步内先按2.2.2.1节给出的一维RC热扩散方程更新瞬时温度$$T(t)$$，再按2.2.2.2节的Bloch和修正Callen–Callen标度计算$$M_s(T(t))$$与$$K_i(T(t))$$，并将温度修正后的材料参数代入有效场各分量。在自热关闭的对照仿真中$$T\equiv300\,\mathrm{K}$$恒定，所有其他设置 (脉冲波形、初始角度抽样种子、数值积分步长) 严格一致，因此两条轨迹之间的差异完全由温度反馈引入，无须额外的噪声平均即可定量分离自热效应的贡献。

![自热反馈对纯SOT翻转轨迹的影响](figs/Chapter02_local_06.png)

**图2.6** 纯SOT写入操作点 ($$V_{\mathrm{MTJ}}=0\,\mathrm{V}$$、$$I_{\mathrm{SOT}}=-1500\,\mu\mathrm{A}$$，对应沟道电压$$V_{\mathrm{SOT}}\approx 1.16\,\mathrm{V}$$、3 ns写入脉冲 + 5 ns弛豫) 下自热反馈对磁化轨迹的影响。(a)$$m_z(t)$$对比：自热关闭 (实线，黑) 与自热开启 (虚线，红) ；脉冲在$$t=3\,\mathrm{ns}$$处关断后两条轨迹均完成$$-1\!\to\!+1$$方向翻转；插图显示二者之差$$m_z^{\mathrm{ON}}-m_z^{\mathrm{OFF}}$$ (放大100×) 在脉冲前段最大达6%–8%、弛豫阶段振荡幅度达10%，是自热引致进动周期变化的可视化。(b)自热开启情形下$$T(t)$$的时变轨迹：脉冲启动后$$T(t)$$以$$\tau_{\mathrm{th}}\approx17.5\,\mathrm{ps}$$的指数刚性上升至闭式解预测的稳态$$T_{\mathrm{eq}}=359.2\,\mathrm{K}$$ (琥珀色虚线，$$\Delta T_{\mathrm{eq}}\approx 59\,\mathrm{K}$$)，脉冲关断后以同样的时间常数指数衰减回环境温度。(c)$$R_{\mathrm{MTJ}}(t)$$从$$R_{AP}$$跃迁至$$R_P$$；两条曲线在过渡区附近的细微相位差来源于自热开启情形下$$H_{\mathrm{PMA}}$$的轻微削弱使进动频率略有降低。(d)自由层材料参数的相对漂移：峰值温度对应的$$\Delta M_s/M_s\approx-5.15\%$$、$$\Delta K_i/K_i\approx-10.89\%$$，按Callen–Callen指数2.18与Bloch指数1.5的差异$$K_i$$比$$M_s$$衰减更快，最终有效PMA场漂移$$\Delta H_{\mathrm{PMA}}/H_{\mathrm{PMA}}\approx-5.74\%$$。

仿真给出两点对后续概率建模具有方法论意义的结论：其一，$$\tau_{\mathrm{th}}\!\sim\!17.5\,\mathrm{ps}$$远小于纳秒级写入脉冲宽度，MgO热阻主导的一维热路在任何实际脉冲内部都迅速进入稳态，故后续概率仿真可使用稳态温升对各操作点进行单点温度修正，无须维持完整的瞬态$$T(t)$$；其二，即便在深超阈值工作点上$$\Delta K_i/K_i\approx-11\%$$，$$m_z(t)$$主轨迹差异仍仅在$$10^{-2}$$量级，表明单器件级自热反馈不构成概率曲线展宽的主导机制，相关讨论在2.3.5节展开。

---

### 2.2.3 数值求解方法

在前述章节中，磁化动力学已经通过包含自旋轨道力矩、VCMA调制以及热噪声项的扩展Landau–Lifshitz–Gilbert (sLLG) 方程给出，并在2.2.2节中进一步引入了温度依赖材料参数与自热反馈。由于sLLG方程属于典型的强非线性随机微分方程 (SDE)，且在整个时间演化过程中须严格满足磁化矢量模长守恒 ($$|\mathbf{m}|=1$$)，同时温度状态变量在每一步内与磁化状态共同更新，传统的显式数值积分方法 (如显式Euler或常规Runge-Kutta) 往往面临严重的数值漂移和稳定性问题。本节推导适合大规模Monte Carlo仿真的保模长几何积分算法，并说明其在Stratonovich随机微积分框架下的收敛性质。

---

#### 2.2.3.1 sLLG方程的统一形式与刚性问题

为与前文物理模型严格对齐，热噪声作为等效磁场的一部分而非独立力矩项参与进动与阻尼过程。将包含热噪声的总有效场记为$$\mathbf{H}_{\mathrm{tot}} = \mathbf{H}_{\mathrm{eff}} + \mathbf{H}_{\mathrm{TH}}$$，其中$$\mathbf{H}_{\mathrm{TH}}$$在每个离散时间步$$\Delta t$$内按2.2.1.2节给出的离散采样关系重新采样，采样幅度依赖于当前时间步的器件温度$$T$$。结合前述LLS显式形式，sLLG方程统一表示为

$$
\frac{d\mathbf{m}}{dt}
= -\frac{\gamma}{1+\alpha^2}
\left[
\mathbf{m} \times \mathbf{H}_{\mathrm{tot}}
+ \alpha\,\mathbf{m} \times (\mathbf{m} \times \mathbf{H}_{\mathrm{tot}})
\right]
+ \boldsymbol{\tau}_{\mathrm{SOT}}
$$

其中SOT力矩项$$\boldsymbol{\tau}_{\mathrm{SOT}}$$由2.2.1.1节给出的分解式给出。利用矢量恒等式$$\mathbf{m} \times (\mathbf{m} \times \mathbf{H}) = (\mathbf{m} \cdot \mathbf{H})\mathbf{m} - \mathbf{H}$$ (基于$$|\mathbf{m}|=1$$的前提)，上述方程可改写为广义旋转动力学形式：

$$
\frac{d\mathbf{m}}{dt} = \mathbf{w}(\mathbf{m}, t) \times \mathbf{m}
$$

其中$$\mathbf{w}(\mathbf{m}, t)$$是等效旋转角速度矢量，综合了有效场、阻尼以及全部自旋力矩的贡献。

由于LLG方程中进动项$$\mathbf{m} \times \mathbf{H}_{\mathrm{tot}}$$的特征频率$$\sim\gamma H_{\mathrm{tot}}$$ (量级为GHz–THz) 远大于阻尼弛豫速率$$\sim\alpha\gamma H_{\mathrm{tot}}$$ ($$\alpha \ll 1$$)，系统具有显著的刚性 (stiffness) [^ref-ascher-petzold]。对刚性常微分方程而言，显式Euler步$$\mathbf{m}_{n+1} = \mathbf{m}_n + \Delta t\,d\mathbf{m}_n/dt$$要求$$\Delta t < 2/(\gamma H_{\mathrm{tot}})$$以保证稳定，而这一条件在典型PMA器件的进动频率下意味着皮秒量级的步长限制。更严重的是，显式步沿切线方向推进，新状态$$\mathbf{m}_{n+1}$$必然偏离单位球面 ($$|\mathbf{m}_{n+1}| > 1$$)。若通过手动归一化$$\mathbf{m} \leftarrow \mathbf{m}/|\mathbf{m}|$$强制修正，不仅在长时间积分中引入截断误差的累积，还会改变Stratonovich随机积分中热噪声的实际统计权重，导致仿真温度偏离物理设定值。

在实际仿真循环中，每个时间步$$n$$的完整执行顺序如下：首先以第$$n$$步的电学状态 ($$V_{\mathrm{MTJ}}$$，$$V_{\mathrm{SOT}}$$) 驱动2.2.2.1节给出的温度演化方程更新器件温度$$T_{n+1}$$；然后以$$T_{n+1}$$代入2.2.2.2节的温度依赖关系更新材料参数$$M_s$$、$$K_i$$、$$\eta$$；再以更新后的参数重新计算有效场各分量$$\mathbf{H}_{\mathrm{PMA}}$$、$$\mathbf{H}_{\mathrm{VCMA}}$$、$$\mathbf{H}_{\mathrm{D}}$$并重新采样热噪声场$$\mathbf{H}_{\mathrm{TH}}$$；最后以完整的$$\mathbf{H}_{\mathrm{tot}}$$执行2.2.3.2节所述的保模长Cayley步，将$$\mathbf{m}_n$$推进至$$\mathbf{m}_{n+1}$$。这一顺序确保了热扩散方程与磁化动力学在同一时间步内的物理自洽，是2.2.2节耦合反馈机制在数值实现层面的直接对应。

---

#### 2.2.3.2 保模长的隐式中点法与Cayley变换

为避免频繁归一化带来的误差，一种直观的替代方案是将$$\mathbf{m}$$转换至球坐标系$$(\theta,\phi)$$下求解。然而，由于PMA器件的稳定态位于$$m_z \approx \pm 1$$ (即极点$$\theta = 0,\pi$$附近)，球坐标系下的运动方程包含$$1/\sin\theta$$的坐标奇点，在极点附近引发数值发散，因此球坐标更新法不适用于PMA-MRAM的可靠性仿真[^note-dev-cayley]。

更为严谨的方案是在笛卡尔坐标系下采用保结构的几何积分器[^ref-weinan-wang]。对上述广义旋转动力学形式，在时间区间$$[t_n, t_{n+1}]$$内以中点处的状态近似旋转矢量：

$$
\frac{\mathbf{m}_{n+1} - \mathbf{m}_n}{\Delta t}
= \mathbf{w}_n \times \left(\frac{\mathbf{m}_{n+1} + \mathbf{m}_n}{2}\right)
$$

其中$$\mathbf{w}_n$$用当前步状态$$\mathbf{m}_n$$和已知外加场评估。将上式展开并移项，分离已知量与未知量：

$$
\mathbf{m}_{n+1} - \frac{\Delta t}{2}(\mathbf{w}_n \times \mathbf{m}_{n+1})
= \mathbf{m}_n + \frac{\Delta t}{2}(\mathbf{w}_n \times \mathbf{m}_n)
$$

设$$\mathbf{w}_n = (w_x, w_y, w_z)^{\mathrm{T}}$$，引入叉乘的反对称矩阵表示：

$$
\mathbf{\Omega}_n
= \begin{pmatrix}
0 & -w_z & w_y \\
w_z & 0 & -w_x \\
-w_y & w_x & 0
\end{pmatrix}
$$

则上式化为线性代数方程组：

$$
\left(\mathbf{I} - \frac{\Delta t}{2}\mathbf{\Omega}_n\right)\mathbf{m}_{n+1}
= \left(\mathbf{I} + \frac{\Delta t}{2}\mathbf{\Omega}_n\right)\mathbf{m}_n
$$

解析求解得到更新公式：

$$
\mathbf{m}_{n+1}
= \underbrace{\left(\mathbf{I} - \frac{\Delta t}{2}\mathbf{\Omega}_n\right)^{-1}
\left(\mathbf{I} + \frac{\Delta t}{2}\mathbf{\Omega}_n\right)}_{\displaystyle\mathbf{A}_n}
\mathbf{m}_n
$$

矩阵$$\mathbf{A}_n$$称为**Cayley变换矩阵**。根据线性代数的基本性质，对任意实反对称矩阵$$\mathbf{\Omega}_n$$，其Cayley变换$$\mathbf{A}_n$$必为正交矩阵 ($$\mathbf{A}_n^{\mathrm{T}}\mathbf{A}_n = \mathbf{I}$$，$$\det\mathbf{A}_n = 1$$)，即该更新步等价于对$$\mathbf{m}_n$$施加一次纯三维旋转[^ref-iserles-lie-group]。因此，无论时间步长$$\Delta t$$取何值，都在数学上严格保证$$|\mathbf{m}_{n+1}| = |\mathbf{m}_n| = 1$$，彻底消除了坐标奇点并排除了手动归一化的需要。

**收敛性分析。** 在确定性ODE (即令$$\mathbf{H}_{\mathrm{TH}} = 0$$) 情形下，隐式中点法是经典的二阶对称Runge-Kutta方法，其局部截断误差为$$O(\Delta t^3)$$，全局误差为$$O(\Delta t^2)$$。对于包含Stratonovich白噪声的sLLG方程，收敛阶的分析须区分强收敛与弱收敛两个意义[^ref-kloeden-platen]：强收敛衡量逐条轨迹的精度，中点法的强收敛阶为1.0；弱收敛衡量统计量 (期望值、矩) 的精度，弱收敛阶为2.0，与ODE情形相同。在实际Monte Carlo仿真中，关注的目标量是翻转概率$$P_{\mathrm{sw}}$$ (统计均值)，因此弱收敛阶决定了所需步长精度，这也是在Cayley方法大步长下仍能保持统计精度的理论依据。García-Palacios和Lázaro以及d'Aquino等人[^ref-daquino-midpoint]进一步严格证明，对于Stratonovich意义下的sLLG方程，隐式中点/Cayley方法能够自然收敛于正确的Boltzmann热平衡分布，而无需引入Itô-Stratonovich修正项，这一性质在其他显式随机积分方案中通常需要额外添加噪声修正项才能满足。

由上述Cayley更新公式可见，每一步的计算仅需构造$$3 \times 3$$矩阵$$\mathbf{A}_n$$的一次求逆与矩阵-向量乘积。对于$$3 \times 3$$矩阵，逆矩阵可由解析公式直接给出而无需迭代，因此单步计算代价固定且极低，使在常规CPU平台上对$$10^5$$以上独立轨迹进行快速并行Monte Carlo扫描成为可能。

---

#### 2.2.3.3 随机过程的Monte Carlo统计计算

由于热噪声场$$\mathbf{H}_{\mathrm{TH}}$$赋予系统本征随机性，单次磁化轨迹的终态不足以描述器件的宏观概率行为。为获取基于物理底层的翻转概率，须通过大规模Monte Carlo方法进行独立采样统计。

在固定驱动配置 (外加电压$$V$$、电流密度$$J$$、脉冲宽度$$t_w$$) 下重复执行$$N$$次独立仿真，每次仿真使用不同的随机种子以确保热噪声轨迹统计独立。对第$$i$$次仿真，以$$m_z$$的终态判定是否发生翻转 (以反平行态到平行态翻转为例，要求$$m_z$$在脉冲结束后弛豫至$$m_z < -0.5$$的稳定态)，定义伯努利随机变量：

$$
s_i = \begin{cases} 1, & \text{器件发生翻转} \\ 0, & \text{器件未翻转} \end{cases}
$$

由大数定律，翻转概率的无偏估计为

$$
P_{\mathrm{sw}}(V, J, t_w) = \frac{1}{N}\sum_{i=1}^{N}s_i
$$

由于$$s_i \sim \mathrm{Bernoulli}(P_{\mathrm{sw}})$$，估计量的标准差为$$\sigma_{P} = \sqrt{P_{\mathrm{sw}}(1-P_{\mathrm{sw}})/N}$$，对应的95%置信区间半宽为

$$
\delta P = 1.96\sqrt{\frac{P_{\mathrm{sw}}(1-P_{\mathrm{sw}})}{N}}
$$

为使95%置信区间半宽$$\delta P$$不超过目标量级的一半，所需样本数须满足

$$
N > \frac{(1.96)^2 P_{\mathrm{sw}}(1-P_{\mathrm{sw}})}{\delta P^2}
$$

以分辨$$10^{-3}$$量级为例 (即$$\delta P = 5 \times 10^{-4}$$，$$P_{\mathrm{sw}} \approx 10^{-3}$$)，需要$$N \gtrsim 1.5 \times 10^4$$；若目标精度提升至$$10^{-4}$$量级，则$$N$$须增至$$1.5 \times 10^6$$，这对仿真效率提出了严格要求，也是Cayley变换方案相对于传统小步长显式积分在大步长宽容度上的核心优势得以充分发挥的实际应用背景。

---

---

### 2.2.4 开源仿真平台：vgsot-sim

前述各节建立了从sLLG方程、温度依赖材料参数到Monte Carlo统计的完整物理与算法框架。本节介绍配套的开源Python仿真平台vgsot-sim，对其分层架构、核心物理通道与参数标定方法分别说明。

vgsot-sim以Python为核心实现语言，通过PyPI发布并以pip安装，对外暴露命令行接口与Python API两种调用形式：前者用于复现单组实验或批量扫描，后者用于与系统级仿真框架与深度学习框架集成。仿真器的核心功能是对单个MTJ器件执行sLLG方程的时域积分，在每条轨迹结束后判定翻转状态，并通过大规模重复运行统计翻转概率。平台采用三层架构以保证物理模型的可扩展性与实验配置的灵活性，整体结构的组合关系可形式化表示为

$$
\text{Simulation} = \text{Kernel} \circ \text{Config} \circ \text{IO}
$$

最底层为物理内核 (Kernel)，负责实现sLLG方程的Cayley变换求解、有效场构建、热噪声生成以及温度状态更新等核心计算逻辑；中间层为实验配置 (Config)，用于声明具体仿真条件，包括脉冲参数、材料参数与扫描范围；顶层为输入输出层 (IO)，负责结果的序列化存储、统计汇总以及与外部系统的数据交换。该分层结构的核心设计原则是将物理模型与实验场景完全解耦：内核函数不携带任何与具体实验相关的状态，配置层仅通过参数对象驱动内核行为，因此同一求解器可在参数空间的不同工作点无修改地复用。IO层通过统一的结果数据结构封装磁化轨迹、翻转标志与统计量，使Monte Carlo汇总、曲线拟合以及后续分析流程均可直接调用。整体架构如图2.7所示。

![vgsot-sim三层架构示意](figs/Chapter02_local_07.png)

**图2.7** vgsot-sim仿真软件框架。用户接口层提供命令行与Python API两条等价调用路径；实验配置层以参数对象形式声明标准化测试场景与扫描范围；物理内核层按单步更新回路组织磁化动力学、有效场、热噪声、电学输运与电阻五个功能模块，各模块间的数据流在每个时间步内完成一次磁化状态、电阻与温度的自洽更新，最终对外输出磁化轨迹与翻转概率统计。

物理内核层中，磁化动力学求解模块以离散时间步为基本单元，每步内依次完成以下操作：先由有效场模块按2.2.1.2节给出的分量合成计算当前$$\mathbf{H}_{\mathrm{eff}}$$；再由热噪声模块按

$$
\langle H_{\mathrm{th},i}(t)\, H_{\mathrm{th},j}(t')\rangle
= \delta_{ij}\,\delta(t-t')\,\frac{2\alpha k_B T}{\mu_0 M_s \gamma V},
\qquad
H_{\mathrm{th}} \sim \mathcal{N}\!\left(0,\;\frac{2\alpha k_B T}{\mu_0 M_s \gamma V \Delta t}\right)
$$

采样满足涨落-耗散定理的高斯随机场；电学输运模块将外加电压映射为SOT电流密度与MTJ隧穿电流，提供驱动力矩与Joule热源；热力学模块以当前电学状态按2.2.2.1节的一维RC方程更新器件温度，并按2.2.2.2节的Bloch / Callen-Callen关系反馈至材料参数$$M_s(T)$$、$$K_i(T)$$；TMR模块根据当前磁化方向$$m_z$$与偏压$$V_{\mathrm{MTJ}}$$计算器件瞬时电阻。单步磁化更新的完整形式为

$$
\mathbf{m}_{n+1}
= \mathcal{F}\!\left(\mathbf{m}_n,\;\mathbf{H}_{\mathrm{eff},n},\;\mathbf{H}_{\mathrm{th},n}\right),
$$

其中$$\mathcal{F}$$对应2.2.3.2节所述的Cayley变换步，在数学上严格保证$$|\mathbf{m}_{n+1}|=1$$。

为使仿真结果能够直接服务于不同写入机制的对比分析，平台预置四类标准化的仿真场景：纯SOT基线场景关闭VCMA调制 ($$V_{\mathrm{MTJ}}=0$$)，其翻转概率仅由SOT电流密度与脉冲宽度决定，是建立基准曲线的出发点；VCMA辅助场景在SOT电流基础上施加MTJ偏置电压，通过$$\Delta(V)=\Delta_0-\beta_{\mathrm{VCMA}}V$$动态调低有效能垒以降低写入电流；优化双脉冲场景按2.1.3节的SOT-VCMA联合驱动模型设计两段脉冲序列，先以VCMA脉冲降低能垒、再以SOT脉冲完成翻转，用于量化能效优化收益；SER蒙特卡罗场景对每个驱动参数工作点执行$$N$$次独立轨迹并统计写错误率$$\mathrm{SER}=1-\frac{1}{N}\sum_i s_i$$或等价的翻转概率$$P_{\mathrm{sw}}=1-\mathrm{SER}$$。$$N$$的默认值随目标置信度自适应调整 (参见2.2.3.2节)。

由于行为级紧凑模型的设计目标是端口级输出与实测对齐而非保留全部材料原生常数，平台对若干参数按实测特征量进行标定。具体地，TMR$$_0$$、电阻面积积R·A与有效自旋霍尔角$$\theta_{\mathrm{SH}}$$三者联合调整为$$(1.00,\,16.6\,\Omega\!\cdot\!\mu\mathrm{m}^2,\,0.04)$$，使仿真给出的$$R_P\!\approx\!5\,\mathrm{k}\Omega$$、$$R_{AP}\!\approx\!10\,\mathrm{k}\Omega$$与$$V_{\mathrm{th}}(0.75\,\mathrm{ns})\!\approx\!903\,\mathrm{mV}$$，与2.3.2节同批次实验$$R_P,R_{AP},V_{\mathrm{th}}=(4.9\,\mathrm{k}\Omega,10\,\mathrm{k}\Omega,894\,\mathrm{mV})$$相互吻合至1%以内。需要强调的是，标定后的$$\theta_{\mathrm{SH}}\!\approx\!0.04$$是lump掉若干本模型未显式建模的耗散通道 (Néel–Edelstein界面项、自旋记忆损失、寄生串联电阻、电流方向与$$\hat{\sigma}$$轴局部偏离等) 的端口级有效值[^note-dev-thetacalib]。

---

### 2.2.5 仿真流程与可输出观测量

基于2.2.4节所述vgsot-sim平台，可对本文所研究SOT-MTJ器件在典型写入脉冲下的随机翻转行为进行时域仿真。仿真所需的器件几何、磁学、输运与热学参数取自2.2.2节各小节列出的80 nm级基准参数集；每条磁化轨迹由2.2.3.2节的Cayley变换保模长算法推进，并按2.2.3.2节的Bernoulli统计方法汇总多条独立轨迹，以获得翻转概率的无偏估计。

平台直接输出的核心观测量包括两类。其一是单条轨迹的时域演化量，即自由层归一化磁化分量$$m_z(t)$$及其通过2.2.2.3节给出的TMR与电输运模型换算得到的MTJ瞬时电阻$$R_{\mathrm{MTJ}}(t)$$，用以刻画在给定驱动脉冲下器件磁化状态随时间的随机演化以及读出端的瞬态电阻响应。典型单次事件的输出形式如图2.8所示，从上至下依次为$$m_z(t)$$、$$R_{\mathrm{MTJ}}(t)$$与用于驱动的$$I_{\mathrm{SOT}}(t)$$波形；通过在相同脉冲宽度下扫描不同$$I_{\mathrm{SOT}}$$幅值，可直接观察到磁化翻转发生与否、翻转时刻以及电阻跳变幅度等关键行为特征。

![单次m_z与R_MTJ演化事件](figs/Chapter02_local_08.png)

**图2.8** vgsot-sim在$$t_w = 0.75\,\mathrm{ns}$$写入脉冲下的单次轨迹输出。(a)归一化磁化分量$$m_z(t)$$。(b)由TMR模型换算的瞬时MTJ电阻$$R_{\mathrm{MTJ}}(t)$$ ($$R_P\!\approx\!5\,\mathrm{k}\Omega$$、$$R_{AP}\!\approx\!10\,\mathrm{k}\Omega$$，与图2.12滞回回线幅度一致)。(c)SOT驱动电流脉冲$$I_{\mathrm{SOT}}(t)$$。初始态为PAP=1 ($$m_z\approx-1$$) 、热噪声NON=1、自热反馈开启 (详见2.2.2.5节) ；仿真使用2.2.4节校准至Device A P→AP @ 0.75 ns实验阈值的有效$$\theta_{\mathrm{SH}}=0.04$$，并以代表性RNG种子使各$$I_{\mathrm{SOT}}$$级展现其在SER MC分布中的最可能行为。四条$$I_{\mathrm{SOT}}\in\{-600,-1100,-1300,-2000\}\,\mu\mathrm{A}$$跨越亚阈值、临界、刚翻转与确定性翻转四种情形：600 µA下$$m_z$$维持在$$-1$$不动 ($$R_{\mathrm{MTJ}}\!\approx\!R_{AP}$$) ；1100 µA接近阈值但脉冲关断后仍沿$$-z$$方向回落；1300 µA处出现成功跨越赤道并落入$$+z$$基态的翻转事件 ($$R_{\mathrm{MTJ}}$$跃迁至$$R_P$$) ；2000 µA给出更快的赤道达到时刻。临界电流1100–1300 µA区间与实验$$I_{\mathrm{th}}(0.75\,\mathrm{ns})\approx 1152\,\mu\mathrm{A}$$ ($$V_{\mathrm{th}}\approx 894\,\mathrm{mV}$$) 量纲匹配。

为补充$$m_z(t)$$标量视图，平台同时记录每步的极角$$\theta(t)$$与方位角$$\phi(t)$$并据此重构完整磁化矢量$$\mathbf{m}(t)=(\sin\theta\cos\phi,\,\sin\theta\sin\phi,\,\cos\theta)$$。图2.9在单位球面上给出一条典型超阈值翻转事件的三维轨迹，左面板以时间为色标显示磁化矢量从反平行极 ($$m_z=-1$$) 沿赤道附近螺旋进动并最终收敛至平行极 ($$m_z=+1$$) 的全过程，右面板同步呈现三个笛卡儿分量的时域演化。三维球面视图直接揭示了亚纳秒SOT写入下磁化翻转所特有的螺旋进动结构：相干进动周期与界面各向异性场$$H_k$$给出的Larmor频率一致，进动衰减包络受Gilbert阻尼控制，赤道附映射均以近的随机抖动来自Brown热涨落，三者共同构成sLLG动力学的完整可视化。

![单次磁化矢量在单位球面上的三维轨迹](figs/Chapter02_local_09.png)

**图2.9** vgsot-sim在$$I_{\mathrm{SOT}}=-2000\,\mu\mathrm{A}$$、$$t_w=0.75\,\mathrm{ns}$$确定性翻转条件下的单次磁化矢量轨迹。(a)在单位球面上以时间为色标显示完整$$\mathbf{m}(t)$$演化路径，蓝色圆点表示初态反平行极，金色五角星表示终态平行极，紫色细线为球面网格仅作几何参考。轨迹在初始的SOT驱动阶段 (0–0.75 ns) 沿赤道附近螺旋上升，进动周期约0.2 ns，与$$H_k$$对应的Larmor频率量级一致；脉冲关断后磁化在剩余3.25 ns弛豫窗口内沿Gilbert阻尼通道收敛至上极。(b)同次仿真的笛卡儿分量$$m_x(t)$$、$$m_y(t)$$、$$m_z(t)$$时域演化，赤道附近的高频振荡周期与三维视图所呈现的螺旋间距一致，$$m_z$$穿越零点的时刻对应轨迹跨越赤道；该视角与图2.8的多电流情形形成互补，前者刻画进动几何，后者刻画统计行为。

平台的第二类输出量是相同驱动配置下多次独立仿真的统计汇总结果，即翻转概率$$P_{\mathrm{sw}}$$ (或互补的写错误率$$\mathrm{SER}\equiv 1-P_{\mathrm{sw}}$$) 随SOT驱动电流$$I_{\mathrm{SOT}}$$、脉冲宽度$$t_w$$以及VCMA辅助电压$$V_{\mathrm{MTJ}}$$的变化曲线，用以刻画概率翻转窗口在不同工作点上的宏观响应特性。本节以$$P_{\mathrm{sw}}$$表述以便与2.3.2节 同批次实验Sigmoid直接同量纲对比；仿真器以 `--metric=psw|ser` 开关切换两种表示。图2.10给出一组典型的Monte Carlo扫描结果，即0.75 ns写入脉冲宽度下$$P_{\mathrm{sw}}$$随$$I_{\mathrm{SOT}}$$的变化曲线，并在同一图上对比了自热反馈关闭与开启两种情形，是将器件模型与后续阵列级概率计算评估相连接的直接接口。

![Monte Carlo Psw扫描结果](figs/Chapter02_local_10.png)

**图2.10** $$t_w = 0.75\,\mathrm{ns}$$写入脉冲 + 3.25 ns弛豫窗口下输出的$$P_{\mathrm{sw}}$$–$$|I_{\mathrm{SOT}}|$$蒙特卡罗扫描结果。(a)宽范围扫描 (300–3500 µA，每点80条独立轨迹，Wilson 95% 置信区间以阴影带给出)，蓝色实线为自热反馈关闭、红色虚线为自热开启；青色虚线标示2.3.2节Sigmoid拟合给出的实验阈值$$I_{\mathrm{th}}=V_{\mathrm{th}}/R_W=894\,\mathrm{mV}/776\,\Omega\!\approx\!1152\,\mu\mathrm{A}$$。亚阈值区 ($$|I_{\mathrm{SOT}}|\!\lesssim\!900\,\mu\mathrm{A}$$) $$P_{\mathrm{sw}}\!\approx\!0$$；过渡区1000–1300 µA内陡升至约0.80，与实验Sigmoid形态吻合；超阈值区($$\geq 1500\,\mu\mathrm{A}$$)进入由back-hopping上限主导的$$\sim 0.8$$平台。(b)阈值区精扫描inset (800–1400 µA，实验$$I_{\mathrm{th}}=1152\,\mu\mathrm{A}$$附近密集取点、两侧稀疏，每点80条)，琥珀色曲线呈现清晰的Sigmoid型过渡。仿真50%翻转点与实验阈值一致至5%以内。自热开启支在1100 µA工作点比关闭支高0.21 $$P_{\mathrm{sw}}$$ ($$\Delta T_{\mathrm{eq}}\approx 36\,\mathrm{K}$$对应$$\Delta K_i/K_i\approx-7\%$$，使阈值在统计意义下向左偏移) ；超阈值区$$|\Delta P_{\mathrm{sw}}|$$落入Wilson带宽内不可识别。

上述观测量所对应的具体数值结果，以及基于这些结果对概率翻转模型、自热效应、VCMA辅助写入能效以及工艺波动影响等方面开展的量化分析，将在后续章节中结合具体仿真器架构与实验对比进行讨论。

## 2.3 sMTJ器件实验验证与联合写入概率模型

本节通过直接测量sMTJ单器件的电学特性并结合300 mm晶圆工艺数据，建立从物理器件到概率计算单元的可靠映射。前文已从微磁动力学模型与行为级建模角度系统分析了SOT器件在热涨落作用下的随机翻转特性，但仿真方法难以全面反映真实器件中的工艺波动、材料缺陷、界面粗糙度以及寄生电阻网络等非理想因素，这些因素均会直接影响翻转概率分布及其稳定性。本节通过实验测量对行为级模型进行参数校准，并进一步评估工艺波动对阵列级概率一致性的影响。

在基于MRAM的概率计算体系中，热激活驱动的随机磁化翻转过程构成伯努利分布的物理随机源，器件的翻转概率$$P_{\mathrm{sw}}(V, t)$$作为可控计算资源被主动利用[^ref-borders-factorization]。该计算范式对器件提出特定工作要求，例如，在亚阈值或近阈值偏置条件下翻转概率应呈现稳定、连续且可精确建模的电压响应，器件还应具备足够的耐久性与晶圆级工艺一致性以支撑阵列规模下的统计均匀性。

---

![SOT-MTJ器件实验表征综合图](figs/Chapter02_local_11.png)

**图2.11** SOT-MTJ器件实验平台与统计表征综合视图。(a)高速测试系统原理框图：超快电压脉冲经功率分配器分为两路 (上路可选$$-6\,\mathrm{dB}$$衰减)，射频偏置器合并高频脉冲与10 mV直流偏置后施加于器件顶层或底层电极，定向耦合器监测信号状态，SMU在顶层电极处采集电流响应。(b)重金属钨(W)层电阻率与自旋霍尔角随厚度的变化关系，紫色阴影区域(4–5 nm)标示了兼顾自旋电荷转化效率与沟道电阻的最佳厚度窗口。(c)高周疲劳耐久性测试结果，$$R_{AP}$$、$$R_P$$和$$R_{\mathrm{SOT}}$$在超过$$10^{11}$$次翻转循环后仍保持极高的稳定性，无明显退化。(d)实物照片：芯片样品、器件阵列光学显微镜照片及探针台测试系统。

### 2.3.1高速测试系统架构与器件信息

本研究所采用的器件为驰拓(HIKSTOR)在300 mm晶圆工艺平台上实现的三端SOT-MTJ，MTJ柱采用top-pinned堆叠 ($$\mathrm{CoFeB}/\mathrm{MgO}/\mathrm{CoFeB}/\mathrm{spacer}/\mathrm{SAF}$$) 、标称直径80 nm；SOT通道为$$\beta$$-W薄膜；器件几何与隧穿/输运参数的完整列表参见2.1.1节T型电路定义与2.2.2节仿真参数表 (含$$R\!\cdot\!A=36\,\Omega\!\cdot\!\mu\mathrm{m}^2$$、$$\mathrm{TMR}=100\text{--}120\%$$、$$\theta_{\mathrm{SH}}=0.25$$等标称量)，本节不再重复。下面侧重测试平台架构与器件耐久性、阵列均值等实测特征。

为实现纳秒尺度器件动态行为的精确测量，本实验构建了一套高速电学测试平台，核心组件包括超快脉冲电压源、高带宽示波器、射频探针台及高精度源表。脉冲电压源可产生最小宽度达亚纳秒级的写入脉冲，示波器实时捕捉测试链路中的电压与电流波形。射频链路引入功率分配器、衰减器及定向耦合器实现信号调制与测量隔离，并通过偏置器将直流偏置与高频脉冲叠加以支持多模式测试。

测试链路中超快电压脉冲经功率分配器分为两路，上路可选择性接入$$-6\,\mathrm{dB}$$衰减器；射频偏置器(RF Bias Tee)将高频脉冲与来自Keysight B2901源表(SMU)的10 mV直流偏置合并后施加于被测器件的顶层或底层电极，定向耦合器实时监测链路信号状态，SMU的电流测量端则接于顶层电极以捕获器件响应。射频线缆、连接器与探针存在频率响应，输入脉冲在传输过程中会发生幅值衰减与波形畸变，因此需在断开探针条件下测量脉冲源输出与示波器接收信号之间的转换关系，提取电压校正系数以恢复实际施加于器件端口的有效电压，保证测量结果的定量准确性。

写入路径中电流经SOT通道注入，利用自旋霍尔效应在自由层中产生垂直自旋流从而驱动磁矩翻转。为保证翻转方向的确定性，沿电流方向施加约200 Oe量级的外部面内磁场以打破自由层面内对称性并控制翻转手性，该方法在三端SOT-MRAM实验研究中已被广泛证明能够显著改善写入的方向一致性[^ref-grimaldi-sot-mtj]。

**表2.5** 器件阵列测试统计均值

| 特征 | 测量数值 |
|:-----|:--------|
| 矫顽场$$\mu_0 H_c$$ | 75 mT |
| 偏置场$$\mu_0 H_{\mathrm{offset}}$$ | 0.36 mT |
| MTJ平行态电阻$$R_P$$ | 10.89 kΩ |
| SOT沟道电阻$$R_{\mathrm{SOT}}$$ | 776 Ω |

需要指出，2.2.2节表中所给$$R\!\cdot\!A=36\,\Omega\!\cdot\!\mu\mathrm{m}^2$$与表2.5所列$$R_P=10.89\,\mathrm{k}\Omega$$之间存在表观差异：若简单以物理直径$$D_{\mathrm{phys}}=80\,\mathrm{nm}$$推算几何面积$$A_{\mathrm{phys}}=\pi D_{\mathrm{phys}}^2/4\approx5.03\times10^{-3}\,\mu\mathrm{m}^2$$，对应平行态电阻为$$R\!\cdot\!A/A_{\mathrm{phys}}\approx7.16\,\mathrm{k}\Omega$$，与实测均值低约34%。两者的差异主要源于电学有效直径与物理直径不重合：离子铣刻、再沉积及侧壁损伤使MTJ柱靠近边缘约$$5\text{--}10\,\mathrm{nm}$$环带的电学有效性显著降低，按2.2.2.4节定义可记电学有效直径$$D_{\mathrm{elec}}\approx D_{\mathrm{phys}}-2\delta_{\mathrm{edge}}$$。反演表2.5的实测值，所需的电学有效面积为$$A_{\mathrm{elec}}=R\!\cdot\!A/R_P\approx3.31\times10^{-3}\,\mu\mathrm{m}^2$$，对应$$D_{\mathrm{elec}}\approx64.9\,\mathrm{nm}$$，与上述刻蚀损伤估计自洽。因此，本工作中Brinkman--Dynes--Rowell模型与$$R_{\mathrm{MTJ}}(m_z)$$映射均以$$D_{\mathrm{elec}}$$代入计算，以保持电学量与几何参数之间的自洽。

---

### 2.3.2 sMTJ器件写入特性测量结果

本节基于同一测量窗口下在Device A与Device B两个器件上取得的完整数据集，依次给出临界电压随脉宽的对数依赖、单次写入能耗，以及$$t_w = 0.75\,\mathrm{ns}$$条件下四组概率翻转曲线。

**临界写入电压的脉宽依赖。** 在固定外加面内磁场$$H_x = 200\,\mathrm{Oe}$$条件下，施加脉冲宽度$$t_w = 0.75\,\mathrm{ns}$$、$$1\,\mathrm{ns}$$、$$2\,\mathrm{ns}$$、$$5\,\mathrm{ns}$$的写入脉冲，在每个脉宽下扫描脉冲电压幅值$$V_{\mathrm{SOT}}\in[-1.1, 1.1]\,\mathrm{V}$$，获取完整的电阻-电压滞回回线。每个扫描点先将器件初始化至已知磁化状态，施加单次写入脉冲，随后通过SMU读取MTJ电阻状态以判断翻转是否发生。

实验滞回回线如图2.12(a)所示。器件在正向脉冲作用下由反平行态($$R_{AP}\approx 10\,\mathrm{k\Omega}$$)翻转至平行态($$R_P\approx 4.9\,\mathrm{k\Omega}$$)，在负向脉冲作用下则由平行态翻转至反平行态。随脉冲宽度减小翻转所需的电压阈值单调增大。Device A的5 ns脉冲下正向翻转阈值$$V_{\mathrm{th}+}\approx 0.56\,\mathrm{V}$$，在亚纳秒的0.75 ns条件下升至约0.90 V；反向阈值$$V_{\mathrm{th}-}$$呈现对称变化规律。两方向阈值在绝对值上存在约3%–5%的差异，与阵列级测得的参考层偏置场$$\mu_0 H_{\mathrm{offset}} = 0.36\,\mathrm{mT}$$对两态能垒的非对称调制一致。

由滞回回线的电阻跳变点提取各脉宽下的正负向临界电压后，以$$\ln(t_w)$$为自变量绘制$$V_{\mathrm{th}}$$的依赖关系，得到如图2.12(b)所示的严格线性趋势。对Device A两个方向分别进行线性回归得到拟合关系式

$$
V_{\mathrm{AP\to P}}(t_w) = 0.82 - 0.17\ln(t_w/\mathrm{ns})\ \mathrm{V},\qquad V_{\mathrm{P\to AP}}(t_w) = -0.79 + 0.18\ln(t_w/\mathrm{ns})\ \mathrm{V},
$$

两个方向的决定系数$$R^2$$均在0.995以上，证实在所测0.75–5 ns范围内临界写入电压与脉冲宽度的对数呈严格线性依赖。该对数依赖关系是热激活翻转区的核心特征，是后续反推Néel-Brown模型参数的直接依据。

**写入能耗估算。** SOT-MRAM的写入电流主要流经重金属通道，单次写入能耗由SOT通道上的欧姆耗散主导。以最短脉宽$$t_w = 0.75\,\mathrm{ns}$$、正向临界写入电压$$V_{\mathrm{th}+}\approx 0.90\,\mathrm{V}$$为例，结合SOT通道电阻$$R_{\mathrm{SOT}}\approx 776\,\Omega$$，单次写入能耗可估算为
$$
E_{\mathrm{write}} = \frac{V^2}{R_{\mathrm{SOT}}}\cdot t_w = \frac{(0.90\,\mathrm{V})^2}{776\,\Omega}\times 0.75\,\mathrm{ns}\approx 0.78\,\mathrm{pJ}.
$$

对应电流密度约$$J_c\approx 1.0\times 10^{12}\,\mathrm{A/m^2}$$，与文献中先进SOT-MTJ器件在亚纳秒脉冲下的典型值处于同一量级。该能耗水平满足嵌入式高速存储与边缘概率计算场景对写入功耗的约束。

**$$t_w = 0.75\,\mathrm{ns}$$概率翻转特性。** 固定脉冲宽度为$$t_w = 0.75\,\mathrm{ns}$$，对Device A和Device B各以正、负两个方向扫描写入电压，每个幅值独立重复执行100次写入操作，以成功翻转次数占比定义翻转概率$$P_{\mathrm{sw}}$$。该测量与前述滞回扫描在同一批次内完成 (同一器件、同一连续测试窗口)，构成与NB参数反推直接对应的基准数据。四条$$P_{\mathrm{sw}}(V)$$曲线如图2.13所示，对每条曲线独立进行四参数Sigmoid拟合
$$
P_{\mathrm{sw}}(V) = y_0 + \frac{L}{1+\exp[-(|V|-V_{\mathrm{th}})/k]},
$$

其中$$k$$为尺度参数，对应的logistic斜率参数$$\beta_s = 1/k$$。拟合结果列于表2.6。

**表2.6** $$t_w = 0.75\,\mathrm{ns}$$的四条Sigmoid拟合结果

| 器件 | 方向 | $$V_{\mathrm{th}}$$ (mV) | $$k$$ (mV) | $$\beta_s$$ (V⁻¹) | $$R^2$$ | 备注 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| A | AP→P | 915.5 | 5.45 | 183.5 | 0.969 | $$V\gtrsim 940$$ mV出现回跳平台，拟合受扰 |
| A | P→AP | 894.0 | 22.43 | 44.6 | 0.993 | 干净单段过渡 |
| B | AP→P | 865.9 | 36.54 | 27.4 | 0.968 | 840–860 mV附近出现两段过渡 |
| B | P→AP | 904.5 | 11.26 | 88.8 | 0.995 | 干净单段过渡 |

Device A的AP→P曲线在$$V > 940\,\mathrm{mV}$$出现明显回跳平台，$$P_{\mathrm{sw}}$$在900–1000 mV间维持约0.72直到$$V \gtrsim 1020\,\mathrm{mV}$$才重新上升至1。该现象在高电压SOT翻转中已被报告为back-hopping机制：SOT脉冲结束后磁矩因热辅助再次跃过势垒回到初态，导致名义翻转率饱和在低于1的水平。该非理想行为使单一Sigmoid拟合在高V区偏离数据，$$R^2 = 0.969$$低于其他三条，表中Device A AP→P方向的$$\beta_s$$值偏大且不确定度较高。Device B的AP→P曲线在840–860 mV附近出现约20 mV宽的平台，$$P_{\mathrm{sw}}$$在0.36–0.38之间几乎不变，可能与亚畴先后分步翻转相关。两条P→AP曲线 (Device A与Device B) 均表现出干净的单段过渡，$$R^2$$均达0.993以上，是用于Sigmoid参数提取的主力数据。

综合以上四组曲线的拟合质量与物理一致性，本节采用Device A、P→AP、$$t_w = 0.75\,\mathrm{ns}$$的Sigmoid拟合结果作为后续联合模型的主基准($$V_{\mathrm{th}} = 894\,\mathrm{mV}$$，$$\beta_s = 44.6\,\mathrm{V}^{-1}$$，$$R^2 = 0.993$$)。该选择基于两点考虑。Device A为参数反推依赖的主器件，可保证拟合所得$$\Delta$$、$$V_{c0}$$与Sigmoid参数来自同一物理对象。P→AP方向数据无非理想行为干扰，拟合精度最高。四条曲线的整体图景同时覆盖AP→P与P→AP两方向、Device A与Device B两器件，共同提供了Sigmoid形态的器件间与方向间变异性特征，为后续定量比较与工艺容差分析提供完整证据基础。

---

![sMTJ写入特性与Néel-Brown联合概率模型](figs/Chapter02_local_12.png)

**图2.12** sMTJ器件在不同脉冲宽度下的写入特性与Néel-Brown联合概率模型。(a)$$t_w = 0.75\,\mathrm{ns}$$、$$1\,\mathrm{ns}$$、$$2\,\mathrm{ns}$$、$$5\,\mathrm{ns}$$条件下测得Device A的电阻-脉冲电压滞回回线，随脉冲持续时间缩短翻转电压窗口逐渐展宽。(b)临界翻转电压$$V_{\mathrm{th}\pm}$$随脉冲宽度的对数依赖关系，空心符号为实验数据点，实线为对数线性拟合$$V = a\mp b\ln(t_w/\mathrm{ns})$$；内嵌注释给出由$$\tau_0 = 1\,\mathrm{ns}$$先验反推得到的两方向Néel-Brown参数$$(\Delta, V_{c0}, \tau_{\mathrm{ret}})$$。(c)基于反推NB参数构建的二维联合翻转概率分布$$P_{\mathrm{sw}}(V, t_w)$$热力图(AP→P方向)，紫色等概率轮廓在低概率区可读，白色等概率轮廓在高概率区可读，黑色虚线为50%等概率轨迹即$$V_{\mathrm{th}}(t_w)$$；空心圆(Device A)与三角(Device B)标示两器件的滞回提取点，与50%轨迹吻合。

---

![Sigmoid测量与Néel-Brown外推对比](figs/Chapter02_local_13.png)

**图2.13** $$t_w = 0.75\,\mathrm{ns}$$、$$H_x = 200\,\mathrm{Oe}$$条件下100次重复Sigmoid测量与C2C-修正后的Néel-Brown模型对比。(a)Device A AP→P：实测在$$V \gtrsim 940\,\mathrm{mV}$$出现back-hopping回跳平台。(b)Device A P→AP：干净单段过渡($$R^2 > 0.99$$)，作为主基准曲线。(c)Device B AP→P：840–860 mV附近出现两段过渡。(d)Device B P→AP：干净单段过渡($$R^2 > 0.99$$)。空心符号为实验数据点(误差线为二项分布的Wilson 95%置信区间)，实线为C2C-修正NB曲线 (数学上等价于四参数Sigmoid拟合)，各面板内嵌注释给出$$\eta_c = \beta_s/\beta^{\mathrm{NB}}$$与该曲线的物理特征。Sigmoid拟合对四条曲线的$$V_{\mathrm{th}}$$预测精度均优于+7.4%，但未修正NB预测的斜率(约8 V⁻¹)普遍低于实测数倍，必须以$$\eta_c$$因子作C2C修正方能重现实测分布陡度。

---

![同批次器件间Néel-Brown参数一致性对比](figs/Chapter02_local_14.png)

**图2.14** Device A 与 Device B 两个器件的Néel-Brown参数一致性对比。(a)正向(AP→P)临界翻转电压$$V_{\mathrm{th}+}$$随脉冲宽度$$t_w$$的对数线性依赖；Device A (红色圆点) 与Device B (紫色三角) 数据点近似落在同一条对数直线上，两器件的$$\Delta$$与$$V_{c0}$$拟合值在5%–15%范围内一致。(b)以$$\tau_0 = 1\,\mathrm{ns}$$先验反推得到的两方向热稳定性因子$$\Delta$$柱状图，AP→P (红色) 与P→AP (蓝色) 方向在同一器件上数值接近(器件A的两方向$$\Delta$$分别为5.15与4.91、器件B分别为4.46与4.95)，对应零温临界电压$$V_{c0}$$分别为884 mV (器件A) 与876 mV (器件B)，证实sMTJ作为概率单元具备良好的器件间均匀性，为后续阵列级建模与工艺容差分析提供同批次基线参考。

### 2.3.3 Néel-Brown模型参数的反推

前述$$V_{\mathrm{th}}$$对数线性依赖关系与前文建立的Néel-Brown翻转概率模型在数学上严格等价，可由实验系数反推底层物理参数。取线性势垒近似$$n=1$$，翻转概率满足

$$
P_{\mathrm{sw}}(t_w, V) = 1 - \exp\!\left[-\frac{t_w}{\tau_0}\exp\!\left(-\Delta\left(1-\frac{V}{V_{c0}}\right)\right)\right]
$$

其中$$\tau_0$$为尝试频率倒数、$$\Delta = E_b/k_BT$$为热稳定性因子、$$V_{c0}$$为零温临界电压。令$$P_{\mathrm{sw}}(t_w, V_{\mathrm{th}}) = 0.5$$解得

$$
V_{\mathrm{th}}(t_w) = V_{c0}\left[1 - \frac{1}{\Delta}\ln\!\left(\frac{t_w}{\tau_0\ln 2}\right)\right] = \underbrace{\left[V_{c0} + \frac{V_{c0}}{\Delta}\ln(\tau_0\ln 2)\right]}_{a} - \underbrace{\frac{V_{c0}}{\Delta}}_{b}\ln(t_w)
$$

比对实验拟合$$V_{\mathrm{th}} = a - b\ln(t_w)$$的两个系数可建立$$(\tau_0,\Delta, V_{c0})$$与$$(a,b)$$之间的映射关系。

**模型无关的物理可观测量。** 实验对数线性拟合提供两个独立约束，可直接转换为两个不依赖$$\tau_0$$先验假设的物理量。对数斜率$$b = V_{c0}/\Delta$$给出脉宽每增加一个$$e$$倍所需降低的写入电压，是器件翻转灵敏度的量纲化度量。令$$V_{\mathrm{th}}(t^*) = 0$$即可解得纯热涨落下的50%自发翻转时间$$t^* = \exp(a/b)$$，由此定义零驱动保持时间
$$
\tau_{\mathrm{ret}} \equiv \tau_0\, e^{\Delta} = \frac{t^*}{\ln 2} = \frac{1}{\ln 2}\exp\!\left(\frac{a}{b}\right).
$$

该量仅由拟合截距与斜率共同决定，无需假设$$\tau_0$$。Device A两方向的提取结果如表2.7所示。

**表2.7** Device A模型无关物理可观测量

| 方向 | 对数斜率$$b = V_{c0}/\Delta$$ | $$t^*$$ | 零驱动保持时间$$\tau_{\mathrm{ret}}$$ |
|:---:|:---:|:---:|:---:|
| AP→P | 172 mV | 120 ns | 172 ns |
| P→AP | 175 mV | 94 ns | 135 ns |

AP→P方向的$$\tau_{\mathrm{ret}}$$较P→AP方向长约27%，AP态自发稳定性略高于P态，与参考层杂散偶极场对两态能垒的差异性调制相一致。两方向$$\tau_{\mathrm{ret}}$$均处于百纳秒量级，远低于传统存储MRAM所要求的年量级保持时间[^note-retention-delta]，是sMTJ作为低势垒概率采样单元的典型物理特征[^ref-camsari-pbits]。

**热稳定性因子与零温临界电压。** Néel-Brown模型含三个自由参数而对数线性拟合仅提供两个约束，系统欠定，存在一族以$$\tau_0$$为参数的合法解，族内所有解给出相同的$$\tau_{\mathrm{ret}}$$和$$V_{c0}/\Delta$$。若要进一步分离$$\Delta$$与$$V_{c0}$$的绝对数值需引入$$\tau_0$$的先验约束。对于CoFeB/MgO基自由层，文献中报道的尝试频率$$f_0 = 1/\tau_0$$处于1–10 GHz量级[^ref-brown-thermal]，本文采用$$\tau_0 = 1\,\mathrm{ns}$$作为标准假设值，得器件Néel-Brown模型参数如表2.8所示。

**表2.8** 反推得到的Device A Néel-Brown模型参数($$\tau_0 = 1\,\mathrm{ns}$$)

| 参数 | AP→P | P→AP | 单位 | 物理含义 |
|:---:|:---:|:---:|:---:|:---|
| $$\tau_0$$ | 1.0 | 1.0 | ns | 尝试频率倒数(先验假设) |
| $$\Delta$$ | 5.15 | 4.91 | — | 热稳定性因子 |
| $$V_{c0}$$ | 884 | 857 | mV | 零温临界电压 |
| $$E_b = \Delta k_BT$$ | 133 | 127 | meV | 能垒高度(300 K) |
| $$\tau_{\mathrm{ret}}$$ | 172 | 135 | ns | 零驱动保持时间 |
| $$\beta_s^{\mathrm{NB}} = 2\Delta\ln 2/V_{c0}$$ | 8.08 | 7.94 | V⁻¹ | NB模型预测的Sigmoid斜率 |

两方向反推得到的$$\Delta$$均约为5，对应能垒$$E_b\approx 130\,\mathrm{meV}$$，处于sMTJ概率工作区间的典型参数窗口，可充分抑制短时自发翻转以保证可编程性，又远低于存储MRAM的势垒(约1.5 eV)以实现纳秒级热激活响应。

**映射至有效各向异性场。** 反推得到的$$\Delta_0$$可进一步映射至器件微观物理量以验证其物理合理性。在PMA构型的宏自旋近似下，自由层能垒满足$$E_b = K_{\mathrm{eff}}V_{\mathrm{mag}} = \frac{1}{2}\mu_0 H_k^{\mathrm{eff}} M_s V_{\mathrm{mag}}$$，其中$$K_{\mathrm{eff}}$$为有效单轴各向异性能密度、$$H_k^{\mathrm{eff}}$$为有效各向异性场、$$M_s$$为饱和磁化强度、$$V_{\mathrm{mag}}$$为自由层磁性体积。反解得

$$
H_k^{\mathrm{eff}} = \frac{2\Delta_0 k_BT}{\mu_0 M_s V_{\mathrm{mag}}}
$$

代入实验参数(以AP→P方向为例)，$$\Delta_0 = 5.15$$对应$$E_b = 133\,\mathrm{meV} = 2.13\times 10^{-20}\,\mathrm{J}$$；自由层几何取标称值$$D = 80\,\mathrm{nm}$$、$$t_f = 1.4\,\mathrm{nm}$$得$$V_{\mathrm{mag}} = \pi D^2 t_f/4 \approx 7.0\times 10^{-24}\,\mathrm{m}^3$$；取典型CoFeB饱和磁化强度$$M_s \approx 1.0\times 10^6\,\mathrm{A/m}$$，得$$H_k^{\mathrm{eff}} \approx 4.8\times 10^3\,\mathrm{A/m} \approx 60\,\mathrm{Oe}$$。

该$$H_k^{\mathrm{eff}}$$值显著低于传统存储型MRAM器件的，与器件被设计为低势垒sMTJ的物理定位一致。$$H_k^{\mathrm{eff}}$$亦低于阵列级测得的准静态矫顽场，原因在于本文提取$$\Delta_0$$所依赖的实验数据在200 Oe面内破对称场条件下测得，而$$H_c$$在无面内偏置的准静态条件下测量，两者描述不同工作配置下的能垒，存在差异属于预期结果。典型CoFeB参数先验下得到的$$H_k^{\mathrm{eff}}$$量级与已报道的低势垒SOT-MTJ器件一致，确认了行为级模型反推参数的物理自洽性。

### 2.3.4 Sigmoid实测与NB外推的定量比较

将前节反推得到的Néel-Brown参数代入NB阈值公式，外推至$$t_w = 0.75\,\mathrm{ns}$$概率测量点，可对NB模型在跨观测量一致性方面进行严格检验。Sigmoid斜率由$$\beta_s^{\mathrm{NB}} = 2\Delta\ln 2/V_{c0}$$给出，仅取决于反推的$$\Delta$$与$$V_{c0}$$而与$$t_w$$无关。四条曲线的NB外推值与Sigmoid拟合结果逐项对比列于表2.9。

**表2.9** $$t_w = 0.75\,\mathrm{ns}$$下NB外推与Sigmoid实测的逐项对比

| 器件 | 方向 | $$V_{\mathrm{th}}^{\mathrm{NB}}$$ | $$V_{\mathrm{th}}^{\mathrm{meas}}$$ | $$V_{\mathrm{th}}$$偏差 | $$\beta_s^{\mathrm{NB}}$$ | $$\beta_s^{\mathrm{meas}}$$ | $$\eta_c = \beta_s^{\mathrm{meas}}/\beta_s^{\mathrm{NB}}$$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| A | AP→P | 870 mV | 915.5 mV | +5.2% | 8.08 V⁻¹ | 183.5 V⁻¹ | 22.7 |
| A | P→AP | 843 mV | 894.0 mV | +6.0% | 7.94 V⁻¹ | 44.6 V⁻¹ | 5.6 |
| B | AP→P | 861 mV | 865.9 mV | +0.6% | 7.06 V⁻¹ | 27.4 V⁻¹ | 3.9 |
| B | P→AP | 842 mV | 904.5 mV | +7.4% | 8.02 V⁻¹ | 88.8 V⁻¹ | 11.1 |

表中体现两类性质不同的对比结果。NB模型对四条曲线的$$V_{\mathrm{th}}$$预测偏差仅+0.6%至+7.4%，处于典型热激活拟合的外推精度之内，由实验上验证NB模型能够正确描述器件的均值翻转行为，即翻转电压对脉宽的$$\ln(t_w)$$标度关系。NB模型对四条曲线的$$\beta_s$$预测均为约8 V⁻¹[^note-nb-slope]，而实测$$\beta_s$$跨四条曲线从27.4到183.5 V⁻¹不等，定义的C2C经验收窄因子$$\eta_c = \beta_s^{\mathrm{meas}}/\beta_s^{\mathrm{NB}}$$范围为3.9–22.7、中位数8.3，呈现NB模型对概率分布宽度的系统性过度估计。同平台的$$t_w = 5\,\mathrm{ns}$$ Sigmoid测量给出$$\beta_s = 56.9\,\mathrm{V}^{-1}$$、$$\eta_c \approx 7.0$$，与中位数一致。$$\eta_c$$在跨脉宽测量中保持同一量级，支持其作为器件级C2C分布形态参数的稳定性。

**NB框架下的$$\Delta$$参数反推。** NB模型表达式$$P_{\mathrm{sw}}(t_w, V) = 1 - \exp[-(t_w/\tau_0)\exp(-\Delta(1-V/V_{c0}))]$$包含三个待定参数$$(\tau_0, \Delta, V_{c0})$$。实验为$$\Delta$$提取提供了两条独立的相空间约束，由这两类约束分别尝试反解$$\Delta$$时所得到的数值不同。

脉宽扫描法以$$V_{\mathrm{th}}$$对数斜率$$b = V_{c0}/\Delta$$与$$\tau_0$$先验联立解出$$\Delta$$，记为$$\Delta_{\mathrm{pulse}}$$。该路径直接验证NB模型最具鉴别力的标度律，即对数时窗依赖$$V_{\mathrm{th}} \propto \ln(t_w)$$，并在所测$$0.75\text{-}5\,\mathrm{ns}$$窗口内由$$R^2 > 0.995$$的严格线性度独立检验通过。所提取的$$\Delta_{\mathrm{pulse}}$$严格对应单畴热激活势垒高度$$E_b/k_BT$$，是有明确热力学含义的物理量。在Device A、P→AP方向得$$\Delta_{\mathrm{pulse}} = 4.91$$，对应$$E_b \approx 127\,\mathrm{meV}$$、零驱动保持时间$$\tau_{\mathrm{ret}} = \tau_0 e^{\Delta} \approx 135\,\mathrm{ns}$$，与器件实际表现的纳秒级热激活响应一致。

$$P_{\mathrm{sw}}$$斜率法将NB预测的Sigmoid斜率$$\beta_s^{\mathrm{NB}} = 2\Delta\ln 2/V_{c0}$$与实测斜率$$\beta_s^{\mathrm{meas}}$$等同，反解得$$\Delta_{\mathrm{slope}} = \beta_s^{\mathrm{meas}} V_{c0}/(2\ln 2)$$。该步骤隐含一项重要附加假设，即实验C2C分布严格服从NB双指数函数对应的Gumbel形态。后文将分析的多种微观机制 (亚畴协同跃迁、热辅助进动翻转过渡、尝试频率的弱电压依赖) 均使C2C分布相对Gumbel基线收窄，故$$\Delta_{\mathrm{slope}}$$实际上是热激活势垒与C2C分布锐化效应的混合参数。仍以Device A、P→AP为例，$$\Delta_{\mathrm{slope}} = 44.6 \times 0.857/(2\ln 2) \approx 27.6$$，若强行解释为热稳定因子将给出$$E_b \approx 715\,\mathrm{meV}$$、$$\tau_{\mathrm{ret}}^{\mathrm{slope}} = \tau_0 e^{27.6} \approx 10^{12}\,\mathrm{ns} \approx 17\,\mathrm{min}$$，相对实际器件的纳秒响应偏离十二个数量级，物理上不可接受。

两路径之比$$\Delta_{\mathrm{slope}}/\Delta_{\mathrm{pulse}} = \beta_s^{\mathrm{meas}}/\beta_s^{\mathrm{NB}} \equiv \eta_c$$即为前文定义的C2C收窄因子，所反映的是NB Gumbel分布相对实验分布的过度展宽幅度。该恒等式给出$$\eta_c$$的另一组等价物理解释：同一组实验数据若分别按NB标度律与NB分布形态解读，得到的$$\Delta$$估计值差$$\eta_c$$倍。这一定量关系使NB模型表观的内部不一致问题转化为可控的双层分解。脉宽法$$\Delta = \Delta_{\mathrm{pulse}}$$描述热激活势垒并保留全部物理内涵，C2C分布形态相对Gumbel基线的偏差则单独通过乘性因子$$\eta_c$$描述。两者在数学上正交，在工程上各司其职：$$\Delta_{\mathrm{pulse}}$$作为热稳定参数进入工艺容差与可靠性分析，$$\eta_c$$作为分布形态参数进入Bernoulli采样精度估计。

由此澄清表观的Δ歧义。脉宽法$$\Delta$$是器件唯一具有热力学定义的稳定因子；斜率倒推所得的$$\Delta_{\mathrm{slope}}$$是数学上自洽但物理上空虚的拟合中间量，不应作为独立物理参数报告。本节后续工艺容差分析与下节联合写入概率模型均以$$\Delta_{\mathrm{pulse}}$$与$$\eta_c$$为两个独立的标定输入；混淆二者将导致对工艺余量、保持时间与热可靠性的根本性误判。

$$\eta_c$$的物理意义需结合两类独立随机性机制进行解释，二者对概率曲线形态的影响方向截然相反。器件间离散(D2D，device-to-device)源于工艺波动，体现为不同器件的$$(\Delta, V_{c0})$$存在晶圆级统计分布。由Jensen不等式及NB双指数函数的凸性可严格证明晶圆平均曲线的等效斜率$$\beta_{\mathrm{eff}} \leq \beta_s^{\mathrm{NB}}$$，D2D离散只能使阵列平均曲线变缓而不会变陡。循环间离散(C2C，cycle-to-cycle)是同一器件重复写入时由热涨落与随机初态引入的固有随机性，直接决定单器件$$P_{\mathrm{sw}}(V)$$曲线的陡度。纯NB单畴动力学对C2C的预测呈Gumbel型分布，在逻辑斯蒂过渡区具有确定的斜率$$\beta_s^{\mathrm{NB}}$$。本节所得$$\eta_c > 1$$是同一器件重复写入测出的斜率，反映的是C2C分布形态相对NB Gumbel的系统性收窄，不可能归因于D2D展宽。D2D与C2C的区分在后节进一步通过工艺容差仿真明确体现。

**$$\eta_c > 1$$的物理来源。** NB单畴Gumbel分布对C2C的过度展宽以及实测$$\eta_c$$随器件与方向变化的散布，反映多种微观机制的综合效应。Néel-Brown模型隐含严格单畴翻转假设，真实sMTJ在临界区可能经历亚畴协同跃迁，后者在电压维度的翻转分布较Gumbel更集中。$$t_w = 0.75\,\mathrm{ns}$$已进入热辅助进动翻转过渡区，介于纯热激活 (长脉宽、低电压) 与纯进动翻转 (短脉宽、高电压) 之间的机制交叉区，翻转事件获得部分相干性从而使C2C分布窄于纯热激活Gumbel极限。该过渡机制在多种SOT-MTJ单次翻转动力学表征中被一致报告，是在纳秒量级脉宽下最可能的主导机制。AP→P方向同时出现回跳平台与两段过渡而P→AP方向基本不出现的事实，提示参考层杂散场对正反两方向翻转动力学的非对称调制。

从工程角度看$$\eta_c \approx 5\text{-}10$$的C2C收窄对概率计算应用构成有利特性。相同电压噪声水平下更大的$$\beta_s$$对应更精确的概率编码，作为硬件Bernoulli采样单元时编码分辨率更高。但C2C分布形态依赖于具体器件设计与操作条件，对新器件结构须以实测Sigmoid为准，不可简单套用NB理论值。

### 2.3.5工艺波动对晶圆平均概率曲线的影响

前节确定了同一器件内部C2C分布相对NB Gumbel基线的收窄因子$$\eta_c$$，该量描述单器件层面的概率响应陡度。概率计算阵列由多个独立的Bernoulli采样单元构成，单器件实测特性能否在阵列层面统计地保持，取决于器件间(D2D)参数离散对晶圆平均概率响应的展宽幅度。从工程视角看，热稳定因子的相对涨落$$\mathrm{CV}(\Delta) \equiv \sigma_\Delta/\mu_\Delta$$是连接foundry失配数据与阵列级随机计算精度的关键传递参数。在固定写入电压下，$$\Delta$$的D2D离散直接转化为阵列内Bernoulli概率的器件间偏离，进而以二阶矩形式进入网络等效噪声方差，影响推断精度上界。

热稳定因子由几何体积、界面各向异性与饱和磁化共同决定，将$$\mathrm{CV}(\Delta)$$按物理来源分解便于明确工艺优化的优先级：若主要方差贡献来自几何刻蚀，则提升线宽控制的边际收益超过磁性材料优化；反之亦然。该分解不依赖于具体器件结构的微观参数细节，仅利用foundry PDK中直接可得的局部失配参数与隧穿电学模型即可完成，因此具有跨平台的迁移性。下文从sMTJ PDK中的电学失配出发，经Brinkman隧穿模型反推为几何与界面失配，进而合成$$\mathrm{CV}(\Delta)$$的物理基线，并定量评估该基线对晶圆平均概率响应的影响幅度。

**基于PDK失配参数的$$\mathrm{CV}(\Delta)$$反推。** 驰拓300 mm sMTJ PDK针对局部失配定义三个独立的Gaussian扰动参数，节选自器件模型process block如下：

```spectre
parameters Rsot_mis = 1
parameters Rp_mis   = 1
parameters TMR_mis  = 1
statistics {
    process {
        vary Rsot_mis dist = gauss std = 0.07
        vary Rp_mis   dist = gauss std = 0.07
        vary TMR_mis  dist = gauss std = 0.04
    }
}
```

扰动以乘法方式作用于标称值，std值即对应的相对标准差，给出$$\mathrm{CV}(R_P) = 7\%$$、$$\mathrm{CV}(R_{\mathrm{SOT}}) = 7\%$$、$$\mathrm{CV}(\mathrm{TMR}) = 4\%$$。

平行态MTJ电阻可写为$$R_P = (R\cdot A)/A$$，其中$$R\cdot A$$为隧穿电阻面积积 (仅取决于势垒特性) 、$$A = \pi D^2/4$$为MTJ面积 (仅取决于几何)。两类涨落相互独立，故$$\mathrm{CV}^2(R_P) = \mathrm{CV}^2(R\cdot A) + [2\mathrm{CV}(D)]^2$$。Brinkman低偏压近似[^ref-brinkman-bdr]给出$$\ln(R\cdot A) = \kappa t_{\mathrm{ox}}\sqrt{\bar\varphi} + \mathrm{const.}$$，其中$$\kappa = 2\sqrt{2m^*}/\hbar \approx 10.25\,\mathrm{nm^{-1}\cdot eV^{-1/2}}$$。对CoFeB/MgO接口取标称值$$t_{\mathrm{ox}} \approx 1\,\mathrm{nm}$$、$$\bar\varphi \approx 0.6\,\mathrm{eV}$$，Brinkman灵敏度系数为$$\partial\ln(R\cdot A)/\partial t_{\mathrm{ox}} = \kappa\sqrt{\bar\varphi} \approx 7.94\,\mathrm{nm^{-1}}$$与$$\partial\ln(R\cdot A)/\partial\bar\varphi = \kappa t_{\mathrm{ox}}/(2\sqrt{\bar\varphi}) \approx 6.62\,\mathrm{eV^{-1}}$$。取局部失配尺度下$$\sigma_{t_{\mathrm{ox}}}/t_{\mathrm{ox}} \approx 0.3\%$$、$$\sigma_{\bar\varphi}/\bar\varphi \approx 0.5\%$$得$$\mathrm{CV}(R\cdot A) \approx 3.1\%$$。从$$\mathrm{CV}(R_P) = 7\%$$中扣除势垒贡献后得几何贡献$$\mathrm{CV}(A) \approx 6.3\%$$，对应$$\mathrm{CV}(D) \approx 3.1\%$$，即$$D = 80\,\mathrm{nm}$$时$$\sigma_D \approx 2.5\,\mathrm{nm}$$，与300 mm高产平台双重曝光与精细刻蚀的线宽控制能力一致。

自由层磁性体积$$V_{\mathrm{mag}} = \pi D^2 t_f/4$$的相对涨落为$$\mathrm{CV}(V_{\mathrm{mag}}) = \sqrt{[2\mathrm{CV}(D)]^2 + \mathrm{CV}^2(t_f)}$$；取$$\mathrm{CV}(t_f) \approx 0.3\%$$ (MBE沉积厚度精度) 得$$\mathrm{CV}(V_{\mathrm{mag}}) \approx 6.3\%$$，几乎完全由$$\mathrm{CV}(D)$$主导而$$t_f$$贡献可忽略。界面各向异性场$$H_k$$与MgO/CoFeB界面质量直接相关，PDK将该界面效应投影到TMR失配上，故取$$\mathrm{CV}(H_k) \approx \mathrm{CV}(\mathrm{TMR}) = 4\%$$作为代理。饱和磁化强度$$M_s$$主要由CoFeB成分决定，取典型值$$\mathrm{CV}(M_s) \approx 2\%$$。

**方差合成与变异系数的几何加和规则。** 由$$\Delta \propto H_k M_s V_{\mathrm{mag}}$$对小相对扰动展开，$$\delta\Delta/\Delta = \delta H_k/H_k + \delta M_s/M_s + \delta V_{\mathrm{mag}}/V_{\mathrm{mag}} + \mathcal{O}(\mathrm{CV}^2)$$。诸误差项相互独立，方差线性可加而标准差按平方和的平方根合成：

$$
\mathrm{CV}^2(\Delta) = \mathrm{CV}^2(H_k) + \mathrm{CV}^2(M_s) + \mathrm{CV}^2(V_{\mathrm{mag}}),
$$

$$
\mathrm{CV}(\Delta) = \sqrt{0.040^2 + 0.020^2 + 0.063^2} = \sqrt{0.00598} \approx 7.7\%.
$$

各物理源对$$\mathrm{CV}^2(\Delta)$$的份额见表2.10。几何体积扰动$$V_{\mathrm{mag}}$$贡献66%总方差，界面各向异性$$H_k$$贡献27%，饱和磁化$$M_s$$贡献7%。该分解直接指向工艺优化优先级：将$$\mathrm{CV}(D)$$从3.1%进一步压缩至2%即可使$$V_{\mathrm{mag}}$$贡献的方差减半，对应$$\mathrm{CV}(\Delta)$$从7.7%降至约5.5%；同等比例改善$$M_s$$或$$H_k$$的影响约为前者的1/4–1/2。

**表2.10** $$\mathrm{CV}(\Delta)$$方差预算分解 (独立来源平方和合成) 

| 贡献源 | 相对标准差CV | $$\mathrm{CV}^2$$贡献 | 占总方差比例 |
|:---|:---:|:---:|:---:|
| 几何体积$$V_{\mathrm{mag}}$$ (由$$\mathrm{CV}(R_P)$$经Brinkman模型分解，含直径与厚度合成扰动) | 6.3% | 0.0040 | 66% |
| 界面各向异性$$H_k$$ (由$$\mathrm{CV}(\mathrm{TMR})$$作为界面质量代理) | 4.0% | 0.0016 | 27% |
| 饱和磁化$$M_s$$ (CoFeB成分涨落，取自文献典型值) | 2.0% | 0.0004 | 7% |
| 平方和合成$$\mathrm{CV}(\Delta) = \sqrt{\sum\mathrm{CV}_i^2}$$ | 7.7% | 0.0060 | 100% |

注：$$V_{\mathrm{mag}}$$的6.3%已由$$\mathrm{CV}(V_{\mathrm{mag}}) = \sqrt{[2\mathrm{CV}(D)]^2 + \mathrm{CV}^2(t_f)}$$合成，几何上由$$\mathrm{CV}(D) \approx 3.1\%$$主导，厚度涨落$$\mathrm{CV}(t_f) \approx 0.3\%$$对$$\mathrm{CV}^2(V_{\mathrm{mag}})$$的相对贡献仅约0.2%，本预算中不再单列。

**翻转概率对$$\Delta$$扰动的解析灵敏度。** 工艺波动对晶圆平均概率响应的影响幅度并不直接由$$\mathrm{CV}(\Delta)$$本身决定，而由$$\mathrm{CV}(\Delta)$$与$$P_{\mathrm{sw}}$$对$$\Delta$$灵敏度的乘积控制。NB双指数函数形式使该灵敏度具有简洁的解析表达，可在Monte Carlo仿真之前给出工艺裕度的解析估计。

记$$f(V) = (t_w/\tau_0)\exp[-\Delta(1-V/V_{c0})]$$，则$$P_{\mathrm{sw}} = 1 - \exp(-f)$$。对$$\Delta$$求偏导得

$$
\frac{\partial P_{\mathrm{sw}}}{\partial \Delta} = (1 - P_{\mathrm{sw}})\, f \cdot \left(\frac{V}{V_{c0}} - 1\right)
$$

该表达式在NB$$P_{\mathrm{sw}} = 0.5$$点取得简洁形式：此处$$\exp(-f) = 1/2$$给出$$f = \ln 2$$、$$1 - P_{\mathrm{sw}} = 0.5$$，并由NB阈值条件$$V_{\mathrm{th}}/V_{c0} = 1 - \ln(t_w/(\tau_0 \ln 2))/\Delta$$得

$$
\left.\frac{\partial P_{\mathrm{sw}}}{\partial \Delta}\right|_{V_{\mathrm{th}}} = -\frac{\ln 2}{2}\cdot\frac{1}{\Delta}\ln\!\left(\frac{t_w}{\tau_0 \ln 2}\right)
$$

灵敏度仅取决于无量纲组合$$\xi(t_w) \equiv \ln(t_w/(\tau_0\ln 2))/\Delta$$，即工作脉宽相对零驱动保持时间$$\tau_{\mathrm{ret}} = \tau_0 e^\Delta$$的对数距离。极限行为有清晰的物理对应。$$t_w \to \tau_0\ln 2$$时$$\xi \to 0$$，灵敏度趋零，此即确定性翻转极限$$V_{\mathrm{th}} \to V_{c0}$$。$$t_w \to \tau_{\mathrm{ret}}$$时$$\xi \to 1$$，灵敏度趋$$\ln 2/2 \approx 0.35$$，此即完全热激活极限。

代入Device A、P→AP、$$t_w = 0.75\,\mathrm{ns}$$、$$\tau_0 = 1\,\mathrm{ns}$$、$$\Delta = 4.91$$：$$\ln(0.75/0.693) \approx 0.0792$$，$$\xi = 0.0792/4.91 \approx 0.0161$$，故

$$
\left.\frac{\partial P_{\mathrm{sw}}}{\partial \Delta}\right|_{V_{\mathrm{th}}} \approx -5.6 \times 10^{-3}
$$

由该灵敏度可解析估计单器件层面工艺扰动引入的$$P_{\mathrm{sw}}$$方差：$$\sigma_{P_{\mathrm{sw}}}|_{V_{\mathrm{th}}} \approx |\partial P_{\mathrm{sw}}/\partial \Delta| \cdot \sigma_\Delta = 5.6 \times 10^{-3} \cdot \mathrm{CV}_\Delta \cdot \Delta$$。在PDK基线$$\mathrm{CV}_\Delta = 7.7\%$$下，单器件$$P_{\mathrm{sw}}$$扰动幅度仅0.21%；即使在极端工艺条件$$\mathrm{CV}_\Delta = 60\%$$下亦不超过1.7%，与图2.15(b)的插图给出的晶圆平均偏离幅度在数值上同阶 (晶圆平均还含Sigmoid曲率效应的二阶贡献，故略大)。

上述灵敏度的极小是器件物理工作点选择的直接结果：亚纳秒脉宽下sMTJ被驱动至接近确定性翻转区($$V_{\mathrm{th}}/V_{c0} \approx 0.984$$)，此时$$P_{\mathrm{sw}}$$对势垒高度的微调几乎不敏感，写入概率主要由电压相对$$V_{c0}$$的位置决定。该特性是低势垒sMTJ作为概率单元的核心优势，以纳秒级脉宽换取阵列级概率一致性，而工艺容差被自动压缩。

将该灵敏度估计与传统存储MRAM对比，若以保持时间$$t_w \sim \tau_{\mathrm{ret}} \sim 10^{17}\,\mathrm{ns}$$ (10年) 作为评估点，则$$\xi \to 1$$、$$\partial P_{\mathrm{sw}}/\partial \Delta \to -\ln 2/2$$，此时$$\mathrm{CV}(\Delta) = 7.7\%$$即可使保持失败概率发生$$\sim 0.35 \times 0.077 \times \Delta = \mathcal{O}(1)$$量级变化，工艺容差极为苛刻。本节器件位于另一极端的纳秒激活区，对$$\Delta$$扰动具有内在低灵敏度，工艺余量天然宽裕。

**Monte Carlo数值验证与工艺裕度。** 对$$\Delta \sim \mathcal{N}(\mu_\Delta, \sigma_\Delta^2)$$的器件集合，晶圆级平均概率曲线为
$$
\bar{P}_{\mathrm{sw}}(V, t_w) = \int P_{\mathrm{sw}}(V, t_w \mid \Delta)\, f_\Delta(\Delta)\,\mathrm{d}\Delta
$$

对该平均曲线重新进行Sigmoid拟合得等效斜率$$\beta_{\mathrm{eff}}$$。定义D2D传递函数$$\mathcal{F}(\mathrm{CV}_\Delta) \equiv \beta_{\mathrm{eff}}/\beta_{\mathrm{NB}}^{\mathrm{fit}}$$，其中$$\beta_{\mathrm{NB}}^{\mathrm{fit}}$$为$$\mathrm{CV}_\Delta = 0$$极限下NB单器件曲线用logistic函数拟合得到的斜率(Device A、P→AP、$$t_w = 0.75\,\mathrm{ns}$$条件下约8.35 V$$^{-1}$$)。对$$\mathrm{CV}_\Delta \in \{0, 3\%, \ldots, 60\%\}$$每个值生成$$N = 2 \times 10^4$$个Gaussian样本计算晶圆平均曲线并拟合，得$$\mathcal{F}(\mathrm{CV}_\Delta)$$的数值函数。

由Jensen不等式与NB双指数函数对$$\Delta$$的凸性可严格证明$$\mathcal{F}(\mathrm{CV}_\Delta) \leq 1$$对所有$$\mathrm{CV}_\Delta \geq 0$$成立。D2D离散单调展宽阵列平均曲线，与前节解析灵敏度$$\partial P_{\mathrm{sw}}/\partial \Delta < 0$$给出的方向一致。考虑前节确定的C2C收窄因子$$\eta_c = 5.34$$[^note-eta-fit]后，阵列级Sigmoid斜率的联合预测为

$$
\beta^{\mathrm{eff}}(\mathrm{CV}_\Delta) = \eta_c \cdot \mathcal{F}(\mathrm{CV}_\Delta) \cdot \beta_{\mathrm{NB}}^{\mathrm{fit}}
$$

仿真结果如图2.15所示。在PDK基线$$\mathrm{CV}_\Delta = 7.7\%$$处$$\mathcal{F} = 0.997$$，联合预测$$\beta^{\mathrm{eff}} = 44.5\,\mathrm{V^{-1}}$$，为实测$$\beta_s = 44.6\,\mathrm{V^{-1}}$$的99.7%。该0.3%的修正幅度可由前节解析灵敏度估计独立验证：在PDK基线下单器件$$\sigma_{P_{\mathrm{sw}}} \approx 0.21\%$$，按晶圆平均的Jensen修正，斜率退化预期约$$(\sigma_{P_{\mathrm{sw}}}/0.5)^2 \cdot \mathcal{O}(1) \sim 0.2\%$$，与MC数值结果同阶。三层证据 (PDK-Brinkman反推$$\mathrm{CV}_\Delta$$、解析灵敏度估计、Monte Carlo数值仿真) 相互支撑，确认双层分解框架的定量自洽性。

工艺裕度评估见表2.11。工艺容差保留度随$$\mathrm{CV}_\Delta$$增大缓慢下降，$$\geq 99\%$$、$$\geq 95\%$$、$$\geq 90\%$$三档保留度目标对应的容差边界分别为15.4%、36.5%、58.6%。PDK基线7.7%已位于强保留度区(99.7%)内部，工艺余量充裕。该宽容差正是前节灵敏度分析揭示的物理后果：纳秒激活工作点$$\xi(0.75\,\mathrm{ns}) = 0.016$$将$$P_{\mathrm{sw}}$$对$$\Delta$$扰动的响应压制到亚百分之级，即便工艺$$\mathrm{CV}_\Delta$$扩大到接近一倍($$\sim 60\%$$)，阵列平均斜率退化仍在10%以内。

**表2.11** 工艺裕度：阵列$$\beta^{\mathrm{eff}}$$相对单器件$$\beta_s$$的保持度与所需$$\mathrm{CV}_\Delta$$上限(Device A、P→AP、$$t_w = 0.75\,\mathrm{ns}$$)

| $$\beta^{\mathrm{eff}}/\beta_s^{\mathrm{meas}}$$目标 | 所需$$\mathrm{CV}_\Delta$$上限 |
|:---:|:---:|
| $$\geq 99\%$$ | 15.4% |
| $$\geq 95\%$$ | 36.5% |
| $$\geq 90\%$$ | 58.6% |

将该结果置于工程语境，在仅将工艺失配投影为$$\Delta$$单参数扰动的NB传递函数中，PDK基线对应单器件$$\beta_s$$保持率99.7%，远高于多数概率计算应用的$$\geq 95\%$$精度门槛。这说明热稳定因子离散本身对该纳秒工作点的晶圆平均斜率退化较弱；若进一步讨论固定端电压写入，则还需把$$R_{\mathrm{SOT}}$$引起的电压到电流换算漂移以及磁性参数对动力学阈值的共同影响纳入器件级求解。故工艺优化不能仅凭$$\Delta$$单参数图判断实测Sigmoid局部形状，仍需结合下述宏自旋失配仿真与C2C分布形态的器件级稳定性 (即$$\eta_c$$的器件间均匀性) 共同评估。

需要指出，本节基于PDK标称失配参数推导得到的$$P_{\mathrm{sw}}$$响应是一条统计意义下的等效Sigmoid曲线，其形状由$$(\mathrm{CV}_\Delta, \eta_c)$$两个标量参数决定，无法刻画实测中观察到的back-hopping平台、两段过渡等单器件级畸变。换言之，实测$$P_{\mathrm{sw}}$$相对PDK推导基线存在的局部漂变要大于PDK单参数化所给出的展宽幅度。完整的非理想性建模——包括back-hopping、亚畴协同跃迁、尝试频率电压依赖等微观机制——已在2.2节通过sLLG动力学求解器、温度依赖材料参数反馈与TMR/电输运非线性等模块完整实现，可在系统仿真前端逐器件注入实测中观察到的所有畸变特征，从而保证阵列级建模与实验的端到端对齐。



---

![工艺波动对sMTJ概率响应的综合影响](figs/Chapter02_local_15.png)

**图2.15** 工艺波动对sMTJ概率响应的综合影响，基于PDK失配的方差预算与Monte Carlo验证，以Device A、P→AP、$$t_w = 0.75\,\mathrm{ns}$$实测为基准。(a)$$\mathrm{CV}(\Delta) = 7.7\%$$方差预算分解，$$V_{\mathrm{mag}}$$中由横向面积失配引入的份额贡献约66%方差、$$H_k$$约27%、$$M_s$$约7%，自由层厚度$$t_f$$单独贡献约0.2%；红色虚线标示按平方和合成法则得到的总$$\mathrm{CV}(\Delta) = 7.7\%$$位置[^note-variance-sum]。(b)不同$$\mathrm{CV}_\Delta$$下的C2C校准晶圆平均Sigmoid曲线族，CV=0按构造等于实测 (青色虚线)，PDK基线$$\mathrm{CV}_\Delta = 7.7\%$$ (琥珀色) 与实测几乎完全重合，CV扩展至60%以体现极端工艺条件下的微弱展宽；插图给出各曲线相对CV=0的偏差$$\Delta P_{\mathrm{sw}}$$ (单位%)，呈现典型的双叶结构 (过渡区前后符号相反，对应Sigmoid斜率减缓)，PDK基线偏差<0.3%、CV=60%偏差达约$$\pm 3\%$$，与解析灵敏度$$\partial P_{\mathrm{sw}}/\partial \Delta \approx -5.6 \times 10^{-3}$$的预测一致。(c)D2D传递函数$$\mathcal{F}(\mathrm{CV}_\Delta) = \beta_{\mathrm{eff}}/\beta_{\mathrm{NB}}^{\mathrm{fit}}$$的Monte Carlo数值，PDK基线处$$\mathcal{F} = 0.997$$ (琥珀色星号)，$$\mathrm{CV}_\Delta = 60\%$$时降至约0.90，单调下降反映Jensen不等式的渐进生效。(d)四组参考的联合对比 (对数y轴)，蓝色虚线为NB单器件拟合斜率 (约8.35 V$$^{-1}$$) 、蓝色方块为晶圆NB预测$$\mathcal{F}\cdot\beta_{\mathrm{NB}}^{\mathrm{fit}}$$、青色虚线为实测$$\beta_s = 44.6\,\mathrm{V^{-1}}$$、红色圆线为联合预测$$\eta_c\mathcal{F}\beta_{\mathrm{NB}}^{\mathrm{fit}}$$；PDK基线处联合预测44.5 V$$^{-1}$$为实测99.7%，验证双层分解框架的定量自洽性。曲线在面板内挤压程度小这一点本身即是物理结果：$$t_w = 0.75\,\mathrm{ns}$$工作点($$V_{\mathrm{th}}/V_{c0} \approx 0.984$$)已接近NB确定性极限，对$$\Delta$$扰动的灵敏度天然较低，因此即便$$\mathrm{CV}_\Delta$$高至60%，阵列平均斜率退化也仅10%量级。

为检验上述$$\Delta$$单参数投影是否遗漏端电压驱动下的器件级展宽，本工作进一步将同一PDK失配预算注入宏自旋求解器：由$$R_P$$残余失配采样$$D$$与$$D_{\mathrm{elec}}$$，由$$t_f$$、$$M_s$$与TMR代理采样自由层厚度、饱和磁化和各向异性参数，并将$$R_{\mathrm{SOT}}$$失配映射到$$R_W$$。该设置在固定$$|V_{\mathrm{SOT}}|$$下尤其关键，因为$$R_W$$变化会同步改变每个样本实际获得的$$|I_{\mathrm{SOT}}|$$与SOT自热功率。图2.16显示标称器件在实测阈值附近仍保留清晰的Sigmoid上升段，而工艺失配平均曲线因D2D阈值横向散布而变浅。由此可见，图2.15给出的“PDK基线下$$\Delta$$扰动展宽较弱”结论成立于NB单参数传递函数；若观察固定端电压下的宏自旋响应，则磁性参数漂移与$$R_{\mathrm{SOT}}$$电压-电流换算漂移仍会对晶圆平均$$P_{\mathrm{sw}}$$曲线产生可见展宽。

---

![宏自旋工艺失配Monte Carlo下的Psw端电压扫描](figs/Chapter02_local_16.png)

**图2.16** 工艺失配下的端电压写入概率Monte Carlo曲线，$$t_p = 0.75\,\mathrm{ns}$$、$$V_{\mathrm{MTJ}} = 0$$。主图将图2.10的$$|I_{\mathrm{SOT}}| = 300\text{-}3500\,\mu\mathrm{A}$$宽谱扫描按标称$$R_W$$映射为固定$$|V_{\mathrm{SOT}}|$$轴，插图放大$$800\text{-}1400\,\mu\mathrm{A}$$阈值窗口对应的电压区间。蓝色圆线为标称宏自旋曲线，红色方线为PDK失配样本的晶圆平均，阴影为Wilson 95%区间，绿色虚线为同批次实测$$V_{\mathrm{th}} = 894\,\mathrm{mV}$$参考线；工艺样本采用6个D2D失配器件、每点8次热噪声试验。插图中标称曲线在阈值附近仍呈Sigmoid上升，失配平均曲线变浅则反映不同器件阈值横向散布以及$$R_{\mathrm{SOT}}$$引起的端电压到驱动电流换算漂移，而非单器件Sigmoid特征消失。

### 2.3.6采样数对概率估计精度的影响

2.3.5至2.3.5节建立的联合写入概率模型在两种尺度上都依赖于采样操作。第一种是2.3.5节用于将PDK失配参数向晶圆平均Sigmoid响应映射的Monte Carlo仿真，每个$$\mathrm{CV}_\Delta$$值下对$$\Delta\sim\mathcal{N}(\mu_\Delta, \sigma_\Delta^2)$$抽取$$N$$个样本计算晶圆平均曲线并拟合$$\beta_{\mathrm{eff}}$$，由此得到D2D传递函数$$\mathcal{F}(\mathrm{CV}_\Delta)$$，此为建模端采样。第二种是sMTJ作为硬件Bernoulli随机源在运行时的物理采样，固定偏置$$(V, t_w)$$下重复执行$$K$$次写入-读取循环，以经验频率$$\hat p_K$$估计$$p = P_{\mathrm{sw}}(V, t_w)$$，此为硬件端采样。两者虽同为采样但物理含义与优化目标截然不同：建模端的$$N$$取决于所需建模精度与仿真预算，与器件运行无关；硬件端的$$K$$取决于每概率输出的能耗-精度预算，每次采样消耗0.78 pJ能量并占用0.75 ns时间。本节分两小节依次给出两种采样数的量化分析。2.3.7.1节以PDK工艺工况为映射锚点，通过Monte Carlo多次独立重复给出$$\hat{\mathcal{F}}(N,\mathrm{CV}_\Delta)$$的偏差-方差预算与推荐$$N$$值；2.3.7.2节对运行时采样则利用Binomial分布可解析的特性直接以精确覆盖率表达式求解最小$$K$$，同时以Monte Carlo路径与CLT近似作为对照验证精确解的必要性。

**建模端：$$\hat{\mathcal{F}}(N, \mathrm{CV}_\Delta)$$估计量对MC采样数$$N$$的敏感度。** 2.3.5节以$$N = 20000$$次Gaussian采样作为建模基线，确保$$\mathcal{F}(\mathrm{CV}_\Delta)$$的一次性标定精度。但若将该双层框架嵌入电路级行为仿真、PBNN训练前端的器件mismatch注入、或硬件在环工艺校准回路，$$\mathcal{F}(\cdot)$$可能需反复求值成百上千次，此时$$N$$直接决定总体仿真时间。以此视角，$$N$$越小越好，但需量化有限$$N$$引入的估计方差是否可接受，并给出与PDK工艺工况匹配的最小采样数建议。

固定操作点为主基准Device A, P→AP, $$t_w = 0.75\,\mathrm{ns}$$；选取三档代表性PDK工况，即$$\mathrm{CV}_\Delta = 7.7\%$$(PDK基线，2.3.5节经Brinkman反推得到的实际工艺水平)、15%(中等工艺恶化，对应2.3.5节表2.11中$$\beta$$保持度≥99%的边界)、30%(严重工艺恶化，对应$$\beta$$保持度约96%)；在$$N\in\{100, 200, 500, 1000, 2000, 5000, 10000\}$$对数网格上以$$R = 40$$个独立随机种子重复执行MC估计，记录每次估计值

$$
\hat{\mathcal{F}}(N,\mathrm{CV}_\Delta) = \frac{\hat\beta_{\mathrm{eff}}(N, \mathrm{CV}_\Delta)}{\beta_{\mathrm{NB}}^{\mathrm{fit}}},\qquad \hat\beta_{\mathrm{eff}} = \mathrm{Sigmoid\,fit\,of}\ \frac{1}{N}\sum_{i=1}^N P_{\mathrm{sw}}(V\mid\Delta_i),\ \Delta_i\sim\mathcal{N}(\mu_\Delta, \sigma_\Delta^2)
$$

以$$N_{\mathrm{ref}} = 50000$$的单次估计作为参考真值$$\mathcal{F}_{\mathrm{ref}}$$，记录$$\hat{\mathcal{F}}$$在不同$$N$$下的均值(偏差指标)与种子间标准差(方差指标)。晶圆平均曲线$$\bar P_{\mathrm{sw}}(V) = N^{-1}\sum_i P_{\mathrm{sw}}(V\mid\Delta_i)$$在每个电压点$$V$$上是$$N$$个独立随机变量的算术平均，由中心极限定理其标准误按$$N^{-1/2}$$衰减；Sigmoid拟合将该逐点噪声映射到斜率参数$$\hat\beta_{\mathrm{eff}}$$，继承相同的$$N^{-1/2}$$标度。采样数敏感性仿真结果汇总于表2.12、图示于图2.17。

**表2.12** 不同$$(\mathrm{CV}_\Delta, N)$$组合下的估计量性能($$R = 40$$次独立种子，$$\mathcal{F}_{\mathrm{ref}}$$由$$N_{\mathrm{ref}} = 50000$$给出)

| $$\mathrm{CV}_\Delta$$ | $$N$$ | $$\mathcal{F}_{\mathrm{ref}}$$ | $$\langle\hat{\mathcal{F}}\rangle$$ | 偏差 | $$\sigma(\hat{\mathcal{F}})$$ | $$\sigma/\langle\hat{\mathcal{F}}\rangle$$ | $$P[\text{rel.err.}<2\%]$$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 7.7% (PDK基线) | 100 | 0.9973 | 0.9983 | +0.10% | 0.0080 | 0.81% | 100.0% |
| 7.7% | 1000 | 0.9973 | 0.9971 | −0.02% | 0.0024 | 0.24% | 100.0% |
| 7.7% | 10000 | 0.9973 | 0.9977 | +0.04% | 0.0009 | 0.09% | 100.0% |
| 15% (中等) | 100 | 0.9901 | 0.9923 | +0.22% | 0.0170 | 1.71% | 80.0% |
| 15% | 200 | 0.9901 | 0.9904 | +0.03% | 0.0117 | 1.18% | 87.5% |
| 15% | 10000 | 0.9901 | 0.9906 | +0.05% | 0.0013 | 0.13% | 100.0% |
| 30% (严重) | 100 | 0.9629 | 0.9549 | −0.82% | 0.0328 | 3.43% | 30.0% |
| 30% | 2000 | 0.9629 | 0.9635 | +0.06% | 0.0066 | 0.69% | 100.0% |
| 30% | 10000 | 0.9629 | 0.9633 | +0.04% | 0.0028 | 0.28% | 100.0% |

关键观察有三项。其一，$$\sigma(\hat{\mathcal{F}})$$在对数-对数坐标下与$$\propto N^{-1/2}$$参考虚线严格重合(图2.17(b))，偏差$$|\langle\hat{\mathcal{F}}\rangle - \mathcal{F}_{\mathrm{ref}}|$$在所有$$(N, \mathrm{CV}_\Delta)$$组合下均小于1%，证实$$\hat{\mathcal{F}}$$为近似无偏估计量。其二，相对标准差$$\sigma(\hat{\mathcal{F}})/\langle\hat{\mathcal{F}}\rangle$$在固定$$N$$下随$$\mathrm{CV}_\Delta$$近似线性增长(图2.17(c))，符合总方差$$\mathrm{Var}(\hat{\mathcal{F}})\propto\mathrm{CV}_\Delta^2/N$$的理论预期；PDK基线工况在$$N\geq 100$$即可持续低于1%相对误差水平。其三，$$N_{\mathrm{rec}}$$的精确求解采用两阶段迭代算法：首先对主扫$$\sigma(N)$$数据做$$\sigma = C\cdot N^{-1/2}$$对数线性回归提取标度系数$$C$$，由覆盖率表达式$$2\Phi(\mathrm{tol}\cdot F_{\mathrm{ref}}/\sigma(N)) - 1 \geq 1-\alpha$$解析反推种子点$$N_{\mathrm{seed}} = (Cz/(\mathrm{tol}\cdot F_{\mathrm{ref}}))^2$$；其次在$$N_{\mathrm{seed}}$$附近以$$R_{\mathrm{verify}} = 150$$个独立种子执行持久单调二分搜索，直接给出整数精度下的最小采样数(图2.17(d))。

推荐$$N$$值在三档工况下近似按$$\mathrm{CV}_\Delta^2$$标度增长，以2%精度行为例57:215:1023约为1:3.8:18，与理论比1:3.8:15吻合；偏离反映Sigmoid拟合在过渡区的轻度非线性放大效应。该标度意味着PDK工艺所处的工艺水平直接决定合适的MC预算：量产工艺维持PDK基线水平时$$N_{\mathrm{rec}} = 57$$即可达2%精度，若工艺恶化至CV≤15%(仍在2.3.5节表2.11允许工艺裕度内)$$N$$需升至约215，仅当工艺显著退化(CV≥30%)时才需$$N\sim 10^3$$量级。相对于2.3.5节的$$N = 20000$$建模基线，PDK基线工况可压缩$$N$$达350倍，由于MC的主要计算量为$$N\times|V_{\mathrm{grid}}|$$规模的元素级指数与均值运算，所得加速比直接线性正比于$$N$$的压缩倍数，可直接转化为电路仿真的实际时间节约。若精度基准放宽至5%相对误差，全部三档工况$$N_{\mathrm{rec}}\leq 143$$[^note-mc-budget]。实用部署中建议在工艺监控中同步跟踪$$\mathrm{CV}_\Delta$$、按上述$$N_{\mathrm{rec}}$$表自适应调整嵌入式MC采样数；若工艺控制稳定在PDK基线附近则可将建模采样预算大幅压缩，为上层应用释放计算资源。

---

![MC采样数对D2D传递函数估计量的影响](figs/Chapter02_local_17.png)

**图2.17** D2D传递函数$$\hat{\mathcal{F}}(N,\mathrm{CV}_\Delta)$$估计量对Monte Carlo采样数$$N$$的敏感度($$R = 40$$独立种子，Device A, P→AP, $$t_w = 0.75\,\mathrm{ns}$$；$$\mathcal{F}_{\mathrm{ref}}$$由$$N_{\mathrm{ref}} = 50000$$给出)。(a)$$\hat{\mathcal{F}}$$的种子平均值与5至95%分位带随$$N$$的变化，三条水平点线为对应$$\mathcal{F}_{\mathrm{ref}}$$，PDK基线(琥珀)分位带最窄、严重工况(红)最宽但均围绕各自参考值收敛。(b)种子间标准差$$\sigma(\hat{\mathcal{F}})$$对数-对数图，三条实线与$$\propto N^{-1/2}$$参考虚线(黑)斜率相同，验证中心极限定理标度。(c)相对标准差$$\sigma(\hat{\mathcal{F}})/\langle\hat{\mathcal{F}}\rangle$$(百分比)映射到$$\beta^{\mathrm{eff}}$$预测的相对误差，两条横向点线标示1%与2%精度目标；PDK基线在$$N\geq 100$$即可持续低于1%水平。(d)不同精度容差(1%, 2%, 5%)下满足95%置信度的最小$$N$$柱状图，数值由$$N^{-1/2}$$标度拟合与持久单调二分搜索精确迭代得到，整数精度；PDK基线下2%精度仅需$$N = 57$$，5%精度全部三档工况均可压缩至$$N\leq 143$$。

**硬件端：sMTJ作为Bernoulli随机源的运行时采样数。** 硬件端的采样与建模端存在本质差异：每次采样对应器件的一次物理写入-读取循环，受限于能耗0.78 pJ/次、延迟0.75 ns/次与器件耐久性。在概率表征阶段$$K$$决定$$P_{\mathrm{sw}}$$曲线拟合的置信度(2.3.2节采用$$K = 100$$)，在随机比特流算术中$$K$$决定数值精度与比特流长度，在PBNN推断中$$K$$决定单次前向传播中每个节点的采样次数，$$K$$的压缩直接映射到每概率输出的能耗与吞吐率。与建模端不同，硬件采样问题存在可解析的精确解：$$\hat p_K = K^{-1}\sum_{i=1}^K X_i$$中的$$X_i\sim\mathrm{Ber}(p)$$独立同分布，$$K\hat p_K\sim\mathrm{Bin}(K, p)$$，覆盖率可由Binomial CDF直接求出。下面并列对比三条分析路径(精确Binomial、Monte Carlo、CLT近似)以明确各自的有效性边界。

对任意容差$$\varepsilon > 0$$与置信度$$1-\alpha$$，Binomial精确覆盖率为

$$
\mathrm{cov}(K, p, \varepsilon) \equiv P\!\left(|\hat p_K - p| < \varepsilon\right) = F_{\mathrm{Bin}(K,p)}(k_{\mathrm{hi}}) - F_{\mathrm{Bin}(K,p)}(k_{\mathrm{lo}} - 1)
$$

其中$$k_{\mathrm{lo}} = \lfloor K(p-\varepsilon)\rfloor + 1$$、$$k_{\mathrm{hi}} = \lceil K(p+\varepsilon)\rceil - 1$$为严格不等式$$|\hat p_K - p| < \varepsilon$$对应的整数界，此路径无任何抽样噪声。Monte Carlo路径则对$$(K, p, \varepsilon)$$执行$$M$$次独立Bernoulli(K, p)采样得$$\hat{\mathrm{cov}}_{\mathrm{MC}}$$，本身为Binomial(M, cov)/M型估计量、RMSE按$$M^{-1/2}$$缩放。CLT近似路径给出$$\mathrm{cov} \approx 2\Phi(\varepsilon\sqrt{K/(p(1-p))}) - 1$$以及样本复杂度$$K_{\mathrm{CLT}} = z_{\alpha/2}^2 p(1-p)/\varepsilon^2$$，$$z_{\alpha/2} = 1.960$$对应95%置信度、最坏情况$$p = 0.5$$时$$K_{\mathrm{CLT}}\leq z^2/(4\varepsilon^2)$$。图2.18(a)通过$$K\in\{5, 20, 100, 500\}$$的Binomial PMF直观展示Sigmoid中段的离散性：$$K = 5$$时可达频率集合$$\{0, 0.2, 0.4, 0.6, 0.8, 1.0\}$$中无任何值落入$$\varepsilon = 0.05$$窗口$$(0.45, 0.55)$$、覆盖率为$$0\%$$，这是CLT连续近似完全失效的极端离散区。

CLT公式仅为渐近近似，在$$K$$较小的区间内Binomial覆盖率呈显著的离散阶梯结构：每当$$K$$跨过某个整数阈值使新的可达频率$$k/K$$进入或离开误差带$$(p-\varepsilon, p+\varepsilon)$$时覆盖率发生阶跃，阶梯局部可能出现短暂倒退(即$$K$$增大反而覆盖率下降1至2个百分点；图2.18(c)在$$\varepsilon = 0.02$$、$$K\in[200]$$范围锯齿振幅最大达$$\pm 10\%$$)。简单的网格搜索或单点阈值判定容易在此类局部涨落处给出偏小的$$K$$估计而在实际部署中无法持续满足精度要求。为获得对后续运行持续可靠的采样数下限，定义持久单调阈值

$$
K_{\mathrm{req}}(p, \varepsilon, \alpha) = \min\!\left\{K^*\in\mathbb{N}: \mathrm{cov}(K', p, \varepsilon)\geq 1-\alpha\quad\forall K'\geq K^*\right\}
$$

即从$$K^*$$起覆盖率始终不再跌破目标置信度$$1-\alpha$$的最小整数。该定义下$$K_{\mathrm{req}}$$是真正的安全下限；求解算法为两步迭代：首先以CLT估计值的四倍为初始上界$$K_{\mathrm{hi}}$$并倍增扩展直至稳定合格，其次对$$[1, K_{\mathrm{hi}}]$$逐点计算覆盖率、自右向左取后缀最小值$$\mathrm{cov}_{\min}(K) \equiv \min_{K\leq K'\leq K_{\mathrm{hi}}}\mathrm{cov}(K')$$，$$\{K: \mathrm{cov}_{\min}(K)\geq 1-\alpha\}$$的最小元即$$K_{\mathrm{req}}$$。基于SciPy `binom.cdf`的精确数值评估下典型工作点总运行时间毫秒级。同一算法可平行作用于$$\hat{\mathrm{cov}}_{\mathrm{MC}}$$得到MC版本$$\hat K_{\mathrm{req}}^{\mathrm{MC}}$$，后者由于MC估计量的随机涨落既可能高估也可能低估真实$$K_{\mathrm{req}}$$。

将上述框架应用于主基准Sigmoid($$V_{\mathrm{th}} = 894\,\mathrm{mV}$$，$$\beta_s = 44.6\,\mathrm{V}^{-1}$$)，选取代表性工作点$$p\in\{0.1, 0.5, 0.9\}$$；由$$V(p) = V_{\mathrm{th}} + \beta_s^{-1}\ln[p/(1-p)]$$得对应偏置电压分别为844.7 mV、894.0 mV、943.3 mV，覆盖Sigmoid的低、中、高概率区。精确、MC、CLT三种估计器的并列求解结果见表2.13。

**表2.13** 不同工作点$$p$$与精度$$\varepsilon$$组合下三种估计器的$$K_{\mathrm{req}}$$对比(95%置信度，MC以$$M = 2\times 10^4$$计算)

| 工作点$$p$$ | 精度$$\varepsilon$$ | $$K_{\mathrm{req}}^{\mathrm{exact}}$$ | $$K_{\mathrm{req}}^{\mathrm{MC}}$$ | $$K_{\mathrm{CLT}}$$ | MC偏差 | CLT偏差 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.1 | 0.02 | **892** | 916 | 864.3 | +2.7% | −3.1% |
| 0.1 | 0.05 | **147** | 147 | 138.3 | 0.0% | −5.9% |
| 0.1 | 0.10 | **41** | 41 | 34.6 | 0.0% | −15.6% |
| 0.5 | 0.02 | **2451** | 2472 | 2400.9 | +0.9% | −2.0% |
| 0.5 | 0.05 | **399** | 410 | 384.1 | +2.8% | −3.7% |
| **0.5** | **0.10** | **106** | 101 | 96.0 | **−4.7%** | −9.4% |
| 0.9 | 0.02 | **901** | 901 | 864.3 | 0.0% | −4.1% |
| 0.9 | 0.05 | **147** | 147 | 138.3 | 0.0% | −5.9% |
| 0.9 | 0.10 | **41** | 41 | 34.6 | 0.0% | −15.6% |

精确算法通过严格单调性核验：以$$(p=0.5, \varepsilon=0.02)$$为例，$$\mathrm{cov}(2450) = 0.9500$$、$$\mathrm{cov}(2451) = 0.9523$$，且2451之后覆盖率始终维持或超过目标。表2.13揭示下列几项规律。其一，精度与$$K_{\mathrm{req}}$$呈严格$$\varepsilon^{-2}$$标度。对固定$$p$$将$$\varepsilon$$从0.02放宽至0.10，$$K_{\mathrm{req}}$$由892(或2451)降至41(或106)，与理论比值25吻合至$$\pm 5\%$$；图2.18(d)以双对数坐标给出$$p = 0.5$$下$$K_{\mathrm{req}}(\varepsilon)$$从$$\varepsilon = 0.01$$延伸至$$\varepsilon = 0.20$$的完整关系，与$$\varepsilon^{-2}$$参考虚线吻合，为系统设计师提供简洁的容差-采样数换算工具。其二，$$p = 0.5$$对采样数最不友好。在同一$$\varepsilon$$下$$K_{\mathrm{req}}(0.5)$$约为$$K_{\mathrm{req}}(0.1)$$或$$K_{\mathrm{req}}(0.9)$$的2.6至2.7倍，与$$p(1-p)$$在$$p = 0.5$$达最大值0.25的理论预期(最坏比值$$0.25/0.09\approx 2.78$$)吻合。该现象在PBNN推断中意味着Sigmoid中段(判决边界)周围节点采样成本显著高于两端饱和区，合理调度采样资源的优先顺序应由工作点位置决定。其三，MC与CLT的有效性边界截然不同。MC方法($$M = 2\times 10^4$$)在所有9个组合上与精确解偏差均在$$\pm 5\%$$以内，但偏差方向不确定；尤其在$$(p, \varepsilon) = (0.5, 0.10)$$一例MC给出$$\hat K_{\mathrm{req}}^{\mathrm{MC}} = 101 < 106 = K_{\mathrm{req}}^{\mathrm{exact}}$$，MC因自身噪声可能给出不安全的低估，若直接用作部署阈值存在实际失效风险。CLT近似在全区间则系统性低估：$$K$$较大时($$K\geq 400$$)偏差小于4%可接受，$$K$$降至约100时偏差达10%、进一步降至约40时达16%，显示CLT在亚百采样区对Binomial离散阶梯的平滑近似不再充分。

2.3.2节对每个电压幅值采用$$K = 100$$次重复测量。由表2.13可推出隐含精度：$$p = 0.5$$时$$K = 100$$对应95% CI半宽约$$\varepsilon\approx 0.100$$；$$p = 0.1$$或$$p = 0.9$$时$$\varepsilon\approx 0.060$$。图2.13实验曲线所叠加的Wilson误差线宽度与该精度预期一致，验证$$K = 100$$作为表征用采样数的合理性。以$$E_{\mathrm{write}} = 0.78\,\mathrm{pJ}$$为基本代价单位，硬件Bernoulli采样的每概率输出总能耗$$E_{\mathrm{total}} = K\cdot E_{\mathrm{write}}$$，针对Sigmoid中段($$p = 0.5$$)，10%, 5%, 2%精度分别对应83 pJ, 311 pJ, 1.91 nJ，两侧饱和区在同一精度下能耗减少约2.6倍。这一数值关系为PBNN推断的异构采样调度提供直接定量依据：决策边界节点需更多采样、饱和节点可显著压缩、整体能耗可在同等精度下获得显著改进。

---

![sMTJ硬件Bernoulli采样可靠性分析](figs/Chapter02_local_18.png)

**图2.18** sMTJ硬件Bernoulli采样可靠性分析(Binomial精确解与MC、CLT对照，Device A, P→AP, $$t_w = 0.75\,\mathrm{ns}$$主基准工作点)。六个子面板按两行三列排列。(a)$$p = 0.5$$工作点下$$\hat p_K$$的Binomial概率质量函数在$$K\in\{5, 20, 100, 500\}$$时的离散分布，琥珀阴影带标示$$\varepsilon = 0.05$$误差带，$$K = 5$$时覆盖率为0%、$$K = 500$$达97.2%接近Gaussian极限。(b)$$\sigma(\hat p_K) = \sqrt{p(1-p)/K}$$对数-对数图，三条实线为Binomial精确值、正方形为$$M = 3000$$次MC实测，$$p = 0.1$$与$$p = 0.9$$曲线因对称性完全重合，MC markers与精确线吻合于$$K\geq 3$$全区间。(c)$$p = 0.5$$下精确覆盖率的离散阶梯曲线(三条$$\varepsilon$$水平)与$$M = 500$$的MC 95%置信带叠加，清晰展示局部非单调的阶梯倒退，五角星标示持久单调算法求解的精确$$K_{\mathrm{req}}$$。(d)MC覆盖率估计器在$$(K, p, \varepsilon) = (100, 0.5, 0.10)$$(真值$$\mathrm{cov}_{\mathrm{exact}} = 0.9431$$)的RMSE随replicate数$$M$$的衰减，200次独立seed测得的RMSE与理论$$\sqrt{\mathrm{cov}(1-\mathrm{cov})/M}$$完全吻合，1% RMSE目标对应$$M\approx 500$$。(e)三工作点$$\times$$三精度下三种估计器的$$K_{\mathrm{req}}$$对比(实心柱为精确Binomial、空心正方为MC $$M = 2\times 10^4$$、空心菱形为CLT近似)，MC与精确解符合至$$\pm 5\%$$以内而CLT全区间系统性低估、在$$K\leq 100$$区偏差达15%以上。(f)$$p = 0.5$$下精确$$K_{\mathrm{req}}(\varepsilon)$$从$$\varepsilon = 0.01$$至0.20的双对数曲线，紫色点线标示$$\propto\varepsilon^{-2}$$参考标度，水平点线标示2.3.2节的$$K = 100$$位置。

## 2.4 本章小结

本章按建模、仿真、实测三段递进，把三端SOT-sMTJ封装为可参数化、可工艺标定，并且与算法语义直接对接的Bernoulli采样接口，由此为第一章所提出的全自旋三位一体架构与时域展开范式提供唯一的物理入口——后续两章对sMTJ的任何调用都不再回到底层物理推导，只经由本章定义的Sigmoid参数$$(u_{\mathrm{th}},\beta_s)$$与五参数行为模型接入。下文综合本章工作的方法学骨架、物理与工艺层面的关键发现、当前限制以及面向后续两类任务的接口约定。

方法学骨架的核心是把T型电路、宏自旋sLLG方程、Néel-Brown热激活模型与工作区Sigmoid近似四个原本分散的层级编织成一条共享参数与驱动量的纵向叙事：T型电路在端口层面解析地解耦$$V_{\mathrm{MTJ}}$$与$$I_{\mathrm{SOT}}$$，sLLG方程在自由层内部刻画磁化矢量在SOT驱动、VCMA压低能垒与热涨落共同作用下的随机演化，Néel-Brown模型把上述微观动力学压缩为参数化翻转概率，Sigmoid近似进一步将复杂指数嵌套形式收敛到$$(u_{\mathrm{th}},\beta_s)$$两个可工艺标定的工程参数。该统一框架避免了文献中常见的同名参数在不同层级取值不一致的混乱，也为下游章节提供了单一接口集。在数值实现上，本工作将自热效应与温度依赖材料参数显式纳入sLLG求解的每步反馈链：每个时间步先以一维RC热扩散方程更新瞬时温度$$T(t)$$，再以Bloch标度律与Callen-Callen幂律修正$$M_s$$、$$K_i$$、$$\eta_{\mathrm{SH}}$$，最终回喂入有效场与SOT等效场，构成完整的$$T\to$$材料参数$$\to\mathbf{H}_{\mathrm{eff}}\to\mathbf{m}$$反馈通道。在器件实测端，本章提出D2D-C2C双层分解框架：Néel-Brown的$$V_{\mathrm{th}}$$预测精确而$$\beta_s$$过估这一表观矛盾被严格分解为正交两层——脉宽法$$\Delta_{\mathrm{pulse}}$$承载单畴热激活势垒的全部热力学含义，循环间收窄因子$$\eta_c\equiv\beta_s^{\mathrm{meas}}/\beta_s^{\mathrm{NB}}$$作为分布形态参数刻画C2C相对Gumbel基线的过度展宽。

物理发现层面，我们通过理论计算发现经典Néel–Brown热激活模型对实验测量得到的sMTJ概率分布宽度存在系统性的过度估计。基于sLLG宏自旋模型的Monte Carlo扫描显示，在$$t_w=0.75\,\mathrm{ns}$$校准点处给出的50%翻转阈值与实验$$I_{\mathrm{th}}\approx 1152\,\mu\mathrm{A}$$测量高度吻合，一致至5%以内。若将器件的自热效应也纳入考量，自热反馈开启支在阈值附近相比关闭支$$P_{\mathrm{sw}}$$高出约0.21，但未观测到明显的展宽现象，证实$$\eta_c$$的物理来源不在静态参数漂移，而必须从sLLG的其它动力学层面 (亚畴协同跃迁、热辅助进动翻转过渡区) 解释。

同时我们还发现sMTJ对$$\Delta$$扰动存在内禀的低灵敏度，这是工作点接近确定性翻转极限的几何必然。解析公式$$\partial P_{\mathrm{sw}}/\partial\Delta|_{V_{\mathrm{th}}}=-(\ln 2/2)\xi(t_w)$$将工艺容差的脉宽依赖压缩为单一无量纲量$$\xi(t_w)=\ln(t_w/(\tau_0\ln 2))/\Delta$$，$$t_w=0.75\,\mathrm{ns}$$下$$\xi\approx 0.016$$使$$P_{\mathrm{sw}}$$对$$\Delta$$几乎绝缘。而对于传统确定性MRAM存储而言，若以保持时间$$t_w \sim \tau_{\mathrm{ret}} \sim 10^{17}\,\mathrm{ns}$$ (10年) 作为评估点，则$$\xi \to 1$$、$$\partial P_{\mathrm{sw}}/\partial \Delta \to -\ln 2/2$$，此时$$\mathrm{CV}(\Delta) = 7.7\%$$即可使保持失败概率发生$$\sim 0.35 \times 0.077 \times \Delta = \mathcal{O}(1)$$量级变化，工艺容差较为苛刻，而sMTJ器件的工作时间集中在纳秒，因而对$$\Delta$$扰动具有内在低灵敏度，工艺余量天然宽裕。

工艺洞见层面，我们将$$R_P$$实测值 (10.89 kΩ) 与按物理直径$$D_{\mathrm{phys}}=80\,\mathrm{nm}$$推算的Brinkman值 (7.16 kΩ) 之间进行比较，由它们之间所存在的34%表观偏差，推导出了电学有效直径的具体数值$$D_{\mathrm{elec}}\approx 65\,\mathrm{nm}$$，这对应刻蚀损伤环带宽度$$\delta_{\mathrm{edge}}\approx 7\,\mathrm{nm}$$，可与SEM/TEM形貌表征实验进行直接比对。后续阵列级建模与可靠性分析均以$$D_{\mathrm{elec}}$$作为电学量的统一基准。

从PDK失配反推的$$\mathrm{CV}(\Delta)=7.7\%$$中，几何体积$$V_{\mathrm{mag}}$$贡献66%总方差，界面各向异性$$H_k$$贡献27%，饱和磁化$$M_s$$贡献7%，这一份额分布可用于直接指导sMTJ工艺优化的优先级：将$$\mathrm{CV}(D)$$从3.1%进一步压缩至2%可使$$V_{\mathrm{mag}}$$贡献的方差减半，对应$$\mathrm{CV}(\Delta)$$从7.7%降至约5.5%；相比之下同等比例改善$$M_s$$或$$H_k$$的边际收益仅约1/4至1/2。

宏自旋工艺失配注入进一步界定了上述结论的适用边界：NB传递函数显示$$\Delta$$单参数扰动在0.75 ns纳秒写入点只会带来很弱的斜率退化，但固定端电压下的Monte Carlo还会叠加磁性参数导致的动力学阈值漂移，以及$$R_{\mathrm{SOT}}$$到$$R_W$$的电压-电流换算漂移。结果表现为单器件Sigmoid上升段仍在，而晶圆平均$$P_{\mathrm{sw}}(V_{\mathrm{SOT}})$$在阈值区被D2D横向散布展宽；因此图2.13中更复杂的实测局部形状仍需结合C2C、方向不对称与back-hopping等机制解释。

然而当前工作依然存在诸多限制，主要集中在三处。首先，本章的实验验证仅基于Device A与Device B两个同批次单器件，尚未覆盖晶圆级多裸片多批次的统计样本，$$\eta_c$$的器件间散布需通过阵列级测试系统对数十至数百颗器件进行$$P_{\mathrm{sw}}(V,t_w)$$批量扫描后才能充分表征；其次，工艺容差分析假设$$H_k$$、$$M_s$$、$$V_{\mathrm{mag}}$$三个误差源相互独立，若实际工艺中$$H_k$$与$$M_s$$通过界面有序度共同关联，将使$$\mathrm{CV}(\Delta)$$被高估或低估10%–30%，完整刻画需要foundry提供$$R_P$$、$$\mathrm{TMR}$$、$$R_{\mathrm{SOT}}$$三者的协方差矩阵；最后，vgsot-sim当前已支持3D磁化轨迹可视化，但生产级Monte Carlo大批量扫描仍以$$m_z$$标量为主要观测量，对亚畴协同跃迁、进动相干性等无法通过$$m_z$$捕捉的微观机制的精细分析需进一步将完整$$\mathbf{m}(t)$$轨迹纳入Monte Carlo后处理通道。

总而言之，本章为后续研究提供了一个跨任务可共用的物理接口。Sigmoid参数化$$(u_{\mathrm{th}},\beta_s)$$是第三章把sMTJ抽象为伊辛自旋、第四章把sMTJ抽象为概率二值权重时唯一的数学入口——两类任务在硬件上之所以能落到同一阵列、由同一调度器在不同时间窗口分别承担，正是因为它们在器件层共用这一接口；vgsot-sim作为可调用的物理引擎，让两类任务在评估器件级C2C波动、写入能耗与Bernoulli采样精度时不必各自重建仿真链；Néel-Brown参数、同批次Sigmoid基准与工艺容差边界等同批次标定数据，则构成第三、四章在器件层接入实验校准时可直接复用的常量集。本章工作的价值因此不在于刻画了某一颗器件的具体行为，而在于把全自旋三位一体架构中所依赖的sMTJ封装成可被两类不同上层任务共同调用的统一对象。

## 注释与文献

[^note-mc-budget]: 此处推荐的$$N$$值以"在统计意义下恢复理想$$P_{\mathrm{sw}}(V)$$曲线全形态"作为评判基准，即要求$$\hat{\mathcal{F}}$$在95%置信度内逼近$$N\rightarrow\infty$$的连续概率响应。在第三章伊辛节点退火与第四章PBNN训练等下游应用中，算法仅依赖工作区内$$P_{\mathrm{sw}}$$的局部线性灵敏度而非全形态精度，所需采样次数可进一步压缩至更小量级，相关定量结果将在对应章节给出。
[^note-retention-delta]: 存储MRAM典型$$\tau_{\mathrm{ret}}>10\,\mathrm{yr}$$对应$$\Delta>60$$。
[^note-nb-slope]: 线性势垒近似下$$\beta_s^{\mathrm{NB}} = 2\Delta\ln 2/V_{c0}$$只取决于$$\Delta$$与$$V_{c0}$$，与$$t_w$$无关。
[^note-eta-fit]: 5.34为Monte Carlo数值拟合所得Device A、P→AP方向值；以表2.8解析参数代入$$\beta_s^{\mathrm{meas}}/\beta_s^{\mathrm{NB,\,analytic}} = 44.6/7.94 = 5.62$$，两者差异源于MC实现对NB拟合的轻度有限$$N$$偏差，不影响下游分析结论。
[^note-variance-sum]: 该合成值小于诸单源CV代数和12.3%，原因在于独立随机变量按方差而非标准差线性叠加：$$\sqrt{\mathrm{Var}(X+Y)} = \sqrt{\mathrm{Var}(X)+\mathrm{Var}(Y)}\leq\sqrt{\mathrm{Var}(X)}+\sqrt{\mathrm{Var}(Y)}$$。

[^note-dev-thetacalib]: 该有效值由一次自下而上的标定试错确定，并非直接取自文献。以文献β-W体系的$$\theta_{\mathrm{SH}}\approx0.25$$起步时，仿真给出的SER 50%阈值仅约$$140\,\mu\mathrm{A}$$ ($$V_{\mathrm{SOT}}\approx109\,\mathrm{mV}$$)，较同批次Device A P→AP实测的$$I_{\mathrm{th}}\approx1.09\,\mathrm{mA}$$ ($$V_{\mathrm{th}}(0.75\,\mathrm{ns})\approx844\,\mathrm{mV}$$) 低约$$4.7$$倍。依2.1.3节临界电流标度$$I_{c0}^{\mathrm{SOT}}\propto1/\theta_{\mathrm{SH}}$$，在$$\theta_{\mathrm{SH}}$$、$$K_i$$、$$M_s$$、$$\alpha$$、$$R_W$$五个候选旋钮按灵敏度分级、每点数十条轨迹的短Monte Carlo试扫中，$$\theta_{\mathrm{SH}}$$被选为主调参 ($$\alpha=0.05$$已达CoFeB典型上限不宜再增)。先调至$$0.07$$使$$t_w=5\,\mathrm{ns}$$点的$$V_{\mathrm{th}}^{\mathrm{sim}}\approx508\,\mathrm{mV}$$与实测$$511\,\mathrm{mV}$$吻合，再细调至$$0.04$$以同时命中$$0.75\,\mathrm{ns}/894\,\mathrm{mV}$$靶点；该终值随后在配套文档中统一替换了若干残留的$$0.07$$旧值。标定流程见`scripts/09_simulation_figures/calibrate_to_experiment.py`。
[^note-dev-cayley]: 此结论源于实现过程中的实测教训而非先验取舍。仿真器最初的积分核 (`dynamic_switching.switching`) 即为球坐标$$(\theta,\phi)$$下的显式Euler步，在PMA稳态$$m_z\approx\pm1$$ (即$$\theta\to0,\pi$$) 附近因运动方程的$$1/\sin\theta$$项触发数值发散；显式切线步还须逐步手动重归一化，反过来扰动了热噪声场的Stratonovich统计权重、使等效仿真温度偏离设定值。为此将内核改写为笛卡尔形式并引入下文的保模长Cayley步 (`dynamic_switching_vector`，`integrator="cayley"`)，并顺带把自旋霍尔极化方向$$\hat{\sigma}_{\mathrm{SH}}$$由原先硬编码的$$-\hat{x}$$改为显式三矢量传入，以支持任意偏置构型。
[^note-dev-hex]: 该垂直关系在标定中被确认为模型的一处方向敏感点：交换偏置场一旦偏离与$$\hat{\sigma}_{\mathrm{SH}}$$ (本仿真器默认$$-\hat{x}$$) 垂直的设置，仿真SER即塌缩为约$$0.5$$的随机比特平台，器件退化为不再受驱动电流极性确定性偏置的无偏硬币 (与上述对称性破缺机制互为印证)。故实现中将$$\mathbf{H}_{\mathrm{EX}}$$严格约束为垂直于$$\hat{\sigma}_{\mathrm{SH}}$$，默认取沿$$-\hat{y}$$的$$-50\,\mathrm{Oe}$$。

[^ref-akerman-tmr]: J. J. Akerman, J. M. Slaughter, R. W. Dave, and I. K. Schuller, "Tunneling criteria for magnetic-insulator-magnetic structures," *Applied Physics Letters*, vol. 79, pp. 3104-3106, 2001. DOI: [10.1063/1.1415412](https://doi.org/10.1063/1.1415412).
[^ref-ascher-petzold]: U. M. Ascher and L. R. Petzold, *Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations*. SIAM, 1998. DOI: [10.1137/1.9781611971392](https://doi.org/10.1137/1.9781611971392).
[^ref-berger-stt]: L. Berger, "Emission of spin waves by a magnetic multilayer traversed by a current," *Physical Review B*, vol. 54, pp. 9353-9358, 1996. DOI: [10.1103/PhysRevB.54.9353](https://doi.org/10.1103/PhysRevB.54.9353).
[^ref-bloch-law]: F. Bloch, "Zur Theorie des Ferromagnetismus," *Zeitschrift fur Physik*, vol. 61, pp. 206-219, 1930. DOI: [10.1007/BF01339661](https://doi.org/10.1007/BF01339661).
[^ref-borders-factorization]: W. A. Borders, A. Z. Pervaiz, S. Fukami, K. Y. Camsari, H. Ohno, and S. Datta, "Integer factorization using stochastic magnetic tunnel junctions," *Nature*, vol. 573, pp. 390-393, 2019. DOI: [10.1038/s41586-019-1557-9](https://doi.org/10.1038/s41586-019-1557-9).
[^ref-brinkman-bdr]: W. F. Brinkman, R. C. Dynes, and J. M. Rowell, "Tunneling conductance of asymmetrical barriers," *Journal of Applied Physics*, vol. 41, pp. 1915-1921, 1970. DOI: [10.1063/1.1659141](https://doi.org/10.1063/1.1659141).
[^ref-brown-thermal]: W. F. Brown Jr., "Thermal fluctuations of a single-domain particle," *Physical Review*, vol. 130, pp. 1677-1686, 1963. DOI: [10.1103/PhysRev.130.1677](https://doi.org/10.1103/PhysRev.130.1677).
[^ref-callen-callen]: H. B. Callen and E. Callen, "The present status of the temperature dependence of magnetocrystalline anisotropy, and the l(l+1)/2 power law," *Journal of Physics and Chemistry of Solids*, vol. 27, pp. 1271-1285, 1966. DOI: [10.1016/0022-3697(66)90012-1](https://doi.org/10.1016/0022-3697(66)90012-1).
[^ref-callen-welton]: H. B. Callen and T. A. Welton, "Irreversibility and generalized noise," *Physical Review*, vol. 83, pp. 34-40, 1951. DOI: [10.1103/PhysRev.83.34](https://doi.org/10.1103/PhysRev.83.34).
[^ref-camsari-pbits]: K. Y. Camsari, R. Faria, B. M. Sutton, and S. Datta, "Stochastic p-bits for invertible logic," *Physical Review X*, vol. 7, 031014, 2017. DOI: [10.1103/PhysRevX.7.031014](https://doi.org/10.1103/PhysRevX.7.031014).
[^ref-daquino-midpoint]: M. d'Aquino, C. Serpico, and G. Coppola, "Midpoint numerical technique for stochastic Landau-Lifshitz-Gilbert dynamics," *Journal of Applied Physics*, vol. 99, 08B905, 2006. DOI: [10.1063/1.2169472](https://doi.org/10.1063/1.2169472).
[^ref-dieny-pma-review]: B. Dieny and M. Chshiev, "Perpendicular magnetic anisotropy at transition metal/oxide interfaces and applications," *Reviews of Modern Physics*, vol. 89, 025008, 2017. DOI: [10.1103/RevModPhys.89.025008](https://doi.org/10.1103/RevModPhys.89.025008).
[^ref-garcia-palacios-sllg]: J. L. Garcia-Palacios and F. J. Lazaro, "Langevin-dynamics study of the dynamical properties of small magnetic particles," *Physical Review B*, vol. 58, pp. 14937-14958, 1998. DOI: [10.1103/PhysRevB.58.14937](https://doi.org/10.1103/PhysRevB.58.14937).
[^ref-gilbert-damping]: T. L. Gilbert, "A phenomenological theory of damping in ferromagnetic materials," *IEEE Transactions on Magnetics*, vol. 40, pp. 3443-3449, 2004. DOI: [10.1109/TMAG.2004.836740](https://doi.org/10.1109/TMAG.2004.836740).
[^ref-grimaldi-sot-mtj]: E. Grimaldi et al., "Single-shot dynamics of spin-orbit torque and spin transfer torque switching in three-terminal magnetic tunnel junctions," *Nature Nanotechnology*, vol. 15, pp. 111-117, 2020. DOI: [10.1038/s41565-019-0607-7](https://doi.org/10.1038/s41565-019-0607-7).
[^ref-ikeda-pma]: S. Ikeda et al., "A perpendicular-anisotropy CoFeB-MgO magnetic tunnel junction," *Nature Materials*, vol. 9, pp. 721-724, 2010. DOI: [10.1038/nmat2804](https://doi.org/10.1038/nmat2804).
[^ref-iserles-lie-group]: A. Iserles, H. Z. Munthe-Kaas, S. P. Norsett, and A. Zanna, "Lie-group methods," *Acta Numerica*, vol. 9, pp. 215-365, 2000. DOI: [10.1017/S0962492900002154](https://doi.org/10.1017/S0962492900002154).
[^ref-julliere-tmr]: M. Julliere, "Tunneling between ferromagnetic films," *Physics Letters A*, vol. 54, pp. 225-226, 1975. DOI: [10.1016/0375-9601(75)90174-7](https://doi.org/10.1016/0375-9601(75)90174-7).
[^ref-kim-mtj-spice]: J. Kim et al., "A technology-agnostic MTJ SPICE model with user-defined dimensions for STT-MRAM scalability studies," *IEEE Custom Integrated Circuits Conference*, pp. 1-4, 2015. DOI: [10.1109/CICC.2015.7338407](https://doi.org/10.1109/CICC.2015.7338407).
[^ref-kittel-domain]: C. Kittel, "Physical theory of ferromagnetic domains," *Reviews of Modern Physics*, vol. 21, pp. 541-583, 1949. DOI: [10.1103/RevModPhys.21.541](https://doi.org/10.1103/RevModPhys.21.541).
[^ref-kloeden-platen]: P. E. Kloeden and E. Platen, *Numerical Solution of Stochastic Differential Equations*. Springer, 1992. DOI: [10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5).
[^ref-krizakova-sot-review]: V. Krizakova, M. Perumkunnil, S. Couet, P. Gambardella, and K. Garello, "Spin-orbit torque switching of magnetic tunnel junctions for memory applications," *Journal of Magnetism and Magnetic Materials*, vol. 562, 169692, 2022. DOI: [10.1016/j.jmmm.2022.169692](https://doi.org/10.1016/j.jmmm.2022.169692).
[^ref-landau-lifshitz]: L. Landau and E. Lifshitz, "On the theory of the dispersion of magnetic permeability in ferromagnetic bodies," in *Perspectives in Theoretical Physics*. Pergamon, pp. 51-65, 1992. DOI: [10.1016/B978-0-08-036364-6.50008-9](https://doi.org/10.1016/B978-0-08-036364-6.50008-9).
[^ref-li-jiang-vcma-sot]: S. Li and Y. Jiang, "Field-free switching model of spin-orbit torque (SOT)-MTJ device with thermal effect based on voltage-controlled magnetic anisotropy (VCMA)," *AIP Advances*, vol. 13, 025030, 2023. DOI: [10.1063/9.0000426](https://doi.org/10.1063/9.0000426).
[^ref-li-zhang-thermal]: Z. Li and S. Zhang, "Thermally assisted magnetization reversal in the presence of a spin-transfer torque," *Physical Review B*, vol. 69, 134416, 2004. DOI: [10.1103/PhysRevB.69.134416](https://doi.org/10.1103/PhysRevB.69.134416).
[^ref-liu-sot-prl]: L. Liu, O. J. Lee, T. J. Gudmundsen, D. C. Ralph, and R. A. Buhrman, "Current-induced switching of perpendicularly magnetized magnetic layers using spin torque from the spin Hall effect," *Physical Review Letters*, vol. 109, 096602, 2012. DOI: [10.1103/PhysRevLett.109.096602](https://doi.org/10.1103/PhysRevLett.109.096602).
[^ref-liu-spin-hall]: L. Liu et al., "Spin-torque switching with the giant spin Hall effect of tantalum," *Science*, vol. 336, pp. 555-558, 2012. DOI: [10.1126/science.1218197](https://doi.org/10.1126/science.1218197).
[^ref-manchon-sot]: A. Manchon and S. Zhang, "Theory of nonequilibrium intrinsic spin torque in a single nanomagnet," *Physical Review B*, vol. 78, 212405, 2008. DOI: [10.1103/PhysRevB.78.212405](https://doi.org/10.1103/PhysRevB.78.212405).
[^ref-miron-sot]: I. M. Miron et al., "Perpendicular switching of a single ferromagnetic layer induced by in-plane current injection," *Nature*, vol. 476, pp. 189-193, 2011. DOI: [10.1038/nature10309](https://doi.org/10.1038/nature10309).
[^ref-moodera-tmr]: J. S. Moodera, L. R. Kinder, T. M. Wong, and R. Meservey, "Large magnetoresistance at room temperature in ferromagnetic thin film tunnel junctions," *Physical Review Letters*, vol. 74, pp. 3273-3276, 1995. DOI: [10.1103/PhysRevLett.74.3273](https://doi.org/10.1103/PhysRevLett.74.3273).
[^ref-nozaki-vcma]: T. Nozaki et al., "Understanding voltage-controlled magnetic anisotropy effect at Co/oxide interface," *Scientific Reports*, vol. 13, 10640, 2023. DOI: [10.1038/s41598-023-37422-4](https://doi.org/10.1038/s41598-023-37422-4).
[^ref-nozaki-vcma-feb]: T. Nozaki et al., "Voltage-induced magnetic anisotropy changes in an ultrathin FeB layer sandwiched between two MgO layers," *Applied Physics Express*, vol. 6, 073005, 2013. DOI: [10.7567/APEX.6.073005](https://doi.org/10.7567/APEX.6.073005).
[^ref-slonczewski-stt]: J. C. Slonczewski, "Current-driven excitation of magnetic multilayers," *Journal of Magnetism and Magnetic Materials*, vol. 159, pp. L1-L7, 1996. DOI: [10.1016/0304-8853(96)00062-5](https://doi.org/10.1016/0304-8853(96)00062-5).
[^ref-stoner-wohlfarth]: E. C. Stoner and E. P. Wohlfarth, "A mechanism of magnetic hysteresis in heterogeneous alloys," *Philosophical Transactions of the Royal Society A*, vol. 240, pp. 599-642, 1948. DOI: [10.1098/rsta.1948.0007](https://doi.org/10.1098/rsta.1948.0007).
[^ref-weinan-wang]: W. E and X.-P. Wang, "Numerical methods for the Landau-Lifshitz equation," *SIAM Journal on Numerical Analysis*, vol. 38, pp. 1647-1665, 2000. DOI: [10.1137/S0036142999352199](https://doi.org/10.1137/S0036142999352199).
[^ref-zhang-vgsot]: K. Zhang, D. Zhang, C. Wang, L. Zeng, Y. Wang, and W. Zhao, "Compact modeling and analysis of voltage-gated spin-orbit torque magnetic tunnel junction," *IEEE Access*, vol. 8, pp. 50792-50800, 2020. DOI: [10.1109/ACCESS.2020.2980073](https://doi.org/10.1109/ACCESS.2020.2980073).

