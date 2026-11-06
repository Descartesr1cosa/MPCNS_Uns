# MarsUns：非结构网格上的十离子火星 Hall-MHD

把结构网格程序 `MPCNS_Mercury`（分支 Mars，`Z5_Mars`）的十离子火星模型移植到非结构胞腔复形上。原程序的公式、系数与参考文档见 `MarsUns/docs/reference/`（原样拷贝，仅供参考）。
代码分三块：`Common/`（与模型无关的 MHD 部件，命名空间 `mhd`）、`MercuryUns/`（水星专有）、`MarsUns/`（火星专有，命名空间 `mars`）；火星与水星各自一个可执行文件，互不依赖。

## 文件
| 内容 | 位置 |
|---|---|
| 十离子、化学网络、大气/光电离、壳场偶极子（与原程序逐字相同，只依赖标量） | `MarsUns/physics/Mars{Species,Chemistry,AtmosphereModel,Dipoles}.*` |
| 求解器（本分支新写，含点隐式源项） | `MarsUns/driver/MarsSolver.{h,cpp}`，入口 `MarsUns/main.cpp` |
| 参数 / 网格 / 偶极子 | `MarsUns/INPUT/mars.txt`，`MarsUns/INPUT/mars_mesh.txt`，`MarsUns/INPUT/Dipoles_4840` |
| 测试 / 脚本 / 工具 | `MarsUns/tests/test_mars.cpp`（`test_mars_np1/np4`），`MarsUns/scripts/run.sh`，`MarsUns/tools/*.py` |
| 共用部件（流体通量、CT、UCT、半隐式 Hall、边界、IO） | `Common/`，见仓库 README |

运行：`MarsUns/scripts/run.sh 4 0.1`（生成无固体壳的立方球网格，内边界 = 100 km 高度 r = 1.02945 R_M，然后 `mpirun -np 4 build/MarsUns/MarsUns MarsUns/INPUT/mars.txt ...`）。

## 与结构网格版的对应
| | MarsZ5 | 本分支 |
|---|---|---|
| 流体 | 10 个五方程 Euler，局部 Rusanov，一阶 | 同；`second_order=1` 可选限制器 MUSCL |
| 电子密度、权重 | `n_e,eff=√(n_e²+10⁻⁴)`，`w_s=n_s/n_e,eff`，`u⁺` 按真实比例 | 同 |
| 理想 EMF | UCT | `induction/IdealEMF`（任意价态边） |
| 双极 EMF / 力 | `−p_scale ∇p_e/n_e,eff`，`−w_s∇p_e` | 同（节点压强体积加权，Gauss 梯度） |
| Hall | 显式边 EMF + Rusanov 型 whistler 耗散，不加 h² 步长限制 | **半隐式**（E 为未知量，GMRES，复用 `induction/ImplicitInduction`，无 taper），无条件稳定，不需要耗散项 |
| 背景场 | `Badd=IMF+壳场`，感应场 `B` 从 0 起 | 同；偶极子对每条边做两点 Gauss 积分 `A·dl`，再取 `d₁`，`div B_add≈1e-15` |
| 内边界 | 光化学平衡离子密度、T=1500 K，H⁺=0.3×来流密度，E=0 | 同，Dirichlet 虚拟态取在**边界面中心**（原程序取内侧单元中心，细网格上两者趋同） |
| 外边界 | 来流 H⁺ / 其余离子底态，`E=−u×B_imf` | 同（`farfield`：按 u·n 分入流/出流） |
| 时间推进 | 前向 Euler，`dt=min(逐离子 CFL, 总等离子体快磁声 CFL)` | 同 |

## 源项的隐式处理（按用户提示新增）
原程序只有光化学是点隐式；离子回旋、离子-中性碰撞拖曳和热松弛仍显式，网格变细/低空碰撞频率高时刚性明显（显式版本在测试网格上 t=0.1 附近 dt 塌缩到 3e-7，H⁺ 最大速度 50 ≫ 来流 15）。现在显式部分只保留 `w_s J×B` 与 `−w_s∇p_e`，其余按单元点隐式（`MarsSolver::implicit_sources_`）：
1. **回旋** `du_s/dt = (q/m_s)(u_s−u⁺)×B`：垂直 B 平面内 `z=u_x+iu_y`，`dz_s/dt=−iΩ_s(z_s−Z)`，`Z`=数密度加权平均（即 `u⁺`）。Crank–Nicolson，十种离子经 `Z` 耦合，**闭式求解** `Z'=P/(1−Q)`。A 稳定、二阶、严格守恒质量动量 `Σρ_s u_s`；动能变化计入总能（热能不变）；平行分量不变。
2. **拖曳** `dm/dt=−νm`，`ν=t_ref Σ_n n_n·2e−10`：后向 Euler，损失的动能留在总能中（变成热，和原程序一致）。
3. **热松弛 + 摩擦加热** `dp/dt=(γ−1)[k(185−T)+S_f]`：对 p 后向 Euler（T<50 K 时按 50 K 显式）。
4. 光化学：原有耦合后向 Euler + 动量/能量损失点隐式（已移植，失败时最多 16 子步）。

测试网格（4 rank，t=0.1，其余相同）：显式源项 4692 步、末期 dt=2.9e-7、`|u_H|max=50.8`；隐式源项 **1560 步**、dt≈9e-5、`|u_H|max=15.4`（≈来流），divB 同为 1e-15 量级。

## 验证（`ctest`，全部 23 项通过）
- A 远场自由流：H⁺ 变化 4e−6（氢冕碰撞/光电离的真实源），感应磁通 2e−9；
- B 单偶极子背景：`div B_add`=2e−15，与解析偶极场 L2 相对误差 1.9e−3；
- C1 低空稠密气体一步 dt=1（νdt≫1）：O₂⁺ 速度 0.5→2e−3，压强为正；C2 Ωdt≫1：质量动量变化 1.35e−3（=氢冕拖曳 νdt），速度有界；
- 整机：np=1 与 np=4 的 300 步结果在打印位数内一致；4840 偶极子（截断 1500 km）可运行；`init_fields=<fields_final.unsb>` 续算可用（背景场按当前设置重建）。

## 未做 / 注意
- 没有与结构网格 MarsZ5 的逐点对比（没有同一套网格与基准输出）；测试网格很粗（首层 0.01 R_M，未分辨 10 km 量级大气标高），只验证稳定性与结构正确性。
- 多重网格未开发（按要求）；稳态/伪时间、电阻、水星专用的源项在 Mars 分支中不使用。
- 原程序的偶极子面通量缓存（`Read_Badd_Dipoles`）未移植：这里每次启动重算（4840 源、4 万条边约 1 s 量级）。
- 化学失败（`SolvePointImplicitChemistry` 16 子步仍不收敛）会抛异常终止，与原程序一致。
