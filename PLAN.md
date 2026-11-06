# MPCNS_Uns — 全非结构 Hall-MHD（水星磁层）方案与进度

最后更新：2026-11-05（M3、M5a、M4 简单版求解器完成，粗解→加密→续算全流程打通；水星算例的物理补全已写入但全模型运行在 t≈0.2–0.5 失稳，见 `docs/MercuryRun.md`，**暂存于此**）。本文件同时是**方案**和**进度记录**，每个里程碑完成后更新第 5 节。

> 多重网格 / 局部伪时间步（稳态模式）的数学、实验与开发记录见 [docs/Multigrid.md](docs/Multigrid.md)。
> 火星十离子分支（`MarsUns/`，含源项点隐式处理）见 [MarsUns/docs/Mars.md](MarsUns/docs/Mars.md)；代码按 `Common/`（模型无关）/ `MercuryUns/` / `MarsUns/` 划分，见 [README.md](README.md)。
> 数值方法的数学表达与参考文献见 [docs/Methods.md](docs/Methods.md)；水星算例的进度与对比见 [docs/MercuryRun.md](docs/MercuryRun.md)。

## 1. 目标

1. 全部网格用带符号的胞腔复形（cell → face → edge → node）描述，不再区分结构块与非结构区，不限定单元类型（首期只实现六面体）。所有场按实体编号存放。
2. 先在稀疏粗网格上算粗解。之后读取粗解，用密度、压强、`h|J|/|B|` 等指示器加最大层级，或手动指定区域，加密磁尾重联区等局部，得到**非共形**加密网格。悬挂点算作粗单元自己的顶点，粗单元变成多面体。加密后网格固定，计算中不做 AMR。
3. Preprocess 重写为通用复形生成与分区工具：读入/生成网格、构建关联表、owner 与别名符号、按进程输出子网格与通信计划。MPI 负载均衡用带权单元图（顶点权 = 单元计算成本，区分流体+电磁与只算电磁的固体壳；边权 = 共享面通信量），并报告各进程负载比与通信量。
4. 离散保持 DEC：磁场只通过 `B −= dt·d₁(E)` 更新，所有 EMF（理想、Hall、电阻、人工项）都只作为边上的 1-cochain 进入；`d₂d₁=0` 只由关联表保证，对任意单元类型、悬挂点和 MPI 界面成立。粗解转移到细网格时 B 的细化保持每个单元 `d₂B=0`。度量只出现在 Hodge 星、面积矢量、体积和重构里，对偶量保证为正。
5. Hall：沿用 `α(J×B)`，边上平均与 whistler 扩散改成按边的对偶面份额加权（任意价态），与理想 EMF 共用同一套边装配和同一个 owner；加入随 h² 收紧的 Hall 时间步限制，必要时子循环或半隐式。**Hall 首期留空并记录**（见 `docs/HallTodo.md`）。

## 2. 已确定的决定

| 事项 | 决定 |
|---|---|
| 仓库 | 单个新 git 仓库在 `MPCNS_Uns/` 根目录；旧代码原样放 `legacy/`（只读参考，不参与构建） |
| 应用名 | 求解器 `MercuryUns`，预处理 `Preprocess`（可执行 `mpcns_pre`） |
| 加密入口 | 两者都要：`mpcns_pre task=refine` 独立步骤，**以及**求解器初始化可调用的库接口 |
| 加密依据 | 流动指示器（密度、压强、`h|J|/|B|`）与手动区域（box/sphere/shell/cylinder）可单独或组合使用 |
| 输出文件 | 二进制。内部交换用自描述 UNSBIN；可视化用 VTK XML 二进制（VTU/PVTU/PVD，ParaView/VisIt）。Tecplot 见 `docs/FileFormats.md` |
| 分区 | 内置 Hilbert/Morton/RCB + 贪心边界精修；METIS 作为可选项（本机没有可用的 METIS 头文件和库） |
| 语言/依赖 | C++20，MPI（openmpi-5.0.11）；PETSc 3.26 仅在隐式电阻/Hall 半隐式阶段才接入 |

## 3. 目录结构

```
MPCNS_Uns/
├─ PLAN.md                    本文件
├─ CMakeLists.txt  cmake/
├─ include/uns/  src/uns/     共享库，按职责分目录，依赖单向
│   ├─ core/        基本类型、Params（setup 文件解析）
│   ├─ complex/     Incidence、CellComplex、ComplexBuilder、Validate（只含拓扑，不含坐标）
│   ├─ geometry/    Geometry：面积矢量、体积、对偶份额 w_{e,c}、Hodge 星、几何检查
│   ├─ generate/    HexMesh → Mesh；盒、立方球壳生成器
│   ├─ partition/   CellGraph（成本模型）、Partitioner、PartitionReport
│   ├─ parallel/    LocalMesh、CommPlan、Distributor（子复形、owner、别名符号、通信计划）
│   ├─ refine/      RefineSpec（区域/指示器）、Lattice（格点键）、RefineMesh（非共形加密）、Indicators、Transfer（B/U 转移）
│   ├─ io/          BinContainer(UNSBIN)、MeshIO、VtkWriter
│   └─ comm/        Exchanger（MPI：前向拷贝、边归约、边装配）
├─ Preprocess/                mpcns_pre：generate / inspect / export / refine / demo_fields
├─ Common/                    与模型无关的 MHD 部件（namespace mhd）：config physics state boundary fluid induction diagnostics io app
├─ MercuryUns/                水星专有：driver/ multigrid/ initial/ physics/(外逸层、回旋、双极) fluid/(Na+ 源项) induction/(内部电阻) INPUT/
├─ MarsUns/                   火星专有：driver/ physics/ INPUT/ docs/ scripts/ tools/ tests/
├─ tests/  tools/             单元与回归测试；check_vtu.py、run_selftest.sh
├─ docs/                      Conventions.md  FileFormats.md  HallTodo.md
├─ Run/                       算例（数据不入库）
└─ legacy/                    MPCNS_Mercury 与旧 Preprocess 的原样拷贝
```

求解器中旧的结构块模块（`1_grid`、`2_topology`、`3_field`、`4_halo`、`7_metric`、`8_dec`）不再沿用：它们以 `(block,i,j,k)` 和 ghost 层为核心，无法表示实体编号存储。可继承的是物理公式（流体 RHS、Hall 限幅与 taper、偶极子背景、边界条件语义）、PETSc 隐式电阻的结构和算例参数风格，从 `legacy/` 读取移植。

## 4. 关键设计（摘要）

- **复形**：`e2n`（−1 起点/+1 终点）、`f2e`（有序闭环）、`c2f`（外向符号）的带符号 CSR，另有转置表。规范定向由全局节点号推出，各进程无需通信即一致；别名符号仍存表，使导入网格也适用。详见 `docs/Conventions.md`。
- **Owner**：实体的 owner = 关联单元所属进程中编号最小者。子复形 = owned 单元 + 一层共点邻居，因此 owned 单元的每条边、每个面都在本地看到全部关联单元。
- **边装配**：每个进程对自己 owned 单元的贡献求部分和，别名边归约到 owner，再广播回所有副本（`Exchanger::assemble_edges`）。理想、Hall、电阻、人工项共用这一条路径，所以同一条边的所有副本拿到相同的 E。
- **CT**：`faraday_update` 对所有 owned 单元的面做 `B_f −= dt Σ f2e·E`。每个面的两个副本用相同的 E 和相同的规范环顺序，得到相同的更新。浮点下 `d₂B` 保持在圆整误差量级（实测 ~1e-15），不是位级为零；`d₂d₁=0` 本身是关联表上的整数恒等式。
- **加密**（M3，已实现）：每个粗六面体带整数格点 `[0,2^L]³`，粗面/粗边上的点用规范坐标系下的整数键去重，不比较浮点坐标。叶面片若对侧存在更小的面片则让位给它们，粗单元把这些小面片都记作自己的面；每条边在其上的所有节点处剖开。粗 B 转细 B：先把粗面通量按重构场在面片中心的法向值分配并重整化到总和等于粗通量，再对每个被加密的粗单元解一个小的 Hodge 投影使内部各细单元 `d₂B=0`。
- **Hall**：见 `docs/HallTodo.md`。

## 5. 里程碑与进度

状态：✅ 完成并通过测试　🟡 部分完成　⬜ 未开始

| # | 内容 | 状态 | 说明 |
|---|---|---|---|
| M0 | 仓库、legacy 拷贝、CMake 骨架、基线提交 | ✅ | 基线提交 = 旧代码原样拷贝 |
| M1 | 复形、构建器、验证器、几何/Hodge 星、盒与立方球壳生成器 | ✅ | `test_complex_geometry` |
| M2 | Preprocess：分区、子复形、owner/别名/通信计划、UNSBIN、VTU、报告 | ✅ | `test_partition` `test_distribute` `test_io` + `check_vtu.py` |
| M2b | MPI 通信层（前向、边归约、边装配） | ✅ | `test_comm_mpi`（np=1/3/6） |
| M4a | MercuryUns 骨架：FieldSet、边装配、`B −= dt d₁E`、∇·B 诊断；Hall 空钩子 | ✅ | `tools/run_selftest.sh`：np=1/2/5/8 的 `d₂B` ≤ 1.4e-15，B 校验和一致到 12 位 |
| M2c | 读入器：旧 plot3d+fvbnd、Gmsh 六面体、原生复形 | ⬜ | 旧 plot3d 读取可从 `legacy/Preprocess` 移植 |
| M2d | METIS 接入（可选） | ⬜ | 需要你提供 METIS 路径；接口已留（`PartMethod::Metis`） |
| M3 | 加密：指示器、2:1 平衡、非共形复形构造、B/U 转移、`task=refine` | ✅ | `test_refine_mesh`、`test_refine_spec`、`tools/run_refine_demo.sh`；求解器初始化调用库接口随 M4 接入 |
| M5a | 理想 EMF 推广到任意价态边，并与旧公式对照 | ✅ | `docs/IdealEMF.md`、`test_ideal_emf`（np=1/4/7）；价态 3–8、悬挂点多面体、分区无关 |
| M4 | 求解器：流体 FV、二阶重构、SSP-RK2、CT、边界、输出/重启 | ✅ | 简单版：`docs/Solver.md`、`test_solver`、`tools/run_pipeline_demo.sh`；未含电离/阻力/双极/电阻/Hall |
| M5 | 与旧 Mercury 的物理对照、电阻固体壳、电离/阻力/双极项 | 🟡 | 代码已写：电阻固体壳（隐式并行 CG）、电离/阻力/热松弛、双极项、偏心椭球网格、参考参数文件；**全模型失稳未解决**，t=25 未跑，流固负载标定任务未写（`docs/MercuryRun.md`） |
| M6 | Hall、人工项、电阻、Hall 时间步 | ⬜ | Hall 按用户要求首期留空 |
| M7 | Hall 子循环（A），必要时半隐式（B） | ⬜ | |
| M9 | 稳态加速：局部伪时间步 + 代数聚合 FAS 多重网格（流体部分，B 在粗层冻结） | 🟡 | 代码、单元测试、粗/细网格对比已做；数学与全部实验记录见 `docs/Multigrid.md`；脚本 `scripts/05_steady_run.sh`；快速模式（R=32+4 层 MG）约少 21 倍迭代但稳态误差 15–70 %，R=4 约少 8 倍迭代、误差几 %；B 的收敛是剩余瓶颈，下一步见 `docs/Multigrid.md` 第 8 节 |
| M8 | 粗解 → 加密 → 续算全流程；Tecplot 输出；文档 | 🟡 | 全流程已在合成算例上打通（`tools/run_pipeline_demo.sh`）；Tecplot 输出、性能与真实算例未做 |

### M3 的实现要点与限制
- 每个粗六面体内是八叉树；叶子按「手动区域（叶子的角点或中心落入区域）」与「粗单元级的指示器最小层级」决定是否细分。粗单元必须是拓扑六面体（否则报错）。
- 悬挂点节点用整数格点键去重，不比较浮点坐标；面 = 两侧面片的公共加密；每条边在其上所有节点处剖开（用中点是否存在逐级二分查找）。粗单元自然成为多面体，实测最多 15–21 个面。
- 面平衡（相邻单元层级差 ≤ 1）迭代到收敛；边/点邻居的层级差不强制限制（拓扑上任意层级差都有效）。
- 新节点坐标：三线性插值；`refine_spherical_snap=1` 时保持径向插值（球壳几何）。
- B：粗面通量按粗单元重构向量 `b·A` 分配到面片，再把差额按面积补齐，使面片通量之和严格等于粗通量；粗单元内部各细单元用最小范数投影（CG）加生成树扫描使 `d₂B=0`。
- U：粗单元的限幅线性（Barth–Jespersen）重构，并扣除体积加权平均增量，使每个粗单元内的细单元积分保持不变（以细网格体积计）。
- 预算：`refine_budget_cells>0` 时按 ×1.5 放宽指示器阈值；手动区域不放宽，若区域本身超预算会给出警告。
- 限制：只支持六面体粗单元；只做各向同性 2×2×2 剖分；粗解的 B 面通量重构用单元级常向量（面片通量分配为一阶精度，但散度严格保持）。

### 已验证（以及怎么验证的）
- 拓扑：`d₁d₀=0`、`d₂d₁=0`（整数精确）、环闭合、面-单元双覆盖符号相反、欧拉数（盒 1，球壳 2）。随机翻转局部边定向后恒等式仍成立。
- 几何：单元闭合 `|ΣA|/Σ|A|` ~1e-17，体积与全部 Hodge 星为正，对偶面片与边方向同向，`Σ_c w_{e,c}=1`。立方球壳含价态 3 的角点。
- 分区：Hilbert/Morton/RCB，np=1…16，负载比 ≤ 1.03；精修不增加割量。
- 子复形：各进程拓扑自洽，owned 数之和等于全局数；发送/接收列表两侧逐项一致；边装配模拟结果等于全局关联单元数。
- MPI：真实 `mpirun` 下的边装配与面前向拷贝逐项正确。
- VTU：独立的 Python 解析器重算多面体体积，与写出的体积数组一致到 1e-15（加密网格含 15 面多面体，4e-16）。
- 加密拓扑：盒（均匀、单个深层单元、随机层级、区域穿过粗单元）和带固体壳的球面壳（含价态 3 的角点）上，`d₁d₀=0`、`d₂d₁=0`、环闭合、面双覆盖符号相反、欧拉数、体积为正、单元闭合均通过；平面几何下加密前后总体积相等到 1e-12；边界面标记保留。
- B 转移：随机 `B=d₁A` 的粗解转移后，每个细单元的 `|d₂B|/Σ|B_f|` ≤ 1.3e-14（演示算例 138k 细单元为 4.2e-14）；面片通量和等于粗通量到 1e-15；均匀场在平面网格上被精确复现（1e-15）。
- U 转移：每个粗单元内细单元的体积平均等于粗值到 1e-14。
- 非共形、8 进程分区的细网格上，MercuryUns 的 CT 自检 `d₂B` 仍为 1.3e-15。
- 求解器（`docs/Solver.md`）：均匀流在悬挂点多面体网格上保持到 1e-16；密度波二阶（阶数 2.11、1.99）；Sod 激波管 L1 误差 3.5e-3；Alfvén 波包速度、幅值、符号正确（0.0249 对 0.0250，位置 2.975 对 3.0）；水星壳网格 np=4 与串行 U 差 4e-16、Phi 差 3e-14；整条流程（粗网格 1728 → 细网格 15322，t=6 续算到 t=8，2452 步，8 进程）∇·B ≤ 3e-14。
- 理想 EMF（`docs/IdealEMF.md`）：规则六面体上随机 u、B 与旧公式的转录逐边相差 ≤ 4e-15（450 条内部边）；均匀 B 的迎风项在价态 3、5、6、8 的边和悬挂点多面体上为 1e-14 以下；均匀 u×B 在平面网格上精确；20 步 CT 后 np=1/4/7 与串行的 B 相差 ≤ 9e-16，∇·B ≤ 2e-15。测试中发现并修掉了 3×3 正则化求解沿边方向灾难性相消导致结果依赖求和顺序的问题（改为边法平面内的 2×2 求解）。

### 尚未验证 / 已知缺口
- VTU 在 ParaView 里的实际打开**没有测过**（本机无 ParaView/VTK）；Tecplot 同样没有测过。
- Hall 还没写（空钩子）；电离、阻力、双极电场、电阻固体壳也没写。`task=selftest` 仍用合成 EMF，真实 EMF 已在 `task=run` 中使用。
- 求解器只验证了数值性质（二阶、激波、Alfvén 波、分区无关、∇·B），**没有**与旧 Mercury 的物理结果对照。
- 理想 EMF 的限制见 `docs/IdealEMF.md`：弯曲面一阶一致、域边界边上的迎风项不严格为零、极不规则边上权重 ω 可能变小/变负（无限幅器）。
- 指示器（密度/压强梯度、`h|J|/|B|`）只用合成的粗解测过管线，没有用真实的粗解校准阈值。
- 求解器初始化时直接调用加密库（而不是独立的 `task=refine` 步骤）还没接入；目前是 `mpcns_pre task=refine` + 求解器 `init_fields=` 两步，库接口 `build_refined_mesh / transfer_B / prolong_cells` 可直接被求解器调用。
- 未做性能测量。

## 6. 构建与运行

```sh
export PATH=/workspace/software/openmpi-5.0.11/bin:$PATH
cmake -S . -B build && cmake --build build -j
cd build && ctest                       # 14 个测试
./build/Preprocess/mpcns_pre Preprocess/INPUT/setup.txt output_dir=Run/demo/geometry proc_num=8
tools/run_selftest.sh build /tmp/uns_selftest      # 生成 + 多进程 CT 自检
tools/run_refine_demo.sh build /tmp/uns_demo       # 粗网格 → 合成粗解 → 加密 → 分区 → 细网格 CT 自检
tools/run_pipeline_demo.sh build /tmp/uns_pipe     # 粗网格 → 粗解求解 → 加密 → 细网格续算（约 30 s，最多 8 进程）
```
ParaView 打开 `Run/demo/geometry/partition.pvtu`，按 `rank` 着色即可检查分区。

## 7. 风险与待定

- 多面体、非平面面的体积与重构质量需在悬挂网格上重点验证（M3 之后）。
- 非共形面上的二阶重构与限幅较敏感。
- Hall 的加权 whistler 扩散需要推导加测试；Hall 步长在加密 ≥ 2 层时预计成为瓶颈。
- 首期加密只做各向同性 2×2×2。
- 加密网格很快变大（演示中 3840 → 2.3 万；区域取 level 3 时 → 13.8 万），需要用预算和区域大小控制。
