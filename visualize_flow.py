import graphviz

def create_data_flow_diagram():
    """
    创建并渲染串联模型详细数据流的有向图，基于 `串联模型研究.md` 的详细描述。
    """
    # 初始化有向图
    dot = graphviz.Digraph('串联模型详细数据流', comment='串联模型详细数据流图', engine='dot')
    dot.attr(label='串联模型详细数据流图 (基于文档2.0-2.5节)', labelloc='t', fontsize='22', fontname='Microsoft YaHei')
    dot.attr(rankdir='TD', splines='ortho') # 从上到下布局，使用正交线以减少交叉

    # 定义节点和边默认属性
    node_attrs = {
        'fontname': 'SimHei',
        'fontsize': '10'
    }
    edge_attrs = {
        'fontname': 'SimHei',
        'fontsize': '8',
        'arrowhead': 'vee',
        'arrowsize': '0.7'
    }
    dot.attr('node', **node_attrs)
    dot.attr('edge', **edge_attrs)

    # 定义不同类型的节点样式
    style_main_input = {'shape': 'folder', 'style': 'filled', 'fillcolor': '#FFDEAD', 'penwidth': '1.5', 'fontsize': '12'} # 橙黄色文件夹代表主要输入参数集合
    style_model = {'shape': 'Mrecord', 'style': 'filled', 'fillcolor': '#ADD8E6', 'penwidth': '1.5', 'fontsize': '11'} # 淡蓝色圆角矩形代表模型
    style_intermediate_data = {'shape': 'ellipse', 'style': 'filled', 'fillcolor': '#FFFFE0', 'penwidth': '1', 'fontsize': '9'} # 淡黄色椭圆代表中间数据/变量
    style_goal = {'shape': 'doubleoctagon', 'style': 'filled', 'fillcolor': '#90EE90', 'penwidth': '2', 'fontsize': '12'} # 绿色双八边形代表最终目标/输出
    style_param_cluster = {'shape': 'box', 'style': 'filled,rounded', 'fillcolor': '#E0FFFF', 'penwidth': '1', 'fontsize': '9'} # 淡青色圆角方框代表模型特定参数簇

    # --- 1. 初始输入参数 ---
    dot.node('SlurryRecipeParams', '浆料配方参数\n(固相、离聚物、溶剂的种类、含量、密度、粒径等)', **style_main_input)

    # --- 2.0 浆料初始粘度模型 ---
    with dot.subgraph(name='cluster_ViscosityModel') as sub_visc:
        sub_visc.attr(label='2.0 浆料初始粘度模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_visc.node('M_Viscosity', '{<in>输入 | 2.0 浆料初始粘度模型 | <out>输出}', **style_model)
        sub_visc.node('P_Viscosity', '模型参数\n(非吸附离聚物浓度 C_non-ads, \nHuggins常数 k_H, 特性粘度 [η], \nK-D/Quemada参数如 φ_m, [η]_particle)', **style_param_cluster)
        sub_visc.node('eta_app_0', '浆料初始表观粘度\nη_app,0', **style_intermediate_data)

    dot.edge('SlurryRecipeParams', 'M_Viscosity:in', label='浆料组分及特性')
    dot.edge('P_Viscosity', 'M_Viscosity:in', label='粘度模型特定参数')
    dot.edge('M_Viscosity:out', 'eta_app_0', label='计算得到')

    # --- 2.1 高速剪切模型 ---
    with dot.subgraph(name='cluster_ShearModel') as sub_shear:
        sub_shear.attr(label='2.1 高速剪切模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_shear.node('M_Shear', '{<in>输入 | 2.1 高速剪切模型 | <out>输出}', **style_model)
        sub_shear.node('P_Shear', '模型参数\n(n_initial, t_span, d_primary, φ_solid, η_solvent, T, \nshear_rate, ultrasound_power, F0_bond, \n凝聚/破碎系数, 孔隙率模型参数等)', **style_param_cluster)
        sub_shear.node('n_final', '最终平均团簇尺寸\nn_final', **style_intermediate_data)
        sub_shear.node('eta_slurry_final_hss', '剪切后浆料粘度\nη_slurry_final', **style_intermediate_data) # HSS模型内部会计算浆料粘度
        sub_shear.node('epsilon_n_final', '最终团簇孔隙率\nε(n_final)', **style_intermediate_data)
        sub_shear.node('h_wet', '初始湿膜厚度\nh_wet', **style_intermediate_data) # 理论上由涂布工艺和浆料粘度决定

    # 从浆料配方中获取HSS所需基础参数
    dot.edge('SlurryRecipeParams', 'P_Shear', label='溶剂粘度, 固含量, 粒径等', style='dashed')
    dot.edge('P_Shear', 'M_Shear:in', label='HSS模型特定参数')
    dot.edge('eta_app_0', 'M_Shear:in', label='初始粘度信息 (影响流动)', style='dashed') # 初始粘度影响初始流动状态和可能的部分参数选择

    dot.edge('M_Shear:out', 'n_final')
    dot.edge('M_Shear:out', 'eta_slurry_final_hss')
    dot.edge('M_Shear:out', 'epsilon_n_final')
    dot.edge('eta_slurry_final_hss', 'h_wet', label='(由涂布工艺确定)', style='dashed') # 假设关系

    # --- 2.2 干燥模型 ---
    with dot.subgraph(name='cluster_DryingModel') as sub_drying:
        sub_drying.attr(label='2.2 干燥模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_drying.node('M_Drying', '{<in>输入 | 2.2 干燥模型 | <out>输出}', **style_model)
        sub_drying.node('P_Drying', '模型参数\n(初始固含量, 溶剂饱和蒸气压, \n溶剂/溶质扩散系数, 离聚物特性, \n干燥条件: T, 湿度, 气流速率)', **style_param_cluster)
        sub_drying.node('h_dry', '干膜厚度\nh_dry', **style_intermediate_data)
        sub_drying.node('epsilon_z', '孔隙率分布\nε(z)', **style_intermediate_data)
        sub_drying.node('phi_i_z', '组分分布\nΦ_i(z)', **style_intermediate_data)

    dot.edge('SlurryRecipeParams', 'P_Drying', label='初始溶剂/固相信息', style='dashed')
    dot.edge('P_Drying', 'M_Drying:in', label='干燥模型特定参数')
    dot.edge('n_final', 'M_Drying:in', label='团簇尺寸 (影响初始孔结构)')
    dot.edge('epsilon_n_final', 'M_Drying:in', label='团簇孔隙率 (影响初始孔结构)')
    dot.edge('h_wet', 'M_Drying:in', label='湿膜厚度')

    dot.edge('M_Drying:out', 'h_dry')
    dot.edge('M_Drying:out', 'epsilon_z')
    dot.edge('M_Drying:out', 'phi_i_z')

    # --- 2.3 孔隙模型 ---
    with dot.subgraph(name='cluster_PoreModel') as sub_pore:
        sub_pore.attr(label='2.3 孔隙模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_pore.node('M_Pore', '{<in>输入 | 2.3 孔隙模型 | <out>输出}', **style_model)
        sub_pore.node('P_Pore', '模型参数\n(颗粒/团簇尺寸分布, 离聚物特性, \n界面能, 材料力学性能)', **style_param_cluster)
        sub_pore.node('pore_structure_detailed', '精细孔结构\n(PSD, τ, ε_eff)', **style_intermediate_data)
        sub_pore.node('ionomer_dist_pore', '离聚物孔壁分布', **style_intermediate_data)
        sub_pore.node('theta_eff', '有效润湿角\nθ_eff', **style_intermediate_data)

    dot.edge('SlurryRecipeParams', 'P_Pore', label='颗粒/离聚物特性', style='dashed')
    dot.edge('P_Pore', 'M_Pore:in', label='孔隙模型特定参数')
    dot.edge('h_dry', 'M_Pore:in', label='干膜厚度')
    dot.edge('epsilon_z', 'M_Pore:in', label='孔隙率分布')
    dot.edge('phi_i_z', 'M_Pore:in', label='组分分布')
    dot.edge('n_final', 'M_Pore:in', label='团簇尺寸 (影响孔结构细节)', style='dashed')


    dot.edge('M_Pore:out', 'pore_structure_detailed')
    dot.edge('M_Pore:out', 'ionomer_dist_pore')
    dot.edge('M_Pore:out', 'theta_eff')

    # --- 2.4 龟裂模型 ---
    with dot.subgraph(name='cluster_CrackModel') as sub_crack:
        sub_crack.attr(label='2.4 龟裂模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_crack.node('M_Crack', '{<in>输入 | 2.4 龟裂模型 | <out>输出}', **style_model)
        sub_crack.node('P_Crack', '模型参数\n(膜的杨氏模量 E, 泊松比 ν, \n断裂韧性 K_IC, 基底约束, 粘附能 G_adh)', **style_param_cluster)
        sub_crack.node('crack_network', '裂纹网络描述\n(密度, 长度, 宽度, 深度)', **style_intermediate_data)
        sub_crack.node('h_eff', '有效催化剂层厚度\nh_eff', **style_intermediate_data)

    dot.edge('SlurryRecipeParams', 'P_Crack', label='材料力学特性', style='dashed') # 部分力学参数可能来自原始材料
    dot.edge('P_Crack', 'M_Crack:in', label='龟裂模型特定参数')
    dot.edge('h_dry', 'M_Crack:in', label='当前膜厚')
    dot.edge('pore_structure_detailed', 'M_Crack:in', label='孔结构 (影响应力集中)')
    dot.edge('ionomer_dist_pore', 'M_Crack:in', label='离聚物分布 (影响局部力学性能)')


    dot.edge('M_Crack:out', 'crack_network')
    dot.edge('M_Crack:out', 'h_eff')

    # --- 2.5 性能模型 ---
    with dot.subgraph(name='cluster_PerformanceModel') as sub_perf:
        sub_perf.attr(label='2.5 性能模型', style='filled', color='lightgrey', fontname='Microsoft YaHei', fontsize='12')
        sub_perf.node('M_Performance', '{<in>输入 | 2.5 性能模型 | <out>输出}', **style_model)
        sub_perf.node('P_Performance', '模型参数\n(交换电流密度, 扩散系数, 电导率, \n工作条件: T, P, 浓度)', **style_param_cluster)
        sub_perf.node('Z_FinalOutput', '最终性能输出\n(极化曲线, 阻抗谱, 功率密度等)', **style_goal)

    dot.edge('P_Performance', 'M_Performance:in', label='性能模型特定参数')
    dot.edge('h_eff', 'M_Performance:in', label='有效厚度')
    dot.edge('crack_network', 'M_Performance:in', label='裂纹网络 (影响传输和活性面积)')
    dot.edge('pore_structure_detailed', 'M_Performance:in', label='孔结构 (影响传输和活性面积)')
    dot.edge('ionomer_dist_pore', 'M_Performance:in', label='离聚物分布 (影响反应位点和传输)')

    dot.edge('M_Performance:out', 'Z_FinalOutput')

    # 渲染图表
    try:
        output_filename = '串联模型详细数据流图'
        dot.render(output_filename, view=False, cleanup=True, format='pdf')
        print(f"图表已成功生成为 '{output_filename}.pdf'")
        dot.render(output_filename, view=False, cleanup=True, format='png')
        print(f"图表已成功生成为 '{output_filename}.png'")
        # print(f"Graphviz 源文件保存在: {output_filename}") # 源文件通常与pdf/png同名，无需额外说明
    except graphviz.backend.execute.CalledProcessError as e:
        print(f"生成图表失败。请确保Graphviz已安装并配置在系统PATH中。错误信息: {e}")
    except Exception as e:
        print(f"发生未知错误: {e}")

if __name__ == '__main__':
    create_data_flow_diagram()

# 确保在文件末尾有一个空行
