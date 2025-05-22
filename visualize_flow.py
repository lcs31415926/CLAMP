import graphviz

def create_data_flow_diagram():
    """
    创建并渲染串联模型数据流的有向图。
    """
    # 初始化有向图，设置图的标题和默认字体以支持中文
    dot = graphviz.Digraph('串联模型数据流', comment='模型参数数据流图', engine='dot')
    dot.attr(label='串联模型参数数据流图 (I/C比核心)', labelloc='t', fontsize='20', fontname='Microsoft YaHei')
    dot.attr(rankdir='TD') # 从上到下布局

    # 定义节点样式 (灵感来自Mermaid classDef)
    node_attrs = {
        'fontname': 'SimHei', # 或者其他支持中文的字体，如 'Microsoft YaHei'
        'fontsize': '10'
    }
    edge_attrs = {
        'fontname': 'SimHei',
        'fontsize': '8'
    }

    # 应用默认节点和边属性
    dot.attr('node', **node_attrs)
    dot.attr('edge', **edge_attrs)

    # 定义不同类型的节点样式
    style_main_input = {'shape': 'ellipse', 'style': 'filled', 'fillcolor': '#FFB6C1', 'penwidth': '2'}
    style_param = {'shape': 'ellipse', 'style': 'filled', 'fillcolor': '#ADD8E6'}
    style_model = {'shape': 'box', 'style': 'filled', 'fillcolor': '#E6E6FA', 'rounded': 'true'}
    style_io = {'shape': 'note', 'style': 'filled', 'fillcolor': '#FFFACD'}
    style_note = {'shape': 'plaintext', 'fontsize': '8', 'fontcolor': '#555555'} # 更像注释
    style_goal = {'shape': 'cds', 'style': 'filled', 'fillcolor': '#90EE90', 'penwidth': '1.5'}
    style_subgoal = {'shape': 'component', 'style': 'filled', 'fillcolor': '#98FB98'}


    # 1. 定义核心输入和参数节点
    dot.node('IC_Ratio', 'I/C 比', **style_main_input)
    dot.node('InitialParticleSize', '初始粒径', **style_param)
    dot.node('CatalystHydrophobicity', '催化剂亲疏水性', **style_param)

    # 2. 定义计算/中间IO节点
    dot.node('InitialViscosity', '浆料初始粘度\n(计算, 尚无具体模型)', **style_io)

    # 3. 定义模型节点
    dot.node('HSM', '高速剪切模型', **style_model)
    dot.node('PoreModel', '孔隙模型', **style_model)
    dot.node('CrackModel', '龟裂模型', **style_model)
    dot.node('DryingModel', '蒸发模型', **style_model)

    # 4. 定义模型输出/输入IO节点
    dot.node('FinalParticleSize', '最终粒径', **style_io)
    dot.node('Porosity_HSM', '总孔隙率\n(来自高速剪切)', **style_io)
    dot.node('SecondaryPoreRatio', '二次孔隙比例', **style_io)
    dot.node('SecondaryPores', '二次孔隙特性\n(如有效孔隙率, 尺寸)', **style_io)
    dot.node('CapillaryForce', '毛细管力', **style_io)
    dot.node('CatalystLayerThickness', '催化剂层厚度', **style_io)
    dot.node('IonomerDistribution', '离聚物偏析分布', **style_io)

    # 5. 定义影响因素/注释节点 (之前Mermaid中的Note节点)
    dot.node('PoreModel_Note1', '影响离聚物填充和孔结构', **style_note)
    dot.node('DryingModel_Note1', '影响溶质浓度和迁移', **style_note)
    dot.node('CapillaryForce_Note', '影响计算', **style_note)
    dot.node('OT1_Note', 'I/C低则离聚物少\n二次孔隙可能更多', **style_note)
    dot.node('OT2_Note', 'I/C低则离聚物少\n薄膜薄', **style_note)
    dot.node('PT1_Note', 'I/C大则离聚物多\n易形成网络', **style_note)
    dot.node('PS1_Note', '影响颗粒间作用和稳定性', **style_note)
    dot.node('PS1_Note2', '影响沉降', **style_note)
    dot.node('PS2_Note', '主要影响因素', **style_note)

    # 6. 定义边 (数据流)
    dot.edge('IC_Ratio', 'InitialViscosity')
    dot.edge('InitialParticleSize', 'HSM')
    dot.edge('InitialViscosity', 'HSM')

    dot.edge('HSM', 'FinalParticleSize')
    dot.edge('HSM', 'Porosity_HSM')

    dot.edge('CatalystHydrophobicity', 'PoreModel')
    dot.edge('Porosity_HSM', 'PoreModel')
    dot.edge('IC_Ratio', 'PoreModel_Note1', style='dashed', arrowhead='none')
    dot.edge('PoreModel_Note1', 'PoreModel', style='dashed', dir='forward')


    dot.edge('PoreModel', 'SecondaryPoreRatio', label='公式46等\n输入: 总孔隙率, I/C, 亲疏水性')
    dot.edge('SecondaryPoreRatio', 'SecondaryPores', label='公式42等\n输入: 总孔隙率')
    dot.edge('SecondaryPores', 'CapillaryForce', label='据公式46')
    dot.edge('CatalystHydrophobicity', 'CapillaryForce_Note', style='dashed', arrowhead='none')
    dot.edge('CapillaryForce_Note', 'CapillaryForce', style='dashed', dir='forward')


    dot.edge('CapillaryForce', 'CrackModel')
    dot.edge('CrackModel', 'CatalystLayerThickness', label='公式15')

    dot.edge('CatalystLayerThickness', 'DryingModel')
    dot.edge('IC_Ratio', 'DryingModel_Note1', style='dashed', arrowhead='none')
    dot.edge('DryingModel_Note1', 'DryingModel', style='dashed', dir='forward')
    dot.edge('DryingModel', 'IonomerDistribution')

    # 定义子图和目标
    with dot.subgraph(name='cluster_OptimizationGoal') as c_opt:
        c_opt.attr(label='核心问题: I/C比的平衡优化', style='filled', color='lightgrey', fontname='SimHei', fontsize='12')
        c_opt.node('OptGoal', 'I/C 比平衡优化', **style_goal)

    dot.edge('IC_Ratio', 'OptGoal', style='dashed', constraint='false') # Constraint false to avoid messing layout too much

    with dot.subgraph(name='cluster_OxygenGoal') as c_o2:
        c_o2.attr(label='目标1: 氧气传输', fontname='SimHei', fontsize='12')
        c_o2.node('OxygenGoal', '最大化氧气传输', **style_goal)
        c_o2.node('OT1', '高二次孔隙的孔隙率\n(I/C 小有利)', **style_subgoal)
        c_o2.node('OT2', '薄离聚物薄膜\n(I/C 小有利)', **style_subgoal)
        dot.edge('SecondaryPores', 'OT1', lhead='cluster_OxygenGoal')
        dot.edge('IonomerDistribution', 'OT2', lhead='cluster_OxygenGoal') # Or CatalystLayerThickness
        dot.edge('CatalystLayerThickness', 'OT2', style='dotted', lhead='cluster_OxygenGoal')

        dot.edge('IC_Ratio', 'OT1_Note', style='dashed', arrowhead='none')
        dot.edge('OT1_Note', 'OT1', style='dashed', dir='forward')
        dot.edge('IC_Ratio', 'OT2_Note', style='dashed', arrowhead='none')
        dot.edge('OT2_Note', 'OT2', style='dashed', dir='forward')


    with dot.subgraph(name='cluster_ProtonGoal') as c_h2:
        c_h2.attr(label='目标2: 质子传输', fontname='SimHei', fontsize='12')
        c_h2.node('ProtonGoal', '最大化质子传输', **style_goal)
        c_h2.node('PT1', '通畅离聚物网络\n(I/C 大有利)', **style_subgoal)
        dot.edge('IonomerDistribution', 'PT1', lhead='cluster_ProtonGoal')

        dot.edge('IC_Ratio', 'PT1_Note', style='dashed', arrowhead='none')
        dot.edge('PT1_Note', 'PT1', style='dashed', dir='forward')


    with dot.subgraph(name='cluster_StabilityGoal') as c_st:
        c_st.attr(label='目标3: 工艺稳定性', fontname='SimHei', fontsize='12')
        c_st.node('StabilityGoal', '优化工艺稳定性', **style_goal)
        c_st.node('PS1', '长沉淀时间\n(I/C 优)', **style_subgoal)
        c_st.node('PS2', '优良浆料流变性\n(I/C 优)', **style_subgoal)
        # InitialViscosity influences PS1 and PS2, IC_Ratio influences them
        dot.edge('InitialViscosity', 'PS1', style='dotted', lhead='cluster_StabilityGoal')
        dot.edge('InitialViscosity', 'PS2', lhead='cluster_StabilityGoal')

        dot.edge('IC_Ratio', 'PS1_Note', style='dashed', arrowhead='none')
        dot.edge('PS1_Note', 'PS1', style='dashed', dir='forward')
        dot.edge('InitialViscosity', 'PS1_Note2', style='dashed', arrowhead='none') # Added this connection
        dot.edge('PS1_Note2', 'PS1', style='dashed', dir='forward') # Added this connection

        dot.edge('IC_Ratio', 'PS2_Note', style='dashed', arrowhead='none')
        dot.edge('PS2_Note', 'PS2', style='dashed', dir='forward')


    # 连接主优化目标到具体目标
    dot.edge('OptGoal', 'OxygenGoal', style='bold')
    dot.edge('OptGoal', 'ProtonGoal', style='bold')
    dot.edge('OptGoal', 'StabilityGoal', style='bold')

    # 渲染图表
    try:
        dot.render('串联模型数据流图', view=False, cleanup=True, format='pdf')
        print("图表已成功生成为 '串联模型数据流图.pdf'")
        dot.render('串联模型数据流图', view=False, cleanup=True, format='png')
        print("图表已成功生成为 '串联模型数据流图.png'")
        print(f"Graphviz 源文件保存在: 串联模型数据流图")

    except graphviz.backend.execute.CalledProcessError as e:
        print(f"生成图表失败。请确保Graphviz已安装并配置在系统PATH中。错误信息: {e}")
    except Exception as e:
        print(f"发生未知错误: {e}")

if __name__ == '__main__':
    create_data_flow_diagram()
