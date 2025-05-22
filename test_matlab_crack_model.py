import numpy as np
import matplotlib.pyplot as plt

# 从新的 matlab_crack_model_py 包中导入各个模块
from matlab_crack_model_py import capillary_pressure_module
from matlab_crack_model_py import critical_coating_thickness_module
from matlab_crack_model_py import griffith_criteria_module
from matlab_crack_model_py import stress_function_module
from matlab_crack_model_py import tirumkudulu_russel_model_module
from matlab_crack_model_py import cracking_degree_module

def run_all_calculations():
    """执行 MATLAB main.m 中的所有计算并打印结果。"""
    print("开始执行 MATLAB main.m 中的计算...\n")

    # --- 计算毛细压力 ---
    gamma_cp = 0.072  # 液体的表面张力，单位为 N/m
    theta_cp = 30     # 接触角，单位为度
    r_p_cp = 8e-6     # 孔隙半径，单位为 m
    p_cap = capillary_pressure_module.capillary_pressure(gamma_cp, theta_cp, r_p_cp)
    print(f'毛细压力 p_cap = {p_cap:.2f} Pa')

    # --- 计算临界断裂应力 (Griffith) ---
    E_gc = 4e9        # 杨氏模量，单位为 Pa
    gamma_gc = 1.0    # 表面能，单位为 J/m^2 (MATLAB注释为 N/m，但 J/m^2更常见)
    a_gc = 0.1e-3     # 裂纹长度，单位为 m
    sigma_c = griffith_criteria_module.griffith_criteria(E_gc, gamma_gc, a_gc)
    print(f'临界断裂应力 sigma_c = {sigma_c / 1e6:.2f} MPa') # 转换为 MPa

    # --- 计算临界干燥应力 (Tirumkudulu & Russel) ---
    r_tr = 8e-6       # 颗粒半径，单位为 m
    h_layer_tr = 10e-6 # 涂层厚度，单位为 m
    G_tr = 0.6e9      # 剪切模量，单位为 Pa
    M_tr = 6          # 无量纲配位数 (MATLAB 中 M，我们 Python 版本中用 M)
    phi_rcp_tr = 0.64 # 随机密堆积的无量纲体积分数
    L_tr = 0.072      # 液体的表面张力，单位为 N/m
    sigma_crit_tr = tirumkudulu_russel_model_module.tirumkudulu_russel_model(
        r_tr, h_layer_tr, G_tr, M_tr, phi_rcp_tr, L_tr
    )
    print(f'临界干燥应力 sigma_crit = {sigma_crit_tr:.2f} Pa')

    # --- 计算临界涂层厚度 ---
    G_cct = 0.6e9      # 剪切模量，单位为 Pa (MATLAB中是0.6 GPa)
    M_cct = 6          # 无量纲配位数
    phi_rcp_cct = 0.64 # 随机密堆积的无量纲体积分数
    L_cct = 0.072      # 液体的表面张力，单位为 N/m
    # MATLAB 版本中 p_cap_max 为 -304 Pa。我们的函数期望 p_cap_max 为负值。
    # 若 MATLAB 脚本中 p_cap_max 是正值并在公式中取负，则这里直接用负值。
    p_cap_max_cct = -304 # 最大毛细压力，单位为 Pa 
    CCT = critical_coating_thickness_module.critical_coating_thickness(
        G_cct, M_cct, phi_rcp_cct, L_cct, p_cap_max_cct
    )
    print(f'临界涂层厚度 CCT = {CCT * 1e6:.2f} μm') # 转换为 μm

    # --- 应力计算 (stress_function) ---
    d_t_sf = 0.01       # 时间相关的距离或长度，单位为 m
    delta_m_t_sf = 0.001 # 时间相关的质量变化，单位为 kg
    # g_sf = 9.81       # 重力加速度 (在 Python 函数内部定义)
    L_param_sf = 0.1    # 特征长度，单位为 m (MATLAB 中参数名为 L)
    E_s_sf = 70e9       # 杨氏模量，单位为 Pa
    b_sf = 0.05         # 宽度或线性尺寸，单位为 m
    h_s_sf = 0.001      # 某一层的厚度，单位为 m
    h_c_t_sf = 0.002    # 随时间变化的另一层厚度，单位为 m
    v_s_sf = 0.3        # 泊松比
    sigma_t_sf = stress_function_module.stress_function(
        d_t_sf, delta_m_t_sf, L_param_sf, E_s_sf, b_sf, h_s_sf, h_c_t_sf, v_s_sf
    )
    print(f'随时间变化的应力 sigma(t) = {sigma_t_sf:.2f} Pa')
    print("\n计算完成。")

def plot_cracking_degree_curves():
    """调用 cracking_degree_module 中的函数并绘制曲线图。"""
    print("\n开始绘制层厚与龟裂覆盖率关系图...")
    # 设置中文显示
    plt.rcParams['font.sans-serif'] = ['SimHei'] # 指定默认字体
    plt.rcParams['axes.unicode_minus'] = False   # 解决保存图像是负号'-'显示为方块的问题

    plot_data = cracking_degree_module.get_cracking_degree_data_and_fits()
    
    plt.figure(figsize=(12, 7))
    
    # 从 cracking_degree_module.py 原始数据中获取 datasets 定义，以便获取阶数
    # (这是一个简化的示例，实际中可能需要更优雅地传递这些信息)
    # 或者直接在 get_cracking_degree_data_and_fits 返回值中包含 degree 信息
    # 更新：get_cracking_degree_data_and_fits 现在返回包含 degree 的原始数据结构
    # 这里我们直接用 cracking_degree_module 内部的 dataset keys 和参数

    # 模拟 matlab 脚本中的数据集
    # 这些是在 cracking_degree_module.get_cracking_degree_data_and_fits 中定义的
    datasets_config = {
        'CB1_IC0': {'degree': 2}, 'CB1_IC02': {'degree': 3},
        'CB1_IC035': {'degree': 3}, 'CB1_IC05': {'degree': 3},
        'CB1_IC075': {'degree': 1}, 'CB1_IC1': {'degree': 3}
    }
    dataset_keys = list(datasets_config.keys())

    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 'orange']
    # 保持与MATLAB脚本中线条样式顺序的一致性 (b-, b--, r-, m-, k-, g-)
    # 注意: MATLAB的颜色 'b--' Python中是 color='blue', linestyle='--'
    # 这里我只用了几种基本linestyle，具体对应MATLAB的绘图代码调整
    line_styles_map = {
        'CB1_IC0':   {'color': 'blue',   'linestyle': '-'},
        'CB1_IC02':  {'color': 'blue',   'linestyle': '--'}, # MATLAB 是 b--
        'CB1_IC035': {'color': 'red',    'linestyle': '-'},
        'CB1_IC05':  {'color': 'magenta','linestyle': '-'},
        'CB1_IC075': {'color': 'black',  'linestyle': '-'},
        'CB1_IC1':   {'color': 'green',  'linestyle': '-'}
    }

    legend_labels = []

    for i, key in enumerate(dataset_keys):
        # 绘制原始数据点 (可选, MATLAB脚本中注释掉了)
        # plt.plot(plot_data[f'{key}_layer_thickness'], plot_data[f'{key}_cracking_area'], 
        #          'o', color=line_styles_map[key]['color'], markersize=5, alpha=0.6)
        
        # 绘制拟合曲线
        plt.plot(plot_data[f'{key}_smooth_thickness_for_plot'], 
                 plot_data[f'{key}_fitted_cracking_for_plot'], 
                 linestyle=line_styles_map[key]['linestyle'], 
                 color=line_styles_map[key]['color'], 
                 linewidth=1.5)
        # legend_labels.append(f'{key} 拟合曲线 (poly{plot_data[f'{key}_fit_coeffs'].shape[0]-1})')
        legend_labels.append(f'{key} 拟合曲线 (poly{datasets_config[key]["degree"]})')

    plt.xlabel('层厚度 (μm)') # MATLAB 是 um, 这里用 μm 更标准
    plt.ylabel('龟裂覆盖率 (%)')
    plt.title('层厚与龟裂覆盖率关系 (Python 复现)')
    plt.legend(legend_labels, loc='upper left')
    plt.grid(True)
    plt.show()
    print("绘图完成。")

if __name__ == "__main__":
    run_all_calculations()
    plot_cracking_degree_curves() 