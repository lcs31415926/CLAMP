import numpy as np
import matplotlib.pyplot as plt

# 从新的 matlab_high_speed_shearing_model_py 包导入核心功能
from matlab_high_speed_shearing_model_py import run_high_speed_shearing_model_from_document

# 尝试设置支持中文的字体
try:
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei'] # 优先使用黑体，然后是微软雅黑
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示为方块的问题
except Exception as e:
    print(f"设置中文字体失败: {e}")

def test_and_plot_hss_model():
    """
    测试高速剪切模型 (文档v2) 并绘制 n(t) 和 d_agg(t) 曲线。
    """
    print("开始测试高速剪切模型 (文档v2) - 使用 matlab_ 模块结构...")

    # --- 模型输入参数 (适配新的 run_high_speed_shearing_model_from_document) ---
    d0_test = 2.5e-6      # (m) 初级催化剂颗粒尺寸 (文献中 silica 为 2.5 um)
    phi_slurry_test = 0.15  # (无量纲) 浆料中固体体积分数 (文献 Fig.5,6 条件)
    gamma_dot_test = 1.0  # (1/s) 剪切速率 (文献 Fig.5b, 6b 条件是 1.0 s^-1)
    
    P_ultrasound_test = 0.0 # (W) 文献主要讨论无超声情况，或者说新文档模型中超声是独立项
    V_slurry_test = 100e-6  # (m^3) 100ml, 对应文档描述的实验规模
    T_kelvin_test = 273.15 + 120 # (K) 文献中 EMMA 实验温度 120 C
    eta0_solvent_viscosity_test = 190.5 # (Pa.s) 文献中 EMMA 在 120C 的粘度

    # F0 模型相关参数
    use_dynamic_F0_hasegawa_test = False # 改为False，测试恒定F0
    constant_F0_joules_test = 2.0e-12    # (J) 进一步提高恒定F0值，例如 2.0e-12 J

    # 系数参数
    Cu_efficiency_test = 0.0 # (无量纲) 超声破碎效率系数，P_ultrasound_test为0时此值无影响
                             # 如果 P_ultrasound_test > 0, 设一个较小值如0.01-0.1
    alpha_b_test = 0.6      # (无量纲) 布朗凝聚系数 (文献 Swift & Friedlander 0.6)
    alpha_s_test = 0.58     # (无量纲) 剪切凝聚系数 (文献 Swift & Friedlander 0.58)
    alpha_i_test = 0.0      # (无量纲) 离聚物凝聚系数 (测试中暂设为0)
    phi_p_max_test = 0.63   # (无量纲) 最大堆积体积分数 (Usui经验参数)
    Df_test = 1.8           # (无量纲) 分形维数 (文献中提到链状，Df可能小于2.2, Hasegawa文献取1.8)
    
    # Eilers模型中分母的最大堆积分数
    phi_p_max_eilers_test = 0.74 # (无量纲) Eilers模型参数, 文献常用0.74 for spheres

    # 模拟控制
    n_initial_test = 200.0  # 初始团簇中的初级粒子数 (根据Hasegawa文献 Fig.6 初始值估算)
    max_time_test = 1000.0  # (s) 模拟时间 (匹配Hasegawa文献 Fig.6)


    print("模型输入参数 (文档v2):")
    print(f"  d0: {d0_test:.2e} m, phi_s: {phi_slurry_test:.3f}, gamma_dot: {gamma_dot_test:.1f} 1/s")
    print(f"  P_ultrasound: {P_ultrasound_test} W, V_slurry: {V_slurry_test:.2e} m^3, T: {T_kelvin_test:.2f} K")
    print(f"  eta0: {eta0_solvent_viscosity_test:.2e} Pa.s, Use Dynamic F0: {use_dynamic_F0_hasegawa_test}")
    if not use_dynamic_F0_hasegawa_test:
        print(f"  Constant F0 (Joules): {constant_F0_joules_test:.2e} J")
    print(f"  Cu_eff: {Cu_efficiency_test:.2f}, alpha_b: {alpha_b_test:.2f}, alpha_s: {alpha_s_test:.2f}, alpha_i: {alpha_i_test:.2f}")
    print(f"  phi_p_max: {phi_p_max_test:.2f}, Df: {Df_test:.2f}, phi_p_max_eilers: {phi_p_max_eilers_test:.2f}")
    print(f"  n_initial: {n_initial_test:.1f}, max_time: {max_time_test:.1f}s")

    # 调用新的模型函数
    hss_results = run_high_speed_shearing_model_from_document(
        initial_catalyst_particle_size_d0=d0_test,
        solid_volume_fraction_slurry_phi=phi_slurry_test,
        shear_rate_gamma_dot=gamma_dot_test,
        ultrasound_power_P=P_ultrasound_test,
        slurry_volume_V=V_slurry_test,
        temperature_T_kelvin=T_kelvin_test,
        solvent_viscosity_eta0=eta0_solvent_viscosity_test,
        use_dynamic_F0_hasegawa_model=use_dynamic_F0_hasegawa_test,
        constant_F0_bond_energy_joules=constant_F0_joules_test,
        ultrasonic_efficiency_Cu=Cu_efficiency_test,
        brownian_coag_coeff_alpha_b=alpha_b_test,
        shear_coag_coeff_alpha_s=alpha_s_test,
        ion_coag_coeff_alpha_i=alpha_i_test,
        max_packing_fraction_phi_p_max=phi_p_max_test,
        fractal_dimension_Df=Df_test,
        n_initial_particles_in_cluster=n_initial_test,
        max_simulation_time_seconds=max_time_test,
        phi_p_max_eilers_denominator=phi_p_max_eilers_test
    )

    if hss_results and hss_results.get("success") and hss_results.get("solution_object"):
        print("模型成功运行。正在绘制结果...")
        sol = hss_results["solution_object"]
        
        # 提取时间和 n(t)
        time_points = sol.t
        n_t_raw = sol.y[0] # n 是解的第一个（也是唯一一个）分量
        
        # 确保 n_t 中的值不小于1.0
        n_t = np.maximum(1.0, n_t_raw)

        # 计算 d_agg(t)
        # d_agg(t) = d0 * (n(t)^(1/Df))
        d_agg_t = d0_test * (n_t**(1.0 / Df_test)) # 使用 1.0/Df 确保浮点除法

        # 绘制结果
        fig, ax1 = plt.subplots(figsize=(12, 7))

        # 绘制 n(t) vs t
        color = 'tab:red'
        ax1.set_xlabel('模拟时间 (s)', fontsize=12)
        ax1.set_ylabel('平均初级粒子数/团簇 (n)', color=color, fontsize=12)
        ax1.plot(time_points, n_t, color=color, linestyle='-', marker='.', label='n(t) (模拟值)')
        ax1.tick_params(axis='y', labelcolor=color, labelsize=10)
        ax1.tick_params(axis='x', labelsize=10)
        ax1.grid(True, linestyle='--', alpha=0.7)
        
        # 根据Hasegawa文献Fig.6b (phi=0.15, gamma_dot=1.0 s^-1)的实验数据估算
        # 时间 (s): ~0 (初始), 10, 100, 1000
        # 平均k值: ~200, ~60, ~30, ~20
        exp_time_hasegawa = np.array([0, 10, 100, 1000])
        exp_n_hasegawa = np.array([200, 60, 30, 20]) # 对应图中的Mean cluster size, k
        if gamma_dot_test == 1.0 and phi_slurry_test == 0.15 and P_ultrasound_test == 0.0:
             ax1.plot(exp_time_hasegawa, exp_n_hasegawa, color='black', linestyle='none', marker='s', markersize=8, label='n(t) (Hasegawa et al. 2009, Fig.6b 实验估算)')


        # 创建第二个y轴，共享x轴，用于绘制 d_agg(t)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('平均团簇尺寸 d_agg (μm)', color=color, fontsize=12) # 文献d0是um级别
        ax2.plot(time_points, d_agg_t * 1e6, color=color, linestyle='--', marker='x', label='d_agg(t) (μm, 模拟值)') # 转换为um
        ax2.tick_params(axis='y', labelcolor=color, labelsize=10)

        fig.tight_layout() 
        title_text = (f'高速剪切模型(文档v2): n(t) 和 d_agg(t) 随时间的变化' + 
                      f'\n(d0={d0_test*1e6:.1f}μm, $\\phi_s$={phi_slurry_test:.2f}, ' +
                      r'$\dot{{\gamma}}$=' + f'{gamma_dot_test:.1f}s$^{{-1}}$, ' +
                      f'P$_u$={P_ultrasound_test:.0f}W, n$_0$={n_initial_test:.0f})')
        plt.title(title_text, fontsize=14)
        
        # 添加图例
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=10) # ax1控制图例位置
        
        plt.show()
        print("绘图完成。请查看弹出的窗口。")

    elif hss_results:
        print("模型运行失败或未返回解对象。")
        print(f"  成功标志: {hss_results.get('success')}")
        print(f"  消息: {hss_results.get('message')}")
    else:
        print("模型运行返回了 None，发生未知错误。")

if __name__ == "__main__":
    test_and_plot_hss_model() 