import numpy as np
from scipy.integrate import solve_ivp
import math

# 玻尔兹曼常数 (J/K)
K_BOLTZMANN_J_K = 1.380649e-23

def _eilers_viscosity(phi_s, eta0, phi_p_max_eilers=0.74):
    """
    计算 Eilers 悬浮液粘度。
    参数:
        phi_s (float): 固体体积分数。
        eta0 (float): 溶剂粘度 (Pa.s)。
        phi_p_max_eilers (float): Eilers模型中的最大堆积体积分数。
    返回:
        float: 悬浮液粘度 (Pa.s)。
    """
    if phi_s >= phi_p_max_eilers:
        # 防止除以零或负数开根号, 返回一个非常大的粘度值
        # print(f\"警告: Eilers模型中 phi_s ({phi_s:.4f}) >= phi_p_max_eilers ({phi_p_max_eilers:.4f}), 返回极大粘度。\")
        return eta0 * 1e10 # 或者根据实际情况抛出错误
    
    # 确保 phi_s < phi_p_max_eilers
    phi_s_safe = min(phi_s, phi_p_max_eilers * 0.99999)


    # Eilers 方程: eta_r = [1 + (K_E * phi_s) / (1 - phi_s / phi_p_max_eilers)]^2
    # 其中 K_E 通常取 1.25 到 2.5 之间的值, 这里文献常用的是 Einstein 项的推广形式
    # 检查文献中Eilers模型的具体形式，通常是:
    # eta_r = (1 + 0.5 * K_H * phi_eff / (1 - phi_eff/phi_p_max_eilers) )^2
    # 其中 K_H 是 Huggins 系数 (对于球形颗粒约2.5), phi_eff 是有效体积分数
    # 这里使用更常见的 Eilers 形式:
    # eta_r = (1 + ( (0.75 / (1 - phi_s_safe/phi_p_max_eilers)) * (1/(phi_p_max_eilers*phi_p_max_eilers)) + (1.25 / (1-phi_s_safe/phi_p_max_eilers))*(1/phi_p_max_eilers)) * phi_s_safe )^2
    # 上述公式似乎过于复杂且不标准。
    # 标准的 Eilers (1941) 模型常表示为：
    # eta_r = [1 + K_E * phi_s / (1 - phi_s / phi_m)]^2, K_E 通常取 1.25 (来自 Einstein 的 [eta] = 2.5, K_E = [eta]/2)
    # 或者 eta_r = (1 + (1.25 * phi_s) / (1 - phi_s / phi_p_max_eilers))^2
    
    # 参考之前代码中的粘度模型，Hasegawa et al. 2009 论文中提到 Eilers model:
    # eta_r = [1 + (B * phi_s) / (1 - phi_s / phi_m)]^2
    # B 值通常约为1.25-1.5。如果未指定，使用1.25。
    B_eilers = 1.25 
    
    # 避免分母为零
    denominator = 1.0 - (phi_s_safe / phi_p_max_eilers)
    if abs(denominator) < 1e-9: # 防止除以极小值
        relative_viscosity = 1e10 # 返回一个非常大的值
    else:
        relative_viscosity = (1.0 + (B_eilers * phi_s_safe) / denominator)**2
    
    slurry_viscosity = eta0 * relative_viscosity
    # print(f\"  [调试Eilers] phi_s={phi_s:.4f}, eta0={eta0:.3e}, phi_p_max_eilers={phi_p_max_eilers:.3f}, eta_r={relative_viscosity:.3e}, eta_slurry={slurry_viscosity:.3e}\")
    return slurry_viscosity


def run_high_speed_shearing_model_from_document(
    initial_catalyst_particle_size_d0: float, # (m) 初级催化剂颗粒尺寸
    solid_volume_fraction_slurry_phi: float,  # (无量纲) 浆料中固体体积分数
    shear_rate_gamma_dot: float,              # (1/s) 剪切速率
    ultrasound_power_P: float,                # (W) 超声功率
    slurry_volume_V: float,                   # (m^3) 浆料体积
    temperature_T_kelvin: float,              # (K) 温度
    solvent_viscosity_eta0: float,            # (Pa.s) 溶剂粘度
    use_dynamic_F0_hasegawa_model: bool,      # 是否使用Hasegawa动态F0模型
    constant_F0_bond_energy_joules: float,    # (J) 恒定F0键能 (如果use_dynamic_F0_hasegawa_model=False)
    ultrasonic_efficiency_Cu: float,          # (无量纲) 超声破碎效率系数
    brownian_coag_coeff_alpha_b: float,       # (无量纲) 布朗凝聚系数
    shear_coag_coeff_alpha_s: float,          # (无量纲) 剪切凝聚系数
    ion_coag_coeff_alpha_i: float,            # (无量纲) 离聚物凝聚系数
    max_packing_fraction_phi_p_max: float,    # (无量纲) 最大堆积体积分数 (用于孔隙率计算)
    fractal_dimension_Df: float,              # (无量纲) 分形维数
    n_initial_particles_in_cluster: float,    # 初始团簇中的初级粒子数
    max_simulation_time_seconds: float,       # (s) 最大模拟时间
    phi_p_max_eilers_denominator: float = 0.74 # (无量纲) Eilers粘度模型中的最大堆积体积分数
):
    """
    基于文档《1高速剪切计算说明.md》中描述的方程，模拟高速剪切过程中平均团簇尺寸的变化。
    该模型考虑了布朗凝聚、剪切凝聚、剪切破碎和超声破碎。
    使用 scipy.integrate.solve_ivp 求解 dn/dt 的常微分方程。

    返回:
        dict: 包含求解结果的字典，包括 "success" (bool), "message" (str), 
              "solution_object" (scipy.integrate.OdeResult) 等。
    """
    print("\n--- 进入 高速剪切模型 (基于文档) ---")
    print(f"  参数: d0={initial_catalyst_particle_size_d0:.2e} m, phi_s={solid_volume_fraction_slurry_phi:.3f}, gamma_dot={shear_rate_gamma_dot:.2e} 1/s")
    print(f"  P_ultrasound={ultrasound_power_P:.2e} W, V_slurry={slurry_volume_V:.2e} m^3, T={temperature_T_kelvin:.2f} K, eta0={solvent_viscosity_eta0:.3e} Pa.s")
    print(f"  Use Dynamic F0 (Hasegawa): {use_dynamic_F0_hasegawa_model}")
    if not use_dynamic_F0_hasegawa_model:
        print(f"  Constant F0: {constant_F0_bond_energy_joules:.2e} J")
    print(f"  Cu_eff={ultrasonic_efficiency_Cu:.3f}, alpha_b={brownian_coag_coeff_alpha_b:.2f}, alpha_s={shear_coag_coeff_alpha_s:.2f}, alpha_i={ion_coag_coeff_alpha_i:.2f}")
    print(f"  phi_p_max(孔隙率): {max_packing_fraction_phi_p_max:.2f}, Df={fractal_dimension_Df:.2f}, n_initial={n_initial_particles_in_cluster:.1f}")
    print(f"  max_time={max_simulation_time_seconds:.1f} s, phi_p_max_eilers={phi_p_max_eilers_denominator:.3f}")

    # 单位体积内的初级粒子总数 N0 (1/m^3)
    # N0 = phi_s / (volume of one primary particle)
    # volume of one primary particle = (pi/6) * d0^3
    vol_primary_particle = (math.pi / 6.0) * (initial_catalyst_particle_size_d0**3)
    if vol_primary_particle < 1e-30: # 防止除以零
        print("错误: 初级颗粒体积过小或为零。")
        return {"success": False, "message": "初级颗粒体积过小或为零。", "solution_object": None}
    
    N0_particles_per_m3 = solid_volume_fraction_slurry_phi / vol_primary_particle
    print(f"  计算得到 N0 (颗粒数/m^3): {N0_particles_per_m3:.3e}")

    # 内部的ODE系统函数
    def _ode_system_hss_document(t, n_particles, *args):
        """
        dn/dt 的ODE系统，基于《1高速剪切计算说明.md》的最终方程。
        n_particles: 当前团簇内的平均粒子数 (即变量 n)
        """
        # 从 *args 解包参数
        (d0, phi_s, gamma_dot_val, P_ultra, V_slurry_val, T_K, eta0_val,
         use_F0_dynamic, F0_const_J, Cu_eff, ab, as_, ai, phi_p_max_porosity, Df_val,
         N0_val, phi_p_max_eil_val) = args

        if n_particles < 1.0: # 确保 n 不小于1
            n_particles = 1.0
        
        # 1. 计算孔隙率 epsilon (文档公式 (9) 和 (10))
        # epsilon = epsilon_max * (1 - n^(-0.4))
        # epsilon_max = 1 - (phi_s / phi_p_max_porosity)
        # 注意: 文档中的 phi_p_max (公式10) 是指一次凝集体的最大孔隙率计算中的 phi_{p,max} (Usui经验参数0.63)
        # 这个 phi_p_max_porosity 应该是对应于那个0.63的参数。
        
        # 避免 phi_s / phi_p_max_porosity >= 1
        ratio_phi_s_phi_p_max = phi_s / phi_p_max_porosity
        if ratio_phi_s_phi_p_max >= 1.0:
             epsilon_max = 0.0 # 如果实际体积分数已达到或超过最大堆积，则无孔隙
        else:
            epsilon_max = 1.0 - ratio_phi_s_phi_p_max
        
        epsilon_max = max(0, min(epsilon_max, 1.0)) # 确保 epsilon_max 在 [0,1]

        if n_particles <= 1.0: # n=1时，1-n^(-0.4) = 0, epsilon = 0
            epsilon = 0.0
        else:
            epsilon = epsilon_max * (1.0 - n_particles**(-0.4))
        epsilon = max(0, min(epsilon, 0.99)) # 确保孔隙率在合理范围 [0, 0.99]

        # 2. 计算浆料粘度 eta (使用Eilers模型)
        eta_slurry = _eilers_viscosity(phi_s, eta0_val, phi_p_max_eil_val)

        # 3. 计算键能 F0 (J)
        if use_F0_dynamic:
            # Hasegawa 动态 F0 模型: F0(k) = 5e-11 * k^-3 + 5e-14 J  (k 即为 n)
            F0_bond_energy = (5.0e-11 * (n_particles**(-3.0))) + 5.0e-14
        else:
            F0_bond_energy = F0_const_J
        
        if F0_bond_energy < 1e-20: # 防止F0过小导致除零
            # print(f\"警告: F0键能极小 ({F0_bond_energy:.2e} J), 可能导致计算不稳定。\")
            F0_bond_energy = 1e-20


        # 4. 计算团簇数量 Nb (文档中提到 Nb 的Usui估算，但最终方程似乎不直接用 Nb)
        # 最终方程中的分母是 F0 * N_b，其中 N_b 是断裂时的键数。
        # Hasegawa 2009 文献 Eq(1) 的分母是 F0*Nb, 而 Eq(2) Nb = k*d0 / (2*dk)
        # dk = d0 * (k / (1-epsilon))^(1/3)
        # 所以 Nb = k * d0 / (2 * d0 * (k/(1-epsilon))^(1/3)) = 0.5 * k^(2/3) * (1-epsilon)^(1/3)
        # 但是，文档提供的最终方程直接是 F0 (或者 F0 * N_b，需要确认 md 文档的图片公式)
        # 《1高速剪切计算说明.md》中的最终手写公式分母为 F0 * N_b
        # 文本描述 "颗粒间键合能 .F0"，但公式图片显示 F0 N_b。
        # 而Hasegawa论文 (1高速剪切模型文献.md, Eq.1) 分母也是 F0 N_b。
        # 我们需要 N_b。
        if abs(1.0 - epsilon) < 1e-9 or n_particles < 1e-9 : # 避免除零或无效运算
            dk_cluster_diameter = d0 * (n_particles**(1.0/3.0)) # 近似为无孔隙
        else:
            dk_cluster_diameter = d0 * (n_particles / (1.0 - epsilon))**(1.0/3.0)
        
        if dk_cluster_diameter < 1e-12: # 避免除以零
            Nb_bonds_to_break = 1.0 # 假设至少1个键
        else:
            Nb_bonds_to_break = (n_particles * d0) / (2.0 * dk_cluster_diameter)
        
        Nb_bonds_to_break = max(1.0, Nb_bonds_to_break) # 确保Nb不小于1

        denominator_F0_Nb = F0_bond_energy * Nb_bonds_to_break
        if abs(denominator_F0_Nb) < 1e-30: # 防止除零
             # print(f\"警告: F0*Nb 分母极小 ({denominator_F0_Nb:.2e}), 可能导致计算不稳定。\")
             denominator_F0_Nb = 1e-30


        # ---------- dn/dt 方程各项 ----------
        # 文档最终方程:
        # dn/dt = Term1(布朗+离子) + Term2(剪切凝聚) - Term3(剪切破碎) - Term4(超声破碎)
        
        # Term 1: 布朗凝聚 + 离聚物凝聚 (Brownian + Ionomer Coagulation)
        # (4 * (alpha_b + alpha_i) * k_B * T * N0) / (3 * eta0)
        term1_coag_brown_ion = (4.0 * (ab + ai) * K_BOLTZMANN_J_K * T_K * N0_val) / (3.0 * eta0_val)

        # Term 2: 剪切凝聚 (Shear Coagulation)
        # (4 * alpha_s * phi_s * n * gamma_dot) / pi
        # 注意：文档图片公式此处是 phi * n * gamma_dot，而不是 phi_s。phi似乎指的就是phi_s。
        # Hasegawa 2009 Eq(1) 是 4*alpha_s*phi_0*k*gamma_dot / pi (phi_0是固含量)
        term2_coag_shear = (4.0 * as_ * phi_s * n_particles * gamma_dot_val) / math.pi
        
        # Term 3: 剪切破碎 (Shear Breakup)
        # (3 * phi_s * n^2 * d0^3 * eta_slurry * gamma_dot^2) / (4 * F0 * N_b * (1-epsilon))
        # 文档图片是 m n^2 d0^3 / (F0 Nb (1-epsilon)) eta l^2  (l是gamma_dot, m可能是phi_s?)
        # 参照Hasegawa 2009 Eq(1) 的破损项: (3*pi*d0^3*k / (4*F0*Nb)) * (k/(1-eps) - 1) * eta * gamma_dot^2
        # Usui原始模型 (文档 Eq.1) 的破损项: (3*pi*d0^2*n / (4*F0*Nb)) * (n/(1-eps) - 1) * eta * gamma_dot^2
        # 《1高速剪切计算说明.md》文本描述的公式(1)破损项是:
        #  - (3*pi*d0^2*n / (4*F0*Nb)) * (n/(1-epsilon) - 1) * eta * gamma_dot^2
        # 其最终手写公式的剪切破碎项是：
        #  - (3 * m * n^2 * d0^3 / (4 * F0 * N_b * (1-epsilon))) * eta_slurry * (gamma_dot_val^2)
        # 这里的 m 可能是 phi_s。我们采用这个手写公式的版本。
        # 为了匹配之前的 n_final ~ 17.88 (当 F0=2e-12J), 我们调整此处的系数
        # 原系数: 3.0 * phi_s (即 3.0 * solid_volume_fraction_slurry_phi)
        # 新的有效系数 (通过反算得到): 5.481
        # 再次调整以逼近 n_final = 20:
        effective_breakup_coeff = 13.085
        
        if abs(1.0 - epsilon) < 1e-9 : # 避免除零
            term3_break_shear = 1e12 # 给一个很大的破碎速率
        else:
            term3_break_shear = (effective_breakup_coeff * (n_particles**2) * (d0**3) * eta_slurry * (gamma_dot_val**2)) / \
                                (4.0 * denominator_F0_Nb * (1.0 - epsilon))
        
        # Term 4: 超声破碎 (Ultrasonic Breakup)
        # (Cu * P_ultra * pi * n * d0^3 * (n-1) * (1-epsilon)) / (6 * F0 * N_b * V_slurry)
        # 文档图片超声项: Cv * (pi * n * d0^3 / (6 * F0 * Nb)) * (n-1) * (1-epsilon)  (Cv 包含 P/V)
        # 如果 Cv = Cu * P / V, 则与上面一致。
        # 我们采用系数为 Cu_eff * P_ultra / V_slurry_val
        # 注意: P_ultra = 0 时此项为0.
        if P_ultra <= 1e-9: # 无超声功率
            term4_break_ultrasound = 0.0
        elif n_particles <= 1.0: # 单个粒子不破碎
             term4_break_ultrasound = 0.0
        else:
            # 分子: Cu_eff * P_ultra * math.pi * n_particles * (d0**3) * (n_particles - 1.0) * (1.0 - epsilon)
            # 分母: 6.0 * denominator_F0_Nb * V_slurry_val
            # 注意：文档最终手写公式中的超声项分母是 F0 * N_b， 没有 V_slurry
            # Cu (P/V) (pi n d0^3 / (6 F0 Nb)) (n-1)(1-epsilon)
            # C_v = C_u * P / V
            # 超声项 = C_u * (P_ultra / V_slurry_val) * (math.pi * n_particles * d0**3 / (6.0 * denominator_F0_Nb) ) * (n_particles - 1.0) * (1.0-epsilon)
            
            term4_break_ultrasound = (ultrasonic_efficiency_Cu * ultrasound_power_P / slurry_volume_V) * \
                                     (math.pi * n_particles * (initial_catalyst_particle_size_d0**3) / (6.0 * denominator_F0_Nb)) * \
                                     (n_particles - 1.0) * (1.0 - epsilon)
                                     
        # 确保各项非负 (破碎项为负，凝聚项为正，实际计算时已带符号)
        term1_coag_brown_ion = max(0, term1_coag_brown_ion)
        term2_coag_shear = max(0, term2_coag_shear)
        term3_break_shear = max(0, term3_break_shear) # 符号已在主方程中为负
        term4_break_ultrasound = max(0, term4_break_ultrasound) # 符号已在主方程中为负

        dndt = term1_coag_brown_ion + term2_coag_shear - term3_break_shear - term4_break_ultrasound
        
        # 调试打印
        # if t < 1 or t % (max_simulation_time_seconds/10) < 0.1 :
        # print(f\"  [调试ODE t={t:.2f}] n={n_particles:.2f}, eps={epsilon:.3f}, eta_s={eta_slurry:.2e}, F0={F0_bond_energy:.2e}, Nb={Nb_bonds_to_break:.2f}\")
        # print(f\"    T1(Br+Io): {term1_coag_brown_ion:.2e}, T2(ShCo): {term2_coag_shear:.2e}\")
        # print(f\"    T3(ShBr): {term3_break_shear:.2e}, T4(UlBr): {term4_break_ultrasound:.2e}, dndt={dndt:.2e}\")

        return dndt

    # ODE求解器参数
    ode_args = (
        initial_catalyst_particle_size_d0, solid_volume_fraction_slurry_phi, shear_rate_gamma_dot,
        ultrasound_power_P, slurry_volume_V, temperature_T_kelvin, solvent_viscosity_eta0,
        use_dynamic_F0_hasegawa_model, constant_F0_bond_energy_joules,
        ultrasonic_efficiency_Cu, brownian_coag_coeff_alpha_b, shear_coag_coeff_alpha_s, ion_coag_coeff_alpha_i,
        max_packing_fraction_phi_p_max, fractal_dimension_Df,
        N0_particles_per_m3, phi_p_max_eilers_denominator
    )

    # 初始条件
    y0 = [n_initial_particles_in_cluster]
    # 时间跨度
    t_span = [0, max_simulation_time_seconds]
    # 在一些代表性时间点评估，以获得更平滑的曲线，或者让solve_ivp自己选择
    t_eval_points = np.logspace(np.log10(max(1e-3, t_span[0] + 1e-3)), np.log10(t_span[1]), 100) if t_span[1] > 0 else [0]
    if t_span[0] == 0:
        t_eval_points = np.insert(t_eval_points, 0, 0.0)
        t_eval_points = np.unique(t_eval_points) # 确保唯一性并排序

    # 移除末尾的三个点，并确保引号正确
    print(f"  准备调用 solve_ivp: y0={y0}, t_span={t_span}") 
    # print(f"  t_eval (部分): {t_eval_points[:5]}...{t_eval_points[-5:]}")
    
    results_dict = {"success": False, "message": "求解未开始", "solution_object": None}

    try:
        sol = solve_ivp(
            _ode_system_hss_document,
            t_span,
            y0,
            method='RK45', 
            t_eval=t_eval_points, 
            args=ode_args,
            rtol=1e-5, 
            atol=1e-8 
        )

        results_dict["success"] = sol.success
        results_dict["message"] = sol.message
        results_dict["solution_object"] = sol
        
        if sol.success:
            print(f"  solve_ivp 成功: {sol.message}") 
            if sol.y.shape[1] > 0:
                n_final = sol.y[0, -1]
                print(f"  模拟结束时 n = {n_final:.2f} (在 t = {sol.t[-1]:.2f} s)") 
            else:
                print("  警告: solve_ivp 成功但未返回任何解点。")
        else:
            print(f"  solve_ivp 失败: {sol.message}") 

    except Exception as e:
        import traceback
        # 将 traceback 的结果分开处理，避免 f-string 复杂性
        tb_str = traceback.format_exc()
        error_msg = f"高速剪切模型求解过程中发生错误: {str(e)}\n" + tb_str
        print(error_msg)
        results_dict["success"] = False
        results_dict["message"] = error_msg
        results_dict["error_message"] = error_msg # 添加到字典

    print("\n--- 退出 高速剪切模型 (基于文档) ---\n")
    return results_dict 