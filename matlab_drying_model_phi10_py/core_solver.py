import numpy as np
from scipy.integrate import solve_ivp

# 最大堆积密度，来自原始脚本的隐式定义
PHI_MAX_PACKING = 0.64
DEFAULT_NOMINAL_TEMP_K = 298.15 # 假设一个名义参考温度，例如 25°C

# 原始 K 定义 (基于平均体积分数)
def K_original(phi1, phi2):
    """
    根据 phi_10.py 计算沉积系数 K。
    使用整个数组 phi1 和 phi2 的平均值。
    """
    phi1_mean = np.mean(phi1)
    phi2_mean = np.mean(phi2)
    
    current_sum = 0.0 # 默认值
    # 仅当均值都为有限数时才计算总和
    if np.isfinite(phi1_mean) and np.isfinite(phi2_mean):
        current_sum = phi1_mean + phi2_mean
    
    value = 1.0 - current_sum
    
    # 如果 value <= 0 (即总平均体积分数 >= 1), K 应该非常小 (接近0)
    # 以避免负数开非整数次方导致nan或复数
    if value <= 1e-9: # 使用一个小的阈值
        return (1e-9)**6.55 # 返回一个非常小的正数
    return value**6.55

# 原始 Z 定义 (基于平均体积分数)
def Z_original(phi1, phi2):
    """
    根据 phi_10.py 计算压缩系数 Z。
    使用整个数组 phi1 和 phi2 的平均值。
    """
    phi1_mean = np.mean(phi1)
    phi2_mean = np.mean(phi2)
    value = PHI_MAX_PACKING - phi1_mean - phi2_mean
    # 避免除以零或负数，如果 value <= 0，则 Z 趋于无穷或无定义，这里设为一个大数或根据实际情况处理
    # phi_10.py 中是直接除，可能会产生 inf
    if value <= 1e-9: # 防止除以零
        return PHI_MAX_PACKING / 1e-9
    return PHI_MAX_PACKING / value

# 原始 PDE 系统定义 (来自 phi_10.py)
# 注意: Pe1 和 Pe2 参数现在是有效Peclet数 Pe1_eff, Pe2_eff
def pde_system_phi10_formulation(t, y, Pe1_eff, Pe2_eff, xi, N_points):
    """
    PDE 系统，严格按照 phi_10.py 的公式。
    Pe1_eff, Pe2_eff 是考虑了温度和蒸发影响的有效Peclet数。
    """
    phi1 = y[:N_points]
    phi2 = y[N_points:]

    # 梯度计算，对应 phi_10.py 中的 np.gradient(phi, xi)
    # 使用 edge_order=2 以获得更准确的边界梯度，尽管 phi_10.py 未指定，但这是 np.gradient 的默认行为之一
    dphi1_dxi = np.gradient(phi1, xi, edge_order=2)
    dphi2_dxi = np.gradient(phi2, xi, edge_order=2)

    # 为避免除以零，对梯度比值进行处理
    safe_dphi1_dxi = np.where(dphi1_dxi == 0, 1e-9, dphi1_dxi) # 避免 dphi1_dxi 为0
    safe_dphi2_dxi = np.where(dphi2_dxi == 0, 1e-9, dphi2_dxi) # 避免 dphi2_dxi 为0
    
    # 防止 Pe1_eff 或 Pe2_eff 为零导致除零错误
    safe_Pe1_eff = Pe1_eff if Pe1_eff != 0 else 1e-9
    safe_Pe2_eff = Pe2_eff if Pe2_eff != 0 else 1e-9

    # term1 和 term2 的计算，严格遵循 phi_10.py
    # 对于 phi1
    # 确保 (phi2/phi1) * (dphi1_dxi / dphi2_dxi) 中的 phi1 和 dphi2_dxi 不为零
    ratio_grad_phi1 = (phi2 / np.where(phi1 == 0, 1e-9, phi1)) * (dphi1_dxi / safe_dphi2_dxi)
    # 确保 (phi1/phi2) * (dphi2_dxi / dphi1_dxi) 中的 phi2 和 dphi1_dxi 不为零
    ratio_grad_phi2_for_phi1_denominator = (phi1 / np.where(phi2 == 0, 1e-9, phi2)) * (dphi2_dxi / safe_dphi1_dxi)


    term1_phi1_numerator = (K_original(phi1, phi2) * phi1 * 
                            (phi1 * (1 - phi1) * ratio_grad_phi1 - 
                             phi1 * phi2 * ((safe_Pe1_eff / safe_Pe2_eff)**3 + ratio_grad_phi1) + 
                             phi2 * (1 - phi2) * (safe_Pe1_eff / safe_Pe2_eff)**3))
    
    term1_phi1_denominator = ((1 - phi1) * (phi1**2) * ratio_grad_phi1 + 
                              2 * phi1 * phi2 * (safe_Pe1_eff / safe_Pe2_eff)**3 + 
                              (phi2**2) * (safe_Pe1_eff / safe_Pe2_eff)**6 * ratio_grad_phi2_for_phi1_denominator)
    
    term1_phi1 = term1_phi1_numerator / np.where(term1_phi1_denominator == 0, 1e-9, term1_phi1_denominator)
    term2_phi1 = (phi1 + (safe_Pe1_eff / safe_Pe2_eff)**3 * phi2) * Z_original(phi1, phi2)

    # 对于 phi2
    # 确保 (phi1/phi2) * (dphi2_dxi / dphi1_dxi) 中的 phi2 和 dphi1_dxi 不为零
    ratio_grad_phi2 = (phi1 / np.where(phi2 == 0, 1e-9, phi2)) * (dphi2_dxi / safe_dphi1_dxi)
    # 确保 (phi2/phi1) * (dphi1_dxi / dphi2_dxi) 中的 phi1 和 dphi2_dxi 不为零
    ratio_grad_phi1_for_phi2_denominator = (phi2 / np.where(phi1 == 0, 1e-9, phi1)) * (dphi1_dxi / safe_dphi2_dxi)

    term1_phi2_numerator = (K_original(phi1, phi2) * phi2 * 
                            (phi2 * (1 - phi2) * ratio_grad_phi2 - 
                             phi2 * phi1 * ((safe_Pe2_eff / safe_Pe1_eff)**3 + ratio_grad_phi2) + 
                             phi1 * (1 - phi1) * (safe_Pe2_eff / safe_Pe1_eff)**3))
    
    term1_phi2_denominator = ((1 - phi2) * (phi2**2) * ratio_grad_phi2 + 
                              2 * phi2 * phi1 * (safe_Pe2_eff / safe_Pe1_eff)**3 + 
                              (phi1**2) * (safe_Pe2_eff / safe_Pe1_eff)**6 * ratio_grad_phi1_for_phi2_denominator)
    
    term1_phi2 = term1_phi2_numerator / np.where(term1_phi2_denominator == 0, 1e-9, term1_phi2_denominator)
    term2_phi2 = (phi2 + (safe_Pe2_eff / safe_Pe1_eff)**3 * phi1) * Z_original(phi1, phi2)

    dterm1_phi1_dxi = np.gradient(term1_phi1, xi, edge_order=2)
    dterm2_phi1_dxi = np.gradient(term2_phi1, xi, edge_order=2)
    
    dterm1_phi2_dxi = np.gradient(term1_phi2, xi, edge_order=2)
    dterm2_phi2_dxi = np.gradient(term2_phi2, xi, edge_order=2)

    # PDE 方程，与 phi_10.py 保持一致
    # 避免 (1-t)**2 为零
    one_minus_t_sq = (1 - t)**2
    if one_minus_t_sq == 0: one_minus_t_sq = 1e-9

    dphi1_dt = (1 / (safe_Pe1_eff * one_minus_t_sq)) * dterm1_phi1_dxi * dterm2_phi1_dxi - (xi / (1 - t if t != 1 else 1e-9)) * dphi1_dxi
    dphi2_dt = (1 / (safe_Pe2_eff * one_minus_t_sq)) * dterm1_phi2_dxi * dterm2_phi2_dxi - (xi / (1 - t if t != 1 else 1e-9)) * dphi2_dxi
    
    # 边界条件，与 phi_10.py 保持一致
    dphi1_dt[0] = 0
    dphi2_dt[0] = 0
    # 尝试在 xi=1 处也设置 dphi/dt = 0，模拟无通量，尽管 phi_10.py 未显式如此做，但梯度计算可能受益
    # dphi1_dt[-1] = 0 
    # dphi2_dt[-1] = 0

    return np.concatenate((dphi1_dt, dphi2_dt))

# 事件函数也需要接收有效Peclet数，尽管在这个特定事件中未使用它们，但为了保持solve_ivp中args的一致性
def event_reached_phi_max_original(t, y, Pe1_eff, Pe2_eff, xi, N_points):
    """
    当任何一个组分的最大体积分数之和达到 phi_max 时，终止计算。
    Pe1_eff, Pe2_eff 是为了与pde_system的参数签名一致。
    """
    phi1 = y[:N_points]
    phi2 = y[N_points:]
    # phi_10.py 使用 np.max(phi1) + np.max(phi2)
    # 这可能导致如果一个组分在某点很高，另一个组分在另一点很高，它们的和也可能触发事件
    # 更合理的可能是检查同一点的 phi1+phi2 是否超限，或者整体平均值是否超限
    # 但为保持一致，这里使用 np.max(phi1) + np.max(phi2)
    if np.max(phi1) + np.max(phi2) >= PHI_MAX_PACKING:
        return 0  # 事件触发
    return 1      # 事件未触发
event_reached_phi_max_original.terminal = True
event_reached_phi_max_original.direction = -1 # 从正到负时触发

def calculate_stratification(phi2_final_profile, phi2_initial_profile, N_points):
    """计算组分2在顶部10%区域的富集百分比。"""
    if phi2_final_profile is None or phi2_initial_profile is None or N_points == 0:
        return np.nan
    
    n_top_10_percent = max(1, int(np.ceil(0.1 * N_points)))
    
    # 确保最终剖面是1D的
    if phi2_final_profile.ndim > 1:
        if phi2_final_profile.shape[1] > 0:
             phi2_final_slice = phi2_final_profile[:, -1] # 取最后一个时间点
        else:
            return np.nan # 没有最终时间点数据
    else:
        phi2_final_slice = phi2_final_profile

    if len(phi2_final_slice) != N_points:
        # 可能由于提前终止，最终剖面点数与N_points不一致，此时可能无法直接比较
        # 或者如果sol.y直接被使用且未正确切片
        print(f"警告: calculate_stratification 中 phi2_final_slice长度 ({len(phi2_final_slice)}) 与 N_points ({N_points}) 不匹配")
        # 尝试截取或填充，或者直接返回nan。为简单起见，如果长度不匹配则返回nan
        return np.nan
        
    phi2_top_avg_final = np.mean(phi2_final_slice[-n_top_10_percent:])
    phi2_avg_initial = np.mean(phi2_initial_profile)

    if phi2_avg_initial < 1e-9: # 避免除以零
        if phi2_top_avg_final > 1e-9:
            return np.inf # 初始为0，最终不为0，视为无限大富集
        else:
            return 0.0 # 初始为0，最终也为0
            
    stratification_percent = ((phi2_top_avg_final - phi2_avg_initial) / phi2_avg_initial) * 100.0
    return stratification_percent

def solve_pde_phi10_based(
    phi1_0_array, 
    phi2_0_array, 
    Pe1_nominal,  # 名义Pe1
    Pe2_nominal,  # 名义Pe2
    E_parameter, 
    N_points=10,
    tau_min=0.17, 
    tau_max=0.99, 
    num_tau_points=None, 
    t_eval_points=None,
    actual_drying_temp_K=DEFAULT_NOMINAL_TEMP_K, # 实际干燥温度 (K)
    nominal_ref_temp_K=DEFAULT_NOMINAL_TEMP_K, # Pe数定义时的参考温度 (K)
    E_evaporation_rate_factor=1.0       # 蒸发速率调整因子
    ):
    """
    求解基于 `phi_10.py` 核心逻辑的PDE系统，并加入温度和蒸发速率影响。

    返回:
    - tuple: (sol_t, phi1_solutions, phi2_solutions, phi_void_solutions, xi_grid, 
              stratification_phi2_percent, Pe1_eff, Pe2_eff, message, status, E_param_out)
    """
    if len(phi1_0_array) != N_points or len(phi2_0_array) != N_points:
        raise ValueError(f"初始条件数组 phi1_0_array 和 phi2_0_array 的长度必须等于 N_points ({N_points})")

    # 计算有效Peclet数
    if actual_drying_temp_K <= 0 or nominal_ref_temp_K <=0:
        print("警告: 温度参数为零或负，Pe数调整可能不正确。将使用名义Pe数。")
        temp_ratio = 1.0
    else:
        temp_ratio = nominal_ref_temp_K / actual_drying_temp_K
    
    Pe1_eff = Pe1_nominal * E_evaporation_rate_factor * temp_ratio
    Pe2_eff = Pe2_nominal * E_evaporation_rate_factor * temp_ratio

    xi_grid = np.linspace(0, 1, N_points)
    
    # 初始条件
    y0 = np.concatenate((phi1_0_array, phi2_0_array))

    # 时间点
    if t_eval_points is not None:
        t_eval = np.array(t_eval_points)
        # 确保 tau_max 是 t_eval 中的最大值，或者如果 t_eval_points 为空则使用 tau_min
        actual_tau_max = max(np.max(t_eval) if len(t_eval) > 0 else tau_min, tau_max)
        # 确保 tau_min 是 t_eval 中的最小值，或者如果 t_eval_points 为空则使用 tau_min
        actual_tau_min = min(np.min(t_eval) if len(t_eval) > 0 else tau_max, tau_min)

    elif num_tau_points is not None:
        t_eval = np.linspace(tau_min, tau_max, num_tau_points)
        actual_tau_min = tau_min
        actual_tau_max = tau_max
    else: # 默认行为，如果两者都未提供，可能需要一个默认的 num_tau_points
        t_eval = np.linspace(tau_min, tau_max, 10) # 默认10个点
        actual_tau_min = tau_min
        actual_tau_max = tau_max
        print("警告: num_tau_points 和 t_eval_points 均未提供，默认使用10个时间点。")

    # 确保求解区间至少是 [actual_tau_min, actual_tau_max]
    # 如果 t_eval 存在，solve_ivp 的 t_span 应覆盖 t_eval 的范围
    solve_span_min = actual_tau_min
    solve_span_max = actual_tau_max
    if t_eval is not None and len(t_eval) > 0:
        solve_span_min = min(actual_tau_min, np.min(t_eval))
        solve_span_max = max(actual_tau_max, np.max(t_eval))
    
    # phi_10.py 使用了事件
    events_to_pass = event_reached_phi_max_original

    sol = solve_ivp(
        pde_system_phi10_formulation,
        [solve_span_min, solve_span_max], # 求解区间
        y0,
        method='RK45', # 与 phi_10.py 一致
        t_eval=t_eval, # 在这些特定时间点输出结果
        args=(Pe1_eff, Pe2_eff, xi_grid, N_points),
        events=events_to_pass,
        dense_output=True # 与 phi_10.py 一致
    )

    phi1_solutions = sol.y[:N_points, :]
    phi2_solutions = sol.y[N_points:, :]

    # 新增: 限制 phi1 和 phi2 为非负值，以保证物理意义
    # 这会影响返回的 phi1 和 phi2 以及基于它们的计算（如偏析量和空隙率）
    # 这是期望的行为，因为负的体积分数没有物理意义。
    phi1_solutions = np.maximum(0, phi1_solutions)
    phi2_solutions = np.maximum(0, phi2_solutions)

    # 计算空隙体积分数
    if phi1_solutions.size > 0 and phi2_solutions.size > 0: # 确保有解
        phi_total_solids = phi1_solutions + phi2_solutions # 使用已修正的 phi1, phi2
        phi_void_solutions = 1.0 - phi_total_solids
        phi_void_solutions = np.maximum(0, phi_void_solutions) # 再次确保非负 (例如，如果 phi1+phi2 > 1)
        # 由于 phi1 和 phi2 已被修正为非负，phi_total_solids >= 0。
        # 因此 1.0 - phi_total_solids <= 1.0。
        # 所以 phi_void_solutions 也会 <= 1.0。
    else: # 如果没有解，则返回空的或形状匹配的nan数组
        # 确保返回的 phi_void_solutions 形状与期望一致，即使是空的
        # 如果 sol.t 为空，则 sol.y 也将是 (2*N_points, 0) 的形状
        # 那么 phi1_solutions 和 phi2_solutions 都是 (N_points, 0)
        phi_void_solutions = np.empty((N_points, phi1_solutions.shape[1]))
        phi_void_solutions[:] = np.nan

    # 计算偏析量
    stratification_phi2_percent = np.nan # 默认为nan
    if sol.t.size > 0 and phi2_solutions.shape[1] > 0:
         # 传递 phi2_solutions 本身，让 calculate_stratification 处理最后一个时间点
        stratification_phi2_percent = calculate_stratification(phi2_solutions, phi2_0_array, N_points)
    
    return (
        sol.t, phi1_solutions, phi2_solutions, phi_void_solutions, xi_grid, 
        stratification_phi2_percent, Pe1_eff, Pe2_eff, sol.message, sol.status, E_parameter
    )

# 注意：原始的 if __name__ == '__main__': 块将被移除，因为它属于测试/示例运行的范畴，
# 现在这个文件是模块的一部分。 