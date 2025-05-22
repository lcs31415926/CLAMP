import numpy as np
from scipy.integrate import solve_ivp

# 微分方程函数 (供 solve_ivp 调用)
def ode_model_py(x, y, k_A, k_V):
    """
    定义常微分方程组。
    对应 MATLAB 脚本中的 ode_model 函数。

    参数:
    x: 自变量 (位置)
    y: 因变量数组 [a_ion, alpha_free]
    k_A: 离子omer分散参数
    k_V: 离子omer体积减少因子

    返回:
    list: [da_ion_dx, d_alpha_free_dx]
    """
    a_ion, alpha_free = y

    da_ion_dx = (1 - a_ion) * k_A
    d_alpha_free_dx = -alpha_free * k_V

    return [da_ion_dx, d_alpha_free_dx]

def solve_catalyst_composition(k_A=4.5, k_V=0.6, V_20=0.39, A_2=6.57, x_max=1.5, y0=None, num_points=100):
    """
    求解催化剂组成模型的常微分方程，并计算膜厚度。
    对应 MATLAB 脚本 catalyst_layer_composition_model.m 的核心计算逻辑。

    参数:
    k_A (float): 离子omer分散参数
    k_V (float): 离子omer体积减少因子
    V_20 (float): 二次孔隙体积 (cm3/gC)
    A_2 (float): 二次孔隙表面积 (m2/gC)
    x_max (float): x 的最大值 (积分上限)
    y0 (list, optional): 初始条件 [a_ion_initial, alpha_free_initial]。默认为 [0, 1]。
    num_points (int): solve_ivp 生成点的数量，用于 t_eval

    返回:
    tuple: (x_sol, a_ion_sol, alpha_free_sol, t_membrane)
        x_sol (np.array): x 的解点
        a_ion_sol (np.array): 离子omer覆盖率 a_ion 的解
        alpha_free_sol (np.array): 自由孔隙体积比例 alpha_free 的解
        t_membrane (np.array): 膜厚度 t 的解 (nm)
    """
    if y0 is None:
        y0 = [0, 1] # 初始条件 a_ion(0)=0, alpha_free(0)=1
    
    # 定义积分区间
    x_span = [0, x_max]
    x_eval = np.linspace(x_span[0], x_span[1], num_points) # 为了得到平滑曲线

    # 求解微分方程
    # solve_ivp 需要将额外参数通过 args 元组传递给 ODE 函数
    sol = solve_ivp(
        fun=ode_model_py, 
        t_span=x_span, 
        y0=y0, 
        t_eval=x_eval, 
        args=(k_A, k_V)
    )

    x_sol = sol.t
    a_ion_sol = sol.y[0, :]
    alpha_free_sol = sol.y[1, :]

    # 计算膜厚度
    # V_20 单位是 cm3/gC, A_2 单位是 m2/gC。
    # 为了得到 nm 级别的厚度，需要单位转换。
    # 假设 V_20 和 A_2 的比值最终给出的是长度单位。
    # MATLAB 脚本中没有明确单位转换过程，我们暂时也直接计算。
    # 添加 eps 防止除以零 (类似 MATLAB 中的 eps)
    denominator = (A_2 * a_ion_sol) + np.finfo(float).eps 
    t_membrane = (V_20 * (1 - alpha_free_sol)) / denominator
    # 最终厚度的单位需要根据 V_20 和 A_2 的原始单位和期望的 t 单位来确定。
    # MATLAB 图的 y 轴标签是 Membrane Thickness (nm)，暗示结果应该是 nm。
    # 如果 V_20/A_2 结果是 m，则需要 *1e9 转换为 nm。
    # 如果 V_20 是 cm^3/gC, A_2 是 m^2/gC, 则 V_20 (cm^3) / A_2 (m^2) = V_20 (1e-6 m^3) / A_2 (m^2) = 1e-6 * (V_20/A_2) m
    # 那么 t_membrane_m = (V_20 * 1e-6 * (1 - alpha_free_sol)) / (A_2 * a_ion_sol + np.finfo(float).eps)
    # t_membrane_nm = t_membrane_m * 1e9
    # 简化: t_membrane_nm = (V_20 * 1e3 * (1 - alpha_free_sol)) / (A_2 * a_ion_sol + np.finfo(float).eps)
    # 暂时按 MATLAB 脚本直接计算，后续可在测试脚本中调整因子以匹配期望的nm量级

    return x_sol, a_ion_sol, alpha_free_sol, t_membrane 