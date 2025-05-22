import numpy as np

def stress_function(d_t, delta_m_t, L_param, E_s, b, h_s, h_c_t, v_s):
    """
    计算随时间变化的应力 sigma(t)。
    注意：MATLAB 版本将 g 作为参数，此 Python 版本使用内部常数 g。
    L_param 对应 MATLAB 中的 L (梁的长度)。

    参数:
    d_t (float): 梁中点处的挠度 (m)
    delta_m_t (float): 梁的总质量变化 (kg)
    L_param (float): 梁的长度 (m) (在MATLAB脚本中参数名为L)
    E_s (float): 基底的杨氏模量 (Pa)
    b (float): 梁的宽度 (m)
    h_s (float): 基底的厚度 (m)
    h_c_t (float): 涂层随时间变化的厚度 (m)
    v_s (float): 基底的泊松比

    返回:
    float: 随时间变化的应力 (Pa)
    """
    # 常数 g（重力加速度）
    g = 9.81  # m/s^2
    
    # 计算分子部分
    # 注意：MATLAB 版本中 L 是参数 L，这里改为 L_param 以避免与表面张力L混淆
    numerator = (d_t - (3 * delta_m_t * g * L_param**3 / (2 * E_s * b * h_s**3))) * E_s * h_s**3
    
    # 计算分母部分
    denominator = 3 * h_c_t * L_param**2 * (h_s + h_c_t) * (1 - v_s)
    
    # 计算应力
    if denominator == 0:
        # 根据具体情况处理分母为零的情况，这里返回 np.inf 或抛出异常
        return np.inf 
    sigma_t = numerator / denominator
    return sigma_t 