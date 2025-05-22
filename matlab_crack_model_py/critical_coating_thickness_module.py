import numpy as np

def critical_coating_thickness(G, M, phi_rcp, L, p_cap_max):
    """
    计算临界涂层厚度。

    参数:
    G (float): 剪切模量 (Pa)
    M (float): 压缩模量 (Pa)
    phi_rcp (float): 随机密堆积体积分数
    L (float): 表面张力 (N/m)
    p_cap_max (float): 最大毛细管压力 (Pa) (应为负值)

    返回:
    float: 临界涂层厚度 (m)
    """
    # 计算临界涂层厚度
    term1 = (G * M * phi_rcp / (2 * L))**(1/2)
    # MATLAB 中的 p_cap_max 预期是正值，然后在公式中取负
    # Python 版本直接使用，调用时确保 p_cap_max 的符号正确
    term2 = (2 * L / (-p_cap_max if p_cap_max > 0 else p_cap_max))**(3/2) 
    CCT = 0.64 * term1 * term2
    return CCT 