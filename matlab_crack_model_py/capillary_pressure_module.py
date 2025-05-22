import numpy as np

def capillary_pressure(L, theta, r_p):
    """
    计算毛细管压力。

    参数:
    L (float): 表面张力 (N/m)
    theta (float): 接触角 (度)
    r_p (float): 孔隙半径 (m)

    返回:
    float: 毛细管压力 (Pa)
    """
    # 将接触角从度转换为弧度
    theta_rad = np.deg2rad(theta)
    
    # 计算毛细压力
    p_cap = (2 * L * np.cos(theta_rad)) / r_p
    return p_cap 