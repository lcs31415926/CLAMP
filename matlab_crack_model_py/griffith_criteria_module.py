import numpy as np

def griffith_criteria(E, gamma, a):
    """
    计算裂纹生成的临界应力（Griffith 模型）。

    参数:
    E (float): 杨氏模量 (Pa)
    gamma (float): 表面能 (J/m^2)
    a (float): 裂纹半长 (m)

    返回:
    float: 临界断裂应力 (Pa)
    """
    # 计算临界断裂应力
    sigma_c = np.sqrt((2 * E * gamma) / (np.pi * a))
    return sigma_c 