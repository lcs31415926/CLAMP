import numpy as np

def tirumkudulu_russel_model(r, h_layer, G, M, phi_rcp, L_surface_tension):
    """
    计算裂纹生成的临界干燥应力 (Tirumkudulu & Russel 模型)。
    参数 L_surface_tension 对应 MATLAB 版本中的 L (液-气界面表面张力)。

    参数:
    r (float): 孔隙半径 (m)
    h_layer (float): 涂层厚度 (m)
    G (float): 剪切模量 (Pa)
    M (float): 压缩模量 (Pa)
    phi_rcp (float): 随机密堆积体积分数
    L_surface_tension (float): 液-气界面表面张力 (N/m) (在MATLAB脚本中参数名为L)

    返回:
    float: 临界干燥应力 (Pa)
    """
    # 计算临界干燥应力
    term1 = (2 * r / h_layer)**(2/3)
    term2 = (G * M * phi_rcp * r / (2 * L_surface_tension))**(1/3)
    term3 = (2 * L_surface_tension) / r 
    sigma_crit = 0.1877 * term1 * term2 * term3
    return sigma_crit 