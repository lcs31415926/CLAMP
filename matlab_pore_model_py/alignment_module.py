import numpy as np

def calculate_alignment_and_contact_angle(delta_DA0=0.25, delta_DA_inf=0.37, k_DA=0.1, t0=6.5, x_range_max=2, x_step=0.1, theta_HI=30, theta_HO=120):
    """
    计算离子omer的分子排列 (delta_DA) 和接触角 (theta) 随位置的变化。
    对应 MATLAB 脚本 alignment_submodel.m 的核心计算逻辑。

    参数:
    delta_DA0 (float): 初始排列值
    delta_DA_inf (float): 无限厚膜的排列值
    k_DA (float): 衰减常数
    t0 (float): 初始膜厚度 (nm)
    x_range_max (float): 位置 x 的最大值
    x_step (float): 位置 x 的步长
    theta_HI (float): 接触角参数 HI (度)
    theta_HO (float): 接触角参数 HO (度)

    返回:
    tuple: (x, delta_DA, theta, t)
        x (np.array): 位置数组
        delta_DA (np.array): 离子omer分子排列数组
        theta (np.array): 接触角数组 (度)
        t (np.array): 膜厚度数组 (nm)
    """
    # 膜厚度变化位置
    x = np.arange(0, x_range_max + x_step, x_step) 
    # 膜厚度随位置变化
    t = t0 * np.exp(-k_DA * x)

    # 离子omer分子排列
    # 注意 MATLAB 中 t-t0，这里根据公式推测应与x或某种形式的厚度变化有关
    # 原始MATLAB公式: delta_DA = delta_DA_inf + (delta_DA0 - delta_DA_inf) .* exp(-(k_DA * (t - t0)));
    # 如果 t 是厚度，t0是初始厚度，(t-t0) 可能是厚度差。 
    # 但 MATLAB 脚本中 t 也是一个随 x 变化的数组。
    # 假设 MATLAB 的 (t - t0) 是希望表达从初始状态开始的某种变化量，
    # 并且这种变化是由 k_DA * x 驱动的，则可以简化或重新解释。
    # 或者，它可能是一个笔误，直接使用 x 或者与 x 相关的项更为合理。
    # 查阅 龟裂.md 中alignment_submodel的描述，似乎是 exp(-k_DA * x) 或 exp(-k_DA * (normalized_distance_or_time))
    # 这里我们尝试一种解释：exp(-(k_DA * x))，因为 t 本身就是 exp(-k_DA * x) 的函数
    # 另一种解释是 MATLAB 代码中的 (t - t0) 实际上是想用 x (位置) 作为衰减因子。
    # 我们将遵循公式的形式，但需要注意 t 的定义。
    # delta_DA = delta_DA_inf + (delta_DA0 - delta_DA_inf) * np.exp(-(k_DA * (t - t0))) 
    # 重新审视公式，如果 t 是当前厚度，t0 是初始厚度，那 (t-t0) 会是负值，这会导致 exp 的参数为正，可能不是预期的衰减。
    # 另一个角度，如果模型是关于距离表面的衰减，那么使用x更直接：
    delta_DA = delta_DA_inf + (delta_DA0 - delta_DA_inf) * np.exp(-k_DA * x) 

    # 计算接触角
    # np.cos 和 np.acosd 需要角度参数为弧度，但这里MATLAB用cosd/acosd，输入输出是度
    # Python 中 np.cos/np.sin/np.tan 输入为弧度，np.rad2deg, np.deg2rad 用于转换
    cos_theta_HI_rad = np.cos(np.deg2rad(theta_HI))
    cos_theta_HO_rad = np.cos(np.deg2rad(theta_HO))
    
    # delta_DA * (cos(theta_HI) - cos(theta_HO)) + cos(theta_HO) 的结果是 cos(theta) 的值
    cos_theta_val = delta_DA * (cos_theta_HI_rad - cos_theta_HO_rad) + cos_theta_HO_rad
    # 防止 acos 的参数超出 [-1, 1] 范围，这可能由于浮点数精度问题导致
    cos_theta_val = np.clip(cos_theta_val, -1.0, 1.0)
    theta_rad = np.arccos(cos_theta_val)
    theta = np.rad2deg(theta_rad) # 将结果从弧度转换回度

    return x, delta_DA, theta, t 