import numpy as np
from scipy.stats import lognorm, norm

def calculate_pore_and_angle_distributions(
    r_mu1=3, sigma_r1=0.8, r_mu2=80, sigma_r2=1.2, 
    phi1=0.21, phi2=0.79, 
    a_Pt_initial=0.4, a_ion_initial=0.6, 
    theta_Pt_deg=10, theta_ion_deg=90, 
    current_r_for_cad=20, 
    pore_radius_min=0.1, pore_radius_max=200, pore_radius_step=0.1,
    contact_angle_min=0, contact_angle_max=180, contact_angle_step=1,
    cad_sigma_assumed=5
):
    """
    计算孔径分布、接触角分布和它们的联合分布。
    对应 MATLAB 脚本 pore_angle_distribution_model.m 的核心计算逻辑。

    参数:
    r_mu1, r_mu2 (float): 一次/二次孔径的几何平均值 (nm) (对应 log(r) 的均值)
    sigma_r1, sigma_r2 (float): 一次/二次孔径分布的对数标准差
    phi1, phi2 (float): 一次/二次孔隙体积分数
    a_Pt_initial, a_ion_initial (float): Pt/离子omer的初始表面覆盖率
    theta_Pt_deg, theta_ion_deg (float): Pt/离子omer的接触角 (度)
    current_r_for_cad (float): 用于计算有效接触角的当前孔径 (nm)
    pore_radius_min, max, step: 孔径范围参数 (nm)
    contact_angle_min, max, step: 接触角范围参数 (度)
    cad_sigma_assumed (float): 假定的接触角分布的标准差 (度)

    返回:
    dict: 包含所有计算结果的字典，例如：
        'pore_radius_array', 'psd1_normalized', 'psd2_normalized', 
        'pore_volume1', 'pore_volume2', 'contact_angle_array', 
        'theta_pore_effective_deg', 'cad_normalized', 
        'R_meshgrid', 'C_meshgrid', 'Z_joint_distribution'
    """
    results = {}

    # --- 孔径分布 (PSD) ---
    pore_radius_array = np.arange(pore_radius_min, pore_radius_max + pore_radius_step, pore_radius_step)
    results['pore_radius_array'] = pore_radius_array

    # MATLAB lognpdf(X, MU, SIGMA) 中的 MU 和 SIGMA 是 log(X) 的均值和标准差。
    # scipy.stats.lognorm(s, loc, scale) 中 s=SIGMA, scale=exp(MU)。loc 通常为0。
    # MU_matlab = log(r_mu_geometric_mean) -> exp(MU_matlab) = r_mu_geometric_mean
    # sigma_matlab = sigma_log_transformed_data
    
    # 一次孔径分布
    # mu_log1 = np.log(r_mu1) # r_mu1已经是几何平均值，即 exp(mu_log_space)
    psd1 = lognorm.pdf(pore_radius_array, s=sigma_r1, loc=0, scale=r_mu1)
    psd1_normalized = psd1 / (np.sum(psd1) * pore_radius_step) # 乘以步长以近似积分
    results['psd1_raw'] = psd1
    results['psd1_normalized'] = psd1_normalized

    # 二次孔径分布
    # mu_log2 = np.log(r_mu2)
    psd2 = lognorm.pdf(pore_radius_array, s=sigma_r2, loc=0, scale=r_mu2)
    psd2_normalized = psd2 / (np.sum(psd2) * pore_radius_step)
    results['psd2_raw'] = psd2
    results['psd2_normalized'] = psd2_normalized

    # 孔隙体积分数 (按分布加权)
    pore_volume1 = phi1 * psd1_normalized
    pore_volume2 = phi2 * psd2_normalized
    results['pore_volume1'] = pore_volume1
    results['pore_volume2'] = pore_volume2

    # --- 接触角分布 (CAD) ---
    contact_angle_array = np.arange(contact_angle_min, contact_angle_max + contact_angle_step, contact_angle_step)
    results['contact_angle_array'] = contact_angle_array

    # 当前孔径 r (来自参数 current_r_for_cad)
    r_current = current_r_for_cad 

    # 调整后的颗粒半径 (不能为负)
    # MATLAB: rad_Pt = r - 1.5; rad_ion = r - 5;
    # 这些常数 (1.5, 5) 可能是特定颗粒的尺寸或层厚度，这里硬编码了
    rad_Pt_eff = max(0, r_current - 1.5)
    rad_ion_eff = max(0, r_current - 5)

    # Pt 纳米颗粒和离子omer的有效覆盖率 (基于当前孔径 r_current)
    # 如果 r_current 为0或很小，可能导致除零。添加一个小的eps。
    a_Pt_effective = a_Pt_initial * (rad_Pt_eff / (r_current + np.finfo(float).eps))
    a_ion_effective = a_ion_initial * (rad_ion_eff / (r_current + np.finfo(float).eps))
    results['a_Pt_effective'] = a_Pt_effective
    results['a_ion_effective'] = a_ion_effective

    # 有效接触角 (Cassie-Baxter like)
    cos_theta_Pt_rad = np.cos(np.deg2rad(theta_Pt_deg))
    cos_theta_ion_rad = np.cos(np.deg2rad(theta_ion_deg))
    
    cos_theta_pore_effective = a_Pt_effective * cos_theta_Pt_rad + a_ion_effective * cos_theta_ion_rad
    # 确保值在 acos 的有效范围内 [-1, 1]
    cos_theta_pore_effective = np.clip(cos_theta_pore_effective, -1.0, 1.0)
    theta_pore_effective_rad = np.arccos(cos_theta_pore_effective)
    theta_pore_effective_deg = np.rad2deg(theta_pore_effective_rad)
    results['theta_pore_effective_deg'] = theta_pore_effective_deg

    # 接触角分布 (假设以有效接触角为均值，标准差由参数 cad_sigma_assumed 定义)
    cad = norm.pdf(contact_angle_array, loc=theta_pore_effective_deg, scale=cad_sigma_assumed)
    cad_normalized = cad / (np.sum(cad) * contact_angle_step) # 乘以步长以近似积分
    results['cad_raw'] = cad
    results['cad_normalized'] = cad_normalized

    # --- 联合分布 ---
    R_meshgrid, C_meshgrid = np.meshgrid(pore_radius_array, contact_angle_array)
    results['R_meshgrid'] = R_meshgrid
    results['C_meshgrid'] = C_meshgrid

    Z_joint_distribution = np.zeros((len(contact_angle_array), len(pore_radius_array)))
    # MATLAB: Z(j, i) = pore_volume1(i) * cad(j) + pore_volume2(i) * cad(j);
    # 使用归一化后的 cad_normalized 和 pore_volume (已经是加权的)
    for i in range(len(pore_radius_array)):
        for j in range(len(contact_angle_array)):
            # Z_joint_distribution[j, i] = (pore_volume1[i] + pore_volume2[i]) * cad_normalized[j]
            # 或者，如果认为 cad 的计算已经基于特定孔径 r_current, 那么它对所有孔径i都一样
            # P(radius_i, angle_j) = P(radius_i) * P(angle_j | radius_i)
            # P(radius_i) 是 (psd1_normalized[i]*phi1 + psd2_normalized[i]*phi2)
            # P(angle_j | radius_i) 是 cad_normalized[j] (假设它对于不同的 i 是相同的，因为 cad 是基于一个 current_r_for_cad 计算的)
            # 所以 Z(j,i) = (phi1 * psd1_normalized[i] + phi2 * psd2_normalized[i]) * cad_normalized[j]
            # pore_volume1[i] = phi1 * psd1_normalized[i]
            # pore_volume2[i] = phi2 * psd2_normalized[i]
            Z_joint_distribution[j, i] = (pore_volume1[i] + pore_volume2[i]) * cad_normalized[j]

    results['Z_joint_distribution'] = Z_joint_distribution
    
    return results 