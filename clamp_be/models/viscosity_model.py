import math

def calculate_volume_fractions(recipe_params: dict) -> dict:
    """
    计算各种组分的体积分数和离聚物的总质量浓度。

    参数:
    recipe_params (dict): 包含以下键值对的字典
        "w_catalyst": 催化剂质量分数 (例如 0.7),
        "rho_catalyst": 催化剂颗粒密度 (g/cm3),
        "w_carbon": 碳载体质量分数 (例如 0.1), (如果适用, 否则为0)
        "rho_carbon": 碳载体颗粒密度 (g/cm3), (如果适用)
        "w_ionomer": 离聚物质量分数 (例如 0.2),
        "rho_ionomer": 离聚物密度 (g/cm3),
        "rho_solvent": 溶剂密度 (g/cm3)
        (注意: 溶剂质量分数 w_solvent = 1 - w_catalyst - w_carbon - w_ionomer)

    返回:
    dict: 包含以下键值对的字典
        "phi_catalyst": 催化剂体积分数,
        "phi_carbon": 碳载体体积分数,
        "phi_ionomer_total": 总离聚物体积分数,
        "phi_solid": 固相总体积分数 (催化剂 + 碳 + 离聚物),
        "phi_solvent": 溶剂体积分数,
        "C_ionomer_total": 总离聚物质量浓度 (g/cm3, 基于总浆料体积)
    """
    w_c = recipe_params["w_catalyst"]
    rho_c = recipe_params["rho_catalyst"]
    w_s = recipe_params.get("w_carbon", 0) # 使用get以防不存在该键
    rho_s = recipe_params.get("rho_carbon", 1.0) # 默认为1防止除零，实际应提供
    w_i = recipe_params["w_ionomer"]
    rho_i = recipe_params["rho_ionomer"]
    rho_sol = recipe_params["rho_solvent"]

    w_solvent = 1 - w_c - w_s - w_i
    if w_solvent < 0:
        raise ValueError("总质量分数超过100%，请检查输入。")

    # 单个组分的体积，假设总质量为 M_total
    # V_x = (w_x * M_total) / rho_x
    # V_total = V_c + V_s + V_i + V_solvent
    # phi_x = V_x / V_total
    # 为了消除M_total，我们计算每单位质量浆料的体积
    # v_x = w_x / rho_x (比容)
    # v_total_slurry = w_c/rho_c + w_s/rho_s + w_i/rho_i + w_solvent/rho_sol

    v_total_slurry = (w_c / rho_c) + \
                       (w_s / rho_s if w_s > 0 else 0) + \
                       (w_i / rho_i) + \
                       (w_solvent / rho_sol)
    
    if v_total_slurry == 0:
        raise ValueError("总比容为零，请检查密度和质量分数输入。")

    phi_c = (w_c / rho_c) / v_total_slurry
    phi_s = (w_s / rho_s) / v_total_slurry if w_s > 0 else 0
    phi_i_total = (w_i / rho_i) / v_total_slurry
    phi_solvent = (w_solvent / rho_sol) / v_total_slurry

    phi_solid = phi_c + phi_s + phi_i_total # 这里的固相包括了离聚物

    # 总离聚物质量浓度 C_i_total (g/cm3 of total slurry)
    # C_i_total = mass_ionomer / V_total_slurry
    # mass_ionomer = w_i * M_total
    # V_total_slurry = M_total * v_total_slurry
    # C_i_total = w_i / v_total_slurry
    C_i_total = w_i / v_total_slurry if v_total_slurry !=0 else 0


    return {
        "phi_catalyst": phi_c,
        "phi_carbon": phi_s,
        "phi_ionomer_total": phi_i_total,
        "phi_solid": phi_solid, # 定义为催化剂、碳、离聚物的总体积分数
        "phi_solvent": phi_solvent,
        "C_ionomer_total": C_i_total
    }

def calculate_effective_solvent_viscosity(viscosity_params: dict, recipe_params: dict) -> float:
    """
    计算有效溶剂粘度。

    参数:
    viscosity_params (dict): 包含以下键值对
        "C_non_ads": 非吸附离聚物质量浓度 (g/cm3),
        "k_H": Huggins常数,
        "intrinsic_viscosity_ionomer": 离聚物特性粘度 [η] (cm3/g)
    recipe_params (dict): 包含 "eta_solvent" (溶剂粘度 Pa.s)

    返回:
    float: 有效溶剂粘度 (Pa.s)
    """
    C_non_ads = viscosity_params["C_non_ads"]
    k_H = viscosity_params["k_H"]
    eta_intrinsic_ionomer = viscosity_params["intrinsic_viscosity_ionomer"]
    eta_sol = recipe_params["eta_solvent"]

    # Huggins方程: eta_eff_sol = eta_sol * (1 + [eta]*C + k_H*[eta]^2*C^2)
    eta_eff_sol = eta_sol * (1 + eta_intrinsic_ionomer * C_non_ads + \
                             k_H * (eta_intrinsic_ionomer**2) * (C_non_ads**2))
    return eta_eff_sol

def calculate_slurry_viscosity_kd(eta_eff_sol: float, 
                                  solid_eff_vol_fraction: float, 
                                  max_packing_fraction_phi_m: float, 
                                  particle_intrinsic_viscosity: float) -> float:
    """
    使用Krieger-Dougherty模型计算浆料粘度。

    参数:
    eta_eff_sol (float): 有效溶剂粘度 (Pa.s)
    solid_eff_vol_fraction (float): 有效固相体积分数 (phi_solid,eff)
    max_packing_fraction_phi_m (float): 最大堆积体积分数 (phi_m)
    particle_intrinsic_viscosity (float): 颗粒的特性粘度 ([eta]_particle, 对于球形颗粒通常为2.5)

    返回:
    float: 浆料粘度 (Pa.s)
    """
    if max_packing_fraction_phi_m == 0:
        raise ValueError("最大堆积体积分数 phi_m 不能为零。")
    if solid_eff_vol_fraction >= max_packing_fraction_phi_m:
        # 当体积分数接近或超过最大堆积分数时，粘度趋于无穷大
        # 可以返回一个非常大的数值或者抛出错误
        print("警告: 固相体积分数接近或超过最大堆积分数，粘度可能非常高。")
        return float('inf')
        
    term = (1 - solid_eff_vol_fraction / max_packing_fraction_phi_m)
    eta_slurry = eta_eff_sol * (term ** (-particle_intrinsic_viscosity * max_packing_fraction_phi_m))
    return eta_slurry

def calculate_slurry_viscosity_quemada(eta_eff_sol: float, 
                                       solid_eff_vol_fraction: float, 
                                       max_packing_fraction_phi_m: float) -> float:
    """
    使用Quemada模型计算浆料粘度。

    参数:
    eta_eff_sol (float): 有效溶剂粘度 (Pa.s)
    solid_eff_vol_fraction (float): 有效固相体积分数 (phi_solid,eff)
    max_packing_fraction_phi_m (float): 最大堆积体积分数 (phi_m)

    返回:
    float: 浆料粘度 (Pa.s)
    """
    if max_packing_fraction_phi_m == 0:
        raise ValueError("最大堆积体积分数 phi_m 不能为零。")
    if solid_eff_vol_fraction >= max_packing_fraction_phi_m:
        print("警告: 固相体积分数接近或超过最大堆积分数，粘度可能非常高。")
        return float('inf')

    term = (1 - solid_eff_vol_fraction / max_packing_fraction_phi_m)
    eta_slurry = eta_eff_sol * (term ** (-2))
    return eta_slurry

def get_initial_apparent_viscosity(recipe_params: dict, p_viscosity: dict) -> float:
    """
    计算浆料的初始表观粘度。

    参数:
    recipe_params (dict): 浆料配方参数 (见 calculate_volume_fractions)
    p_viscosity (dict): 粘度模型相关参数，包含:
        "C_non_ads": 非吸附离聚物浓度 (g/cm3),
        "k_H": Huggins常数,
        "intrinsic_viscosity_ionomer": 离聚物特性粘度 [η] (cm3/g),
        "phi_m_kd": K-D模型中的最大堆积体积分数 (如果选择K-D模型),
        "intrinsic_viscosity_particle_kd": K-D模型中颗粒的特性粘度 (如果选择K-D模型),
        "phi_m_quemada": Quemada模型中的最大堆积体积分数 (如果选择Quemada模型, 可以与phi_m_kd相同或不同),
        "model_choice": "KriegerDougherty" 或 "Quemada"

    返回:
    float: 浆料初始表观粘度 (eta_app_0, Pa.s)
    """
    print("开始计算初始表观粘度...")
    print(f"配方参数: {recipe_params}")
    print(f"粘度参数: {p_viscosity}")

    # 步骤1: 计算体积分数
    volume_fractions = calculate_volume_fractions(recipe_params)
    phi_solid_eff = volume_fractions["phi_solid"] # 这里的phi_solid包含了所有非溶剂组分
    print(f"计算得到的体积分数: {volume_fractions}")

    # 步骤2: 计算有效溶剂粘度
    eta_eff_sol = calculate_effective_solvent_viscosity(p_viscosity, recipe_params)
    print(f"计算得到的有效溶剂粘度: {eta_eff_sol:.4e} Pa.s")

    # 步骤3: 根据选择的模型计算浆料粘度
    model_choice = p_viscosity.get("model_choice", "KriegerDougherty") # 默认为K-D模型
    print(f"选择的粘度模型: {model_choice}")

    if model_choice == "KriegerDougherty":
        phi_m = p_viscosity["phi_m_kd"]
        intrinsic_viscosity_particle = p_viscosity["intrinsic_viscosity_particle_kd"]
        eta_app_0 = calculate_slurry_viscosity_kd(
            eta_eff_sol,
            phi_solid_eff,
            phi_m,
            intrinsic_viscosity_particle
        )
    elif model_choice == "Quemada":
        # Quemada模型通常也需要phi_m，这里假设与K-D的phi_m相同或在P_Viscosity中单独指定
        phi_m = p_viscosity.get("phi_m_quemada", p_viscosity.get("phi_m_kd")) 
        if phi_m is None:
            raise ValueError("Quemada模型需要 phi_m_quemada 或 phi_m_kd 参数。")
        eta_app_0 = calculate_slurry_viscosity_quemada(
            eta_eff_sol,
            phi_solid_eff,
            phi_m
        )
    else:
        raise ValueError(f"未知的粘度模型选择: {model_choice}")

    print(f"计算得到的初始表观粘度 (eta_app_0): {eta_app_0:.4e} Pa.s")
    return eta_app_0

if __name__ == '__main__':
    print("开始粘度模型测试...")

    # 示例输入参数
    SlurryRecipeParams_example = {
        "w_catalyst": 0.6,  # 催化剂质量分数
        "rho_catalyst": 4.0, # 催化剂颗粒密度 (g/cm3)
        "w_carbon": 0.1,    # 碳载体质量分数
        "rho_carbon": 2.0,  # 碳载体颗粒密度 (g/cm3)
        "w_ionomer": 0.3,   # 离聚物质量分数 (0.3使其总和为1)
        "rho_ionomer": 1.8, # 离聚物密度 (g/cm3)
        "rho_solvent": 0.9, # 溶剂密度 (g/cm3)
        "eta_solvent": 0.001 # 溶剂粘度 (Pa.s)
    }
    # 确保质量分数总和为1
    total_solids_w = SlurryRecipeParams_example["w_catalyst"] + \
                     SlurryRecipeParams_example["w_carbon"] + \
                     SlurryRecipeParams_example["w_ionomer"]
    if total_solids_w > 1.0:
        print(f"警告: 固体质量分数总和 {total_solids_w} 大于1，请调整。")
        # 简单调整，按比例减少，使总和为1（这只是一个示例处理方式）
        SlurryRecipeParams_example["w_catalyst"] /= total_solids_w
        SlurryRecipeParams_example["w_carbon"] /= total_solids_w
        SlurryRecipeParams_example["w_ionomer"] /= total_solids_w
        print(f"调整后的配方: w_c={SlurryRecipeParams_example['w_catalyst']:.2f}, w_s={SlurryRecipeParams_example['w_carbon']:.2f}, w_i={SlurryRecipeParams_example['w_ionomer']:.2f}")


    P_Viscosity_example_kd = {
        "C_non_ads": 0.05, # 非吸附离聚物浓度 (g/cm3) - 假设值
        "k_H": 0.35,       # Huggins常数
        "intrinsic_viscosity_ionomer": 20, # 离聚物特性粘度 (cm3/g)
        "phi_m_kd": 0.60,    # K-D模型中的最大堆积体积分数
        "intrinsic_viscosity_particle_kd": 2.5, # K-D颗粒特性粘度
        "model_choice": "KriegerDougherty"
    }

    P_Viscosity_example_quemada = {
        "C_non_ads": 0.05,
        "k_H": 0.35,
        "intrinsic_viscosity_ionomer": 20,
        "phi_m_quemada": 0.60, # Quemada模型最大堆积体积分数
        "model_choice": "Quemada"
    }
    
    print("\n--- 测试Krieger-Dougherty模型 ---")
    try:
        eta_app_0_kd = get_initial_apparent_viscosity(SlurryRecipeParams_example, P_Viscosity_example_kd)
        print(f"K-D模型 - 最终计算得到的初始表观粘度: {eta_app_0_kd:.4e} Pa.s")
        assert eta_app_0_kd > 0, "K-D粘度计算结果应为正值"
        if math.isinf(eta_app_0_kd):
             print("K-D粘度计算结果为无穷大，可能固相体积分数过高。")
    except ValueError as e:
        print(f"K-D模型计算出错: {e}")
    except Exception as e:
        print(f"K-D模型发生未知错误: {e}")

    print("\n--- 测试Quemada模型 ---")
    try:
        # 为了测试Quemada，可以复用部分P_Viscosity_example_kd的参数，但要明确phi_m
        P_Viscosity_example_quemada_test = P_Viscosity_example_kd.copy()
        P_Viscosity_example_quemada_test["model_choice"] = "Quemada"
        # phi_m_quemada可以与phi_m_kd相同或不同，这里假设相同
        P_Viscosity_example_quemada_test["phi_m_quemada"] = P_Viscosity_example_kd["phi_m_kd"] 

        eta_app_0_quemada = get_initial_apparent_viscosity(SlurryRecipeParams_example, P_Viscosity_example_quemada_test)
        print(f"Quemada模型 - 最终计算得到的初始表观粘度: {eta_app_0_quemada:.4e} Pa.s")
        assert eta_app_0_quemada > 0, "Quemada粘度计算结果应为正值"
        if math.isinf(eta_app_0_quemada):
             print("Quemada粘度计算结果为无穷大，可能固相体积分数过高。")
    except ValueError as e:
        print(f"Quemada模型计算出错: {e}")
    except Exception as e:
        print(f"Quemada模型发生未知错误: {e}")

    print("\n--- 测试高固含量情况 (K-D) ---")
    SlurryRecipeParams_high_solid = {
        "w_catalyst": 0.75, # 更高催化剂含量
        "rho_catalyst": 4.0,
        "w_carbon": 0.05,  
        "rho_carbon": 2.0,
        "w_ionomer": 0.20, # 总固含量为1.0
        "rho_ionomer": 1.8,
        "rho_solvent": 0.9, 
        "eta_solvent": 0.001 
    }
    # 此时 w_solvent = 1 - 0.75 - 0.05 - 0.20 = 0.0. 这是极限情况。
    # 如果 w_solvent < 0, calculate_volume_fractions 会报错
    
    try:
        eta_app_0_kd_high = get_initial_apparent_viscosity(SlurryRecipeParams_high_solid, P_Viscosity_example_kd)
        print(f"K-D模型 (高固含量) - 最终计算得到的初始表观粘度: {eta_app_0_kd_high:.4e} Pa.s")
        assert eta_app_0_kd_high > 0, "K-D粘度计算结果应为正值"
    except ValueError as e:
        print(f"K-D模型 (高固含量) 计算出错: {e}") # 预期可能会在这里出错，例如phi_solid > phi_m
    except Exception as e:
        print(f"K-D模型 (高固含量) 发生未知错误: {e}")

    print("\n粘度模型测试结束.")
