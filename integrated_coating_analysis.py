import numpy as np
import matplotlib.pyplot as plt
import json
from matlab_crack_model_py.capillary_pressure_module import capillary_pressure
from matlab_crack_model_py.cracking_degree_module import get_cracking_degree_data_and_fits # 新增导入

# --- 辅助模块导入 (未来替换为实际模型导入) ---
# from matlab_high_speed_shearing_model_py.core_hss_solver import run_hss_simulation
# from matlab_drying_model_phi10_py.core_solver import run_drying_simulation
# from matlab_pore_model_py.main_pore_analyzer import analyze_pore_structure # 假设的
# from matlab_crack_model_py.coating_crack_analysis import run_crack_analysis # 假设的

# --- 桩函数定义 ---

def run_hss_simulation_stub(initial_slurry_properties, shearing_process_params):
    print("\n--- 调用 高速剪切模型 (桩函数) ---")
    print(f"HSS STUB: 接收到 initial_slurry_properties 键: {list(initial_slurry_properties.keys())}")
    print(f"HSS STUB: 接收到 shearing_process_params: {shearing_process_params}")
    
    # 模拟材料表面特性的传递和可能的变化
    output_surface_char = initial_slurry_properties.get('material_surface_char', {})
    # output_surface_char['shear_modification_tag'] = "processed_by_hss" # 示例：标记已被处理

    return {
        'sheared_slurry_state': {
            'effective_particle_size_distribution': np.array([[5.0, 0.8], [15.0, 0.2]]), # 示例: [直径_um, 比例]
            'slurry_rheology': {'viscosity_final': 0.1, 'yield_stress': 0.5}, # Pa.s, Pa
            'dispersion_quality_index': 0.85, # 0-1
            'temperature_final': shearing_process_params.get('temperature', 25) + 2, # 假设剪切升温2度
            'material_surface_char_modified': output_surface_char
        },
        'processing_summary_hss': {'status': 'HSS stub_success', 'simulated_time_ms': 150}
    }

def run_drying_simulation_stub(initial_coating_params, material_properties_drying, drying_process_params, numerical_params_drying):
    print("\n--- 调用 干燥与偏析模型 (桩函数) ---")
    print(f"DRYING STUB: 接收到 initial_coating_params: {initial_coating_params}")
    print(f"DRYING STUB: 接收到 material_properties_drying 键: {list(material_properties_drying.keys())}")
    # print(f"DRYING STUB: material_properties_drying['material_surface_char_effective']: {material_properties_drying.get('material_surface_char_effective')}")
    print(f"DRYING STUB: 接收到 drying_process_params: {drying_process_params}")
    print(f"DRYING STUB: 接收到 numerical_params_drying: {numerical_params_drying}")

    nz = numerical_params_drying.get('nz', 50)
    phi_s_final = np.linspace(0.6, 0.75, nz) # 模拟一个固相体积分数分布
    
    initial_thickness = initial_coating_params.get('wet_film_thickness', 100e-6)
    # 简化估算干膜厚度：湿膜厚度 * (phi1_initial + phi2_initial) / 平均最终固相含量
    # 这里用一个更简化的方式
    avg_solid_content_initial = initial_coating_params.get('component_concentrations_initial',{}).get('phi1_initial',0.2) + \
                                initial_coating_params.get('component_concentrations_initial',{}).get('phi2_initial',0.05)
    
    dry_film_thickness_simulated = initial_thickness * (avg_solid_content_initial / np.mean(phi_s_final)) if np.mean(phi_s_final) > 0 else initial_thickness * 0.3

    return {
        'dried_film_state': {
            'phi1_profile_final': np.linspace(0.4, 0.5, nz), # 示例数据
            'phi2_profile_final': np.linspace(0.2, 0.25, nz), # 示例数据
            'solid_fraction_profile_final': phi_s_final,
            'dry_film_thickness': dry_film_thickness_simulated, # m
            'stratification_metric_phi1': 0.15,
            'stratification_metric_phi2': -0.10,
            'residual_stress_profile_final': np.random.rand(nz) * 1e6 # Pa, 示例应力分布
        },
        'drying_dynamics_data': {
            'film_thickness_vs_time': np.array([initial_thickness, initial_thickness*0.8, initial_thickness*0.5, dry_film_thickness_simulated]),
        },
        'processing_summary_drying': {'status': 'Drying stub_success', 'total_drying_time_simulated': drying_process_params.get('drying_time_total',1800)}
    }

def analyze_pore_structure_stub(dried_film_input, particle_properties_pore, catalyst_info):
    print("\n--- 调用 孔隙结构模型 (桩函数) ---")
    print(f"PORE STUB: 接收到 dried_film_input 键: {list(dried_film_input.keys())}")
    print(f"PORE STUB: 接收到 particle_properties_pore: {particle_properties_pore}")
    print(f"PORE STUB: 接收到 catalyst_info: {catalyst_info}")

    solid_fraction_profile = dried_film_input.get('solid_fraction_profile', np.array([0.7]*10))
    porosity_profile = 1 - solid_fraction_profile
    average_porosity = np.mean(porosity_profile)
    
    # 基于平均孔隙率和输入的粒径（简化）估算平均孔径
    # 假设颗粒为球形，孔隙为颗粒间的空隙
    # mean_particle_size_um = np.mean(particle_properties_pore.get('particle_size_distribution_final', np.array([[10,1]]))[0,0]) # 简化取第一个的尺寸
    psd = particle_properties_pore.get('particle_size_distribution_final', np.array([[10.0,1.0]]))
    if psd is not None and len(psd) > 0 and len(psd[0]) > 0:
        mean_particle_size_um = psd[0][0] # 简化取第一个的尺寸
    else:
        mean_particle_size_um = 10.0 # 默认值


    # 极其简化的孔径估算：孔隙水力半径 R_h = (孔隙度 * 平均粒径) / ((1-孔隙度)*6) (对于球床)
    # mean_pore_size_um = (average_porosity * mean_particle_size_um) / ((1-average_porosity) * 6) if (1-average_porosity) > 1e-6 else mean_particle_size_um
    # 或者更简单地，孔径与粒径成正比，与堆积密度有关
    mean_pore_size_um = mean_particle_size_um * (average_porosity / (1-average_porosity))**0.5 if (1-average_porosity) > 1e-6 else mean_particle_size_um * 0.2


    return {
        'pore_structure_properties': {
            'porosity_profile': porosity_profile,
            'average_porosity': average_porosity,
            'pore_size_distribution': np.array([[mean_pore_size_um * 0.5, 0.3], [mean_pore_size_um, 0.4], [mean_pore_size_um * 1.5, 0.3]]), # [孔径_um, 比例]
            'mean_pore_size': mean_pore_size_um, # um
            'specific_surface_area_pore': 5.0 / (mean_particle_size_um * 1e-6) if mean_particle_size_um > 0 else 5e5, # m^2/m^3 (体积比表面积)
            'tortuosity': 1.5 / average_porosity if average_porosity > 0 else 5,
            'pore_connectivity_index': 0.75
        },
        'catalyst_distribution_in_pores': {'status': 'Catalyst analysis not implemented in stub'}
    }

def run_crack_analysis_stub(film_mechanical_input, stress_inducing_conditions, polymer_concentration):
    print("\n--- 调用 皲裂演化模型 (桩函数) ---")
    print(f"CRACK STUB: 接收到 film_mechanical_input 键: {list(film_mechanical_input.keys())}")
    print(f"CRACK STUB: 接收到 stress_inducing_conditions: {stress_inducing_conditions}")
    print(f"CRACK STUB: 接收到 polymer_concentration: {polymer_concentration}")

    thickness_m = film_mechanical_input.get('dry_film_thickness_actual', 30e-6) # 单位: m
    thickness_um = thickness_m * 1e6 # 转换为 µm

    cracking_data = get_cracking_degree_data_and_fits()

    # 浓度到数据集基础名称的映射
    concentration_to_dataset_map = {
        0.0: 'CB1_IC0',
        0.02: 'CB1_IC02',
        0.035: 'CB1_IC035',
        0.05: 'CB1_IC05',
        0.075: 'CB1_IC075',
        0.1: 'CB1_IC1',
    }
    
    dataset_base_name = None
    # 注意：直接用浮点数作为字典键进行查找可能因精度问题失败。
    # 更稳健的方式是迭代查找最接近的键或使用字符串键。
    if polymer_concentration in concentration_to_dataset_map:
        dataset_base_name = concentration_to_dataset_map[polymer_concentration]
    
    selected_fit_function_key = None
    if dataset_base_name:
        selected_fit_function_key = f'{dataset_base_name}_fit_function'
        print(f"信息: 根据聚合物浓度 {polymer_concentration}，动态选择了拟合函数: {selected_fit_function_key}")
    
    # 回退逻辑
    if not selected_fit_function_key or selected_fit_function_key not in cracking_data:
        if selected_fit_function_key: # 如果是因为key不在cracking_data中
            print(f"警告: 动态选择的拟合函数 {selected_fit_function_key} 在 cracking_data 中未找到。")
        else: # 如果是因为浓度没有匹配到
            print(f"警告: 未找到与聚合物浓度 {polymer_concentration} 精确匹配的数据集。")
        
        default_key = 'CB1_IC0_fit_function' # 默认使用0%聚合物的拟合
        # 你可以根据实际需求实现更智能的回退逻辑，例如选择最接近的浓度所对应的拟合
        # 例如:
        # closest_concentration = min(concentration_to_dataset_map.keys(), key=lambda k: abs(k - polymer_concentration))
        # default_base_name = concentration_to_dataset_map[closest_concentration]
        # default_key = f'{default_base_name}_fit_function'
        # print(f"将尝试回退到最接近的浓度 ({closest_concentration}) 的拟合函数: {default_key}")

        print(f"将回退到固定默认拟合函数: {default_key}")
        selected_fit_function_key = default_key

    cracking_degree_predicted_actual = 0.0 # 初始化
    if selected_fit_function_key in cracking_data:
        fit_function = cracking_data[selected_fit_function_key]
        
        original_thickness_range_key = selected_fit_function_key.replace('_fit_function', '_layer_thickness')
        if original_thickness_range_key in cracking_data:
            min_h = cracking_data[original_thickness_range_key].min()
            max_h = cracking_data[original_thickness_range_key].max()
            # 可以在此处添加警告，如果 thickness_um 超出 [min_h, max_h] 范围
            # print(f"信息: {selected_fit_function_key} 的原始数据厚度范围 [{min_h}um, {max_h}um]。当前厚度: {thickness_um:.2f}um")
        cracking_degree_predicted_actual = fit_function(thickness_um)
    else:
        # 这种情况理论上不应该发生，除非默认的 default_key 也不在 cracking_data 中 (例如 'CB1_IC0_fit_function' 被意外删除)
        print(f"严重警告: 最终选择的拟合函数 {selected_fit_function_key} 仍然无法在 cracking_data 中找到！开裂程度将为0。")

    # 确保开裂程度在合理范围 [0, 100]
    cracking_degree_predicted_actual = min(max(cracking_degree_predicted_actual, 0.0), 100.0)

    is_cracked_actual = cracking_degree_predicted_actual > 0.1 # 假设开裂程度大于0.1%即为开裂

    # 临界开裂厚度仍然使用桩函数的简化逻辑，因为模块不直接提供
    critical_thickness_predicted = 20e-6 # m 

    return {
        'crack_analysis_results': {
            'critical_cracking_thickness_predicted': critical_thickness_predicted,
            'is_cracked': is_cracked_actual,
            'cracking_degree_predicted': cracking_degree_predicted_actual, # %
            'crack_pattern_description': f'基于拟合 ({selected_fit_function_key}) 的模拟开裂' if is_cracked_actual else f'基于拟合 ({selected_fit_function_key}) 未预测开裂',
            'stress_state_film': {'max_principal_stress': np.max(stress_inducing_conditions.get('drying_stress_estimate',np.array([0]))) + stress_inducing_conditions.get('capillary_pressure_max',0)}
        }
    }


# --- 主控函数 ---
def run_integrated_simulation(global_initial_params, hss_params, drying_params, pore_params, crack_params, material_surface_char_global):
    """
    主函数，按顺序运行四个模型并传递数据。
    """
    print("--- 开始集成模拟 ---")
    all_results = {}

    # === 1. 高速剪切模型 ===
    hss_input_slurry_props = {
        'solid_concentration_vol': global_initial_params.get('solid_concentration_vol_initial'),
        'particle_size_distribution_initial': global_initial_params.get('particle_size_distribution_raw'),
        'polymer_concentration_mass': global_initial_params.get('polymer_concentration_initial'),
        'solvent_properties': global_initial_params.get('solvent_properties_base'),
        'particle_properties': global_initial_params.get('particle_properties_base'),
        'material_surface_char': material_surface_char_global 
    }
    hss_results = run_hss_simulation_stub(initial_slurry_properties=hss_input_slurry_props, shearing_process_params=hss_params)
    all_results['hss'] = hss_results
    print(f"高速剪切模型完成。输出摘要: {hss_results.get('processing_summary_hss', {}).get('status')}")


    # === 2. 干燥与偏析模型 ===
    sheared_slurry_state = hss_results.get('sheared_slurry_state', {})
    
    # 假设: phi1, phi2 初始浓度基于全局参数，但可以被HSS输出的更精确的固含量修正（如果HSS能提供）
    phi1_init_val = global_initial_params.get('phi1_after_mixing_or_hss', 0.2)
    phi2_init_val = global_initial_params.get('phi2_after_mixing_or_hss', 0.05)
    
    drying_input_initial_coating = {
        'wet_film_thickness': drying_params.get('wet_film_thickness_target'), 
        'component_concentrations_initial': { 
            'phi1_initial': phi1_init_val, 
            'phi2_initial': phi2_init_val, 
        },
        'solvent_fraction_initial': 1 - (phi1_init_val + phi2_init_val)
    }
    
    drying_input_material_props = {
        'D1': global_initial_params.get('D1_base'), 'D2': global_initial_params.get('D2_base'),
        'K1_sedimentation_coeff': global_initial_params.get('K1_base'), 
        'K2_sedimentation_coeff': global_initial_params.get('K2_base'),
        'Z_compressibility_factor': global_initial_params.get('Z_base'),
        'particle_density_1': global_initial_params.get('particle_properties_base',{}).get('density1', 2500),
        'particle_density_2': global_initial_params.get('particle_properties_base',{}).get('density2', 2200), # e.g. for polymer
        'solvent_density': global_initial_params.get('solvent_properties_base',{}).get('density', 1000),
        'solvent_evaporation_rate_coeff': drying_params.get('evaporation_coeff'),
        'material_surface_char_effective': sheared_slurry_state.get('material_surface_char_modified', material_surface_char_global) 
    }
    
    drying_results = run_drying_simulation_stub(
        initial_coating_params=drying_input_initial_coating,
        material_properties_drying=drying_input_material_props,
        drying_process_params=drying_params.get('process_conditions', {}), 
        numerical_params_drying=drying_params.get('numerical_settings', {})
    )
    all_results['drying'] = drying_results
    print(f"干燥模型完成。输出摘要: {drying_results.get('processing_summary_drying',{}).get('status')}")
    

    # === 3. 孔隙结构模型 ===
    dried_film_state = drying_results.get('dried_film_state', {})
    pore_model_input_film = {
        'solid_fraction_profile': dried_film_state.get('solid_fraction_profile_final'),
        'dry_film_thickness': dried_film_state.get('dry_film_thickness'),
        'phi1_profile': dried_film_state.get('phi1_profile_final'), # 传递以备用
        'phi2_profile': dried_film_state.get('phi2_profile_final'), # 传递以备用
    }
    pore_model_input_particles = {
        'particle_size_distribution_final': sheared_slurry_state.get('effective_particle_size_distribution'),
        'particle_shape_factor': global_initial_params.get('particle_shape_factor_base')
    }
    pore_results = analyze_pore_structure_stub(
        dried_film_input=pore_model_input_film,
        particle_properties_pore=pore_model_input_particles,
        catalyst_info=pore_params.get('catalyst_details', {})
    )
    all_results['pore'] = pore_results
    print(f"孔隙模型完成。平均孔隙率: {pore_results.get('pore_structure_properties',{}).get('average_porosity', 'N/A'):.3f}")


    # === 4. 皲裂演化模型 ===
    pore_structure_props = pore_results.get('pore_structure_properties', {})
    
    # 辅助函数 (此处内联定义，实际应放在外面或导入)
    def calculate_effective_property(base_value, porosity, exponent=2.0): # 简化版
        if base_value is None or porosity is None: return base_value
        return base_value * ((1 - porosity)**exponent) if porosity < 1 else base_value * 0.01

    film_young_modulus_eff = calculate_effective_property(
        global_initial_params.get('young_modulus_solid_film'),
        pore_structure_props.get('average_porosity'),
        exponent=global_initial_params.get('young_modulus_porosity_exponent', 2.5)
    )
    film_KIC_eff = calculate_effective_property(
        global_initial_params.get('KIC_solid_film'),
        pore_structure_props.get('average_porosity'),
        exponent=global_initial_params.get('KIC_porosity_exponent', 1.5)
    )

    crack_model_input_film_mech = {
        'dry_film_thickness_actual': dried_film_state.get('dry_film_thickness'),
        'substrate_properties': global_initial_params.get('substrate_properties_base'),
        'film_material_properties': { 
            'young_modulus_film': film_young_modulus_eff,
            'poisson_ratio_film': global_initial_params.get('poisson_ratio_solid_film'), 
            'fracture_toughness_KIC': film_KIC_eff,
            'thermal_expansion_coeff_film': global_initial_params.get('thermal_expansion_coeff_solid_film')
        },
        'porosity_average': pore_structure_props.get('average_porosity')
    }
    
    # 使用导入的 capillary_pressure 函数
    mean_pore_size_m = 0
    if pore_structure_props.get('mean_pore_size') is not None and pore_structure_props.get('mean_pore_size') > 0:
        mean_pore_size_m = pore_structure_props.get('mean_pore_size') * 1e-6 # 从 um 转换为 m
    
    capillary_pressure_val = 0
    if mean_pore_size_m > 0: # 只有当孔径有效时才计算
        capillary_pressure_val = capillary_pressure( 
            L=material_surface_char_global.get('surface_tension_solvent', 0.072), # N/m
            theta=material_surface_char_global.get('contact_angle_solvent_particle'), # degrees
            r_p=mean_pore_size_m # m
        )

    crack_model_input_stress = {
        'drying_stress_estimate': dried_film_state.get('residual_stress_profile_final'), 
        'capillary_pressure_max': capillary_pressure_val, 
        'temperature_change': crack_params.get('temperature_change_for_thermal_stress')
    }
    
    # 从全局参数获取聚合物浓度
    current_polymer_concentration = global_initial_params.get('polymer_concentration_initial', 0.05) # 默认为0.05

    crack_results = run_crack_analysis_stub(
        film_mechanical_input=crack_model_input_film_mech,
        stress_inducing_conditions=crack_model_input_stress,
        polymer_concentration=current_polymer_concentration # 新增参数传递
    )
    all_results['crack'] = crack_results
    print(f"皲裂模型完成。预测开裂程度: {crack_results.get('crack_analysis_results',{}).get('cracking_degree_predicted', 'N/A'):.2f}%")

    print("\n--- 所有模型运行完毕 (桩函数版本) ---")
    return all_results


# --- NumPy数组JSON序列化辅助类 ---
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# === 主执行块 ===
if __name__ == '__main__':
    # --- 1. 定义全局初始参数 ---
    GLOBAL_INITIAL_PARAMS = {
        'solid_concentration_vol_initial': 0.25, 
        'particle_size_distribution_raw': np.array([[10.0, 0.5], [20.0, 0.5]]), # [diameter_um, fraction]
        'polymer_concentration_initial': 0.05, # mass fraction in solvent
        'solvent_properties_base': {'viscosity': 0.001, 'density': 1000.0}, # Pa.s, kg/m^3
        'particle_properties_base': {'density1': 2500.0, 'density2': 1200.0}, # kg/m^3 (e.g. polymer)
        
        'phi1_after_mixing_or_hss': 0.20, # 初始组分1体积分数 (HSS后或混合后)
        'phi2_after_mixing_or_hss': 0.05, # 初始组分2体积分数
        
        'D1_base': 1e-13, 'D2_base': 5e-13, # m^2/s (扩散系数)
        'K1_base': 1e-9, 'K2_base': 1e-10, # m/s (沉降系数相关) - 需要根据phi10.py的K函数调整
        'Z_base': 1e6, # Pa (压缩模量相关) - 需要根据phi10.py的Z函数调整

        'particle_shape_factor_base': 1.0, #球形
        
        'substrate_properties_base': {
            'young_modulus_substrate': 70e9, #Pa
            'thermal_expansion_coeff_substrate': 23e-6 # 1/K
        },
        'young_modulus_solid_film': 5e9, # Pa (固相本体杨氏模量)
        'poisson_ratio_solid_film': 0.3,
        'KIC_solid_film': 0.5e6, # Pa.m^0.5 (固相本体断裂韧性)
        'thermal_expansion_coeff_solid_film': 50e-6, # 1/K
        'young_modulus_porosity_exponent': 2.5, # 用于计算有效杨氏模量
        'KIC_porosity_exponent': 1.5, # 用于计算有效断裂韧性
    }

    # --- 2. 定义材料表面特性参数 ---
    MATERIAL_SURFACE_CHAR_GLOBAL = {
        'zeta_potential': -30.0, # mV
        'contact_angle_solvent_particle': 60.0, # degrees (溶剂对颗粒)
        'specific_surface_area': 50.0, # m^2/g
        'surface_tension_solvent': 0.072 # N/m (例如水的表面张力)
    }

    # --- 3. 定义各个模型的工艺参数 ---
    HSS_PROCESS_PARAMS = {
        'shear_rate': 10000.0, # 1/s
        'shearing_time': 60.0, # s
        'temperature': 25.0 # degC
    }

    DRYING_PROCESS_PARAMS = {
        'wet_film_thickness_target': 50e-6, # m (目标湿膜厚度)
        'evaporation_coeff': 5e-6, # (m/s)/Pa or other units per core_solver needs
        'process_conditions': {'drying_time_total': 1800.0, 'temperature_drying': 60.0}, # s, degC
        'numerical_settings': {'nz': 50, 'nt': 1000} # 空间/时间离散点数
    }

    PORE_ANALYSIS_PARAMS = {
        'catalyst_details': {} # 暂时为空
    }

    CRACK_ANALYSIS_PARAMS = {
        # 从干燥温度(60C)冷却到室 ২০ C产生的温度变化
        # 'temperature_change_for_thermal_stress': (GLOBAL_INITIAL_PARAMS['thermal_expansion_coeff_solid_film'] - GLOBAL_INITIAL_PARAMS['substrate_properties_base']['thermal_expansion_coeff_substrate']) * (20.0 - DRYING_PROCESS_PARAMS['process_conditions']['temperature_drying'])
        # 上面的写法是直接计算应力的一部分，但通常模型输入的是温差本身
    }
    # 热失配应力计算通常是 E_film * (alpha_film - alpha_substrate) * (T_final - T_initial_stress_free)
    # 这里的 'temperature_change' 应该只是 delta_T, 在模型内部与其他参数组合计算应力
    # 假设 T_initial_stress_free 是干燥温度
    CRACK_ANALYSIS_PARAMS['temperature_change_for_thermal_stress'] = 20.0 - DRYING_PROCESS_PARAMS['process_conditions']['temperature_drying']


    # --- 4. 运行集成模拟 ---
    print("开始运行集成模拟程序 (integrated_coating_analysis.py)...")
    simulation_all_results = run_integrated_simulation(
        GLOBAL_INITIAL_PARAMS,
        HSS_PROCESS_PARAMS,
        DRYING_PROCESS_PARAMS,
        PORE_ANALYSIS_PARAMS,
        CRACK_ANALYSIS_PARAMS,
        MATERIAL_SURFACE_CHAR_GLOBAL
    )

    print("\n\n--- 集成模拟最终结果摘要 ---")
    # 为了避免NpEncoder问题，我们手动提取一些关键信息打印
    
    hss_output = simulation_all_results.get('hss', {}).get('sheared_slurry_state', {})
    drying_output = simulation_all_results.get('drying', {}).get('dried_film_state', {})
    pore_output = simulation_all_results.get('pore', {}).get('pore_structure_properties', {})
    crack_output = simulation_all_results.get('crack', {}).get('crack_analysis_results', {})

    print("  高速剪切模型:")
    if hss_output:
        print(f"    有效粒径分布 (示例): {hss_output.get('effective_particle_size_distribution')}")
        print(f"    最终浆料粘度 (示例): {hss_output.get('slurry_rheology',{}).get('viscosity_final')} Pa.s")
    else:
        print("    无输出。")

    print("  干燥与偏析模型:")
    if drying_output:
        print(f"    最终干膜厚度: {drying_output.get('dry_film_thickness') * 1e6:.2f} μm") # 转换为微米
        # print(f"    固相体积分数分布 (最终点): {drying_output.get('solid_fraction_profile_final')[-1] if drying_output.get('solid_fraction_profile_final') is not None else 'N/A'}")
    else:
        print("    无输出。")

    print("  孔隙结构模型:")
    if pore_output:
        print(f"    平均孔隙率: {pore_output.get('average_porosity'):.3f}")
        print(f"    平均孔径: {pore_output.get('mean_pore_size'):.3f} μm")
    else:
        print("    无输出。")

    print("  皲裂演化模型:")
    if crack_output:
        print(f"    是否预测开裂: {'是' if crack_output.get('is_cracked') else '否'}")
        print(f"    预测开裂程度: {crack_output.get('cracking_degree_predicted'):.2f} %")
        print(f"    预测临界开裂厚度: {crack_output.get('critical_cracking_thickness_predicted')*1e6:.2f} μm")
    else:
        print("    无输出。")

    # 如果需要完整的JSON输出，可以使用下面的代码（确保NpEncoder类已定义）
    # print("\n完整结果 (JSON格式):")
    # print(json.dumps(simulation_all_results, indent=4, cls=NpEncoder))
    
    print("\n集成模拟程序执行完毕。") 