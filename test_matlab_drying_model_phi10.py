import unittest
import numpy as np
import matplotlib.pyplot as plt

# 从新的 matlab_drying_model_phi10_py 包导入所需的功能
from matlab_drying_model_phi10_py import (
    solve_pde_phi10_based, 
    PHI_MAX_PACKING, 
    DEFAULT_NOMINAL_TEMP_K
)

# 尝试设置matplotlib以支持中文显示 (通用方法)
try:
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
except Exception as e:
    print(f"[测试脚本matplotlib中文字体设置警告]: {e}")

def plot_phi_profiles_with_void(tau_points, phi1_sols, phi2_sols, phi_void_sols, xi_grid, E_param, Pe1_eff, Pe2_eff, N_val, title_suffix=""):
    """
    绘制phi1, phi2 和 phi_void 的剖面图，模拟phi_10.py中的绘图行为并加入空隙分布。
    x轴为 xi * (1-tau)，y轴为 phi_val / (1 - E*tau) (phi1, phi2) 或原始 phi_void。
    """
    num_time_points = len(tau_points)
    if num_time_points == 0:
        print(f"没有时间点可以绘图: {title_suffix}")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False) # 不共享Y轴，因为phi_void范围可能不同
    fig.suptitle(rf'$\phi_1, \phi_2, \phi_{{void}}$ 分布 (Pe1_eff={Pe1_eff:.2f}, Pe2_eff={Pe2_eff:.2f}, E={E_param}, N={N_val}){title_suffix}', fontsize=14)
    
    tau_for_plotting = tau_points 
    color_map = plt.cm.viridis(np.linspace(0, 1, num_time_points))

    for i in range(num_time_points):
        current_tau = tau_for_plotting[i]
        x_transformed = xi_grid * (1 - current_tau)

        # phi1 和 phi2 的 y 轴变换
        denominator_phi12 = (1 - E_param * current_tau)
        if abs(denominator_phi12) < 1e-9: 
            denominator_phi12 = 1e-9 * np.sign(denominator_phi12) if np.sign(denominator_phi12) != 0 else 1e-9
        
        y_phi1_transformed = phi1_sols[:, i] / denominator_phi12
        y_phi2_transformed = phi2_sols[:, i] / denominator_phi12
        y_phi_void = phi_void_sols[:, i] # 空隙不进行变换

        axes[0].plot(x_transformed, y_phi1_transformed, label=rf'$\tau={current_tau:.3f}$', color=color_map[i])
        axes[1].plot(x_transformed, y_phi2_transformed, label=rf'$\tau={current_tau:.3f}$', linestyle='--', color=color_map[i])
        axes[2].plot(x_transformed, y_phi_void, label=rf'$\tau={current_tau:.3f}$', linestyle=':', color=color_map[i])

    axes[0].set_xlabel(r'$z/H_0 = \xi \cdot (1-\tau)$')
    axes[0].set_ylabel(rf'$\phi_1 / (1 - E \cdot \tau)$ (E={E_param})') 
    axes[0].set_title(r'组分1 (小颗粒, $Pe_1$)')
    axes[0].legend(loc='best', fontsize='x-small')
    axes[0].grid(True, linestyle=':', alpha=0.7)
    axes[0].set_ylim(bottom=0)

    axes[1].set_xlabel(r'$z/H_0 = \xi \cdot (1-\tau)$')
    axes[1].set_ylabel(rf'$\phi_2 / (1 - E \cdot \tau)$ (E={E_param})') 
    axes[1].set_title(r'组分2 (大颗粒, $Pe_2$)')
    axes[1].legend(loc='best', fontsize='x-small')
    axes[1].grid(True, linestyle=':', alpha=0.7)
    axes[1].set_ylim(bottom=0)

    axes[2].set_xlabel(r'$z/H_0 = \xi \cdot (1-\tau)$')
    axes[2].set_ylabel(r'$\phi_{void}$') 
    axes[2].set_title(r'空隙体积分数')
    axes[2].legend(loc='best', fontsize='x-small')
    axes[2].grid(True, linestyle=':', alpha=0.7)
    axes[2].set_ylim(0, 1) # 空隙分数在0-1之间

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # 调整布局以适应总标题
    plt.show()

def plot_stratification_trend(parameter_values, stratification_percentages, parameter_name_latex, fixed_params_info_str):
    plt.figure(figsize=(8, 6))
    plt.plot(parameter_values, stratification_percentages, marker='o', linestyle='-')
    plt.xlabel(parameter_name_latex, fontsize=12)
    plt.ylabel(r"组分2顶部10%富集度 (%)", fontsize=12)
    plt.title(f"偏析趋势 vs {parameter_name_latex}\n({fixed_params_info_str})", fontsize=14)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.show()

class TestDryingModelPhi10Based(unittest.TestCase):

    def setUp(self):
        """测试前准备通用参数"""
        self.N_points = 6
        self.E_param = 1.0 # E parameter from phi_10.py, used in plotting adjustment
        self.phi1_initial = np.array([0.114, 0.115, 0.118, 0.120, 0.128, 0.138])
        self.phi2_initial = np.array([0.108, 0.109, 0.111, 0.120, 0.132, 0.152])
        self.tau_eval_points = [0.17, 0.34, 0.51, 0.683] # 这些是 *期望* 的评估点
        self.tau_min_sim = 0.17
        self.tau_max_sim = 0.99
        # 新增：默认温度和蒸发速率因子
        self.default_actual_temp_K = DEFAULT_NOMINAL_TEMP_K
        self.default_nominal_ref_temp_K = DEFAULT_NOMINAL_TEMP_K
        self.default_evap_factor = 1.0

    def run_simulation_and_assert(self, Pe1_nominal, Pe2_nominal, 
                                  actual_temp_K=None, nominal_ref_temp_K=None, evap_factor=None,
                                  test_name_suffix="", plot_profiles=True):
        
        # 使用setUp中的默认值（如果未提供特定值）
        act_temp = actual_temp_K if actual_temp_K is not None else self.default_actual_temp_K
        nom_temp = nominal_ref_temp_K if nominal_ref_temp_K is not None else self.default_nominal_ref_temp_K
        evap_f = evap_factor if evap_factor is not None else self.default_evap_factor

        results = solve_pde_phi10_based(
            phi1_0_array=self.phi1_initial,
            phi2_0_array=self.phi2_initial,
            Pe1_nominal=Pe1_nominal,
            Pe2_nominal=Pe2_nominal,
            E_parameter=self.E_param,
            N_points=self.N_points,
            tau_min=self.tau_min_sim,
            tau_max=self.tau_max_sim, 
            t_eval_points=self.tau_eval_points,
            actual_drying_temp_K=act_temp,
            nominal_ref_temp_K=nom_temp,
            E_evaporation_rate_factor=evap_f
        )
        sol_t, phi1_sols, phi2_sols, phi_void_sols, xi_grid, strat_phi2_perc, Pe1_eff, Pe2_eff, msg, status, E_out = results

        print(f"测试场景: {test_name_suffix} (Pe1_nom={Pe1_nominal}, Pe2_nom={Pe2_nominal}, T_act={act_temp}K, E_fac={evap_f})")
        print(f"  有效Pe: Pe1_eff={Pe1_eff:.3f}, Pe2_eff={Pe2_eff:.3f}")
        print(f"  求解器状态: {status}, 消息: {msg}")
        print(f"  求解时间点: {sol_t}")
        print(f"  组分2偏析量: {strat_phi2_perc:.2f}%")

        self.assertIn(status, [0, 1], f"求解器失败: {msg}")
        self.assertEqual(phi1_sols.shape, (self.N_points, len(sol_t)))
        self.assertEqual(phi2_sols.shape, (self.N_points, len(sol_t)))
        self.assertEqual(phi_void_sols.shape, (self.N_points, len(sol_t)))
        
        self.assertTrue(np.all(phi1_sols >= -1e-6), "phi1包含显著负值")
        if np.any(phi2_sols < -1e-6):
            print(f"警告: phi2 在场景 '{test_name_suffix}' 中出现负值: min(phi2)={np.min(phi2_sols)}")
        self.assertTrue(np.all(phi_void_sols >= -1e-6), "空隙体积分数包含显著负值")
        self.assertTrue(np.all(phi_void_sols <= 1.0 + 1e-6), "空隙体积分数超过1")
        self.assertTrue(np.isscalar(strat_phi2_perc) or np.isnan(strat_phi2_perc), "偏析量应为标量或NaN")

        if plot_profiles and len(sol_t) > 0:
            plot_phi_profiles_with_void(
                sol_t, phi1_sols, phi2_sols, phi_void_sols, xi_grid, E_out, 
                Pe1_eff, Pe2_eff, self.N_points, title_suffix=f" ({test_name_suffix})"
            )
        return results

    def test_figure4_conditions(self):
        """测试对应phi_10.py中fig 4的参数条件 (Pe1=0.18, Pe2=0.37)"""
        self.run_simulation_and_assert(Pe1_nominal=0.18, Pe2_nominal=0.37, test_name_suffix="Fig4")

    def test_figure5_conditions(self):
        """测试对应phi_10.py中fig 5的参数条件 (Pe1=0.7, Pe2=0.14)"""
        self.run_simulation_and_assert(Pe1_nominal=0.7, Pe2_nominal=0.14, test_name_suffix="Fig5")

    def test_figure6_conditions(self):
        """测试对应phi_10.py中fig 6的参数条件 (Pe1=2.8, Pe2=5.6)"""
        self.run_simulation_and_assert(Pe1_nominal=2.8, Pe2_nominal=5.6, test_name_suffix="Fig6")

    def test_identical_pe_numbers(self):
        """测试Pe数相同的情况"""
        self.run_simulation_and_assert(Pe1_nominal=0.5, Pe2_nominal=0.5, test_name_suffix="Pe1=Pe2=0.5")

    def test_large_pe_difference(self):
        """测试Pe数差异较大的情况 (Pe1 >> Pe2)"""
        self.run_simulation_and_assert(Pe1_nominal=5.0, Pe2_nominal=0.1, test_name_suffix="Pe1>>Pe2")
    
    def test_large_pe_difference_swapped(self):
        """测试Pe数差异较大的情况 (Pe2 >> Pe1)"""
        self.run_simulation_and_assert(Pe1_nominal=0.1, Pe2_nominal=5.0, test_name_suffix="Pe2>>Pe1")

    def test_temperature_effect_on_stratification(self):
        Pe1_nom, Pe2_nom = 0.7, 1.4 # 最大偏析条件的名义Pe数
        temps_to_test_K = [DEFAULT_NOMINAL_TEMP_K - 20, DEFAULT_NOMINAL_TEMP_K, DEFAULT_NOMINAL_TEMP_K + 20]
        strat_results = []
        print("\n--- 测试温度对偏析的影响 --- (基准Pe1=0.7, Pe2=1.4)")
        for temp_K in temps_to_test_K:
            results = self.run_simulation_and_assert(
                Pe1_nom, Pe2_nom, 
                actual_temp_K=temp_K, 
                test_name_suffix=f"T={temp_K}K",
                plot_profiles=False # 通常趋势图中不绘制每个剖面
            )
            strat_results.append(results[5]) # results[5] is stratification_phi2_percent
        
        # 预期：温度升高 -> Pe_eff 降低 -> 偏析减弱 (如果Pe>1，影响可能复杂)
        # 这里只是简单检查是否有值，具体趋势的断言比较复杂，先绘图观察
        self.assertEqual(len(strat_results), len(temps_to_test_K))
        print(f"温度对偏析影响结果 (T_K vs Strat%): {list(zip(temps_to_test_K, strat_results))}")
        plot_stratification_trend(temps_to_test_K, strat_results, r"实际干燥温度 (K)", f"Pe1_nom={Pe1_nom}, Pe2_nom={Pe2_nom}")

    def test_evaporation_effect_on_stratification(self):
        Pe1_nom, Pe2_nom = 0.7, 1.4 # 最大偏析条件的名义Pe数
        evap_factors_to_test = [0.5, 1.0, 1.5, 2.0]
        strat_results = []
        print("\n--- 测试蒸发速率因子对偏析的影响 --- (基准Pe1=0.7, Pe2=1.4)")
        for factor in evap_factors_to_test:
            results = self.run_simulation_and_assert(
                Pe1_nom, Pe2_nom, 
                evap_factor=factor, 
                test_name_suffix=f"E_factor={factor}",
                plot_profiles=False
            )
            strat_results.append(results[5])
        
        # 预期：蒸发因子增大 -> Pe_eff 增大 -> 偏析增强 (如果Pe>1，影响可能复杂)
        self.assertEqual(len(strat_results), len(evap_factors_to_test))
        print(f"蒸发因子对偏析影响结果 (E_factor vs Strat%): {list(zip(evap_factors_to_test, strat_results))}")
        plot_stratification_trend(evap_factors_to_test, strat_results, r"蒸发速率因子 ($E_{factor}$)", f"Pe1_nom={Pe1_nom}, Pe2_nom={Pe2_nom}")

    def test_custom_initial_conditions_example_with_new_features(self):
        """测试自定义初始条件，并验证新功能（如图表和偏析计算）。"""
        phi1_top_rich = np.linspace(0.2, 0.1, self.N_points)
        phi2_bottom_rich = np.linspace(0.1, 0.2, self.N_points)
        
        original_phi1 = self.phi1_initial
        original_phi2 = self.phi2_initial

        self.phi1_initial = phi1_top_rich
        self.phi2_initial = phi2_bottom_rich
        print(f"\n运行自定义初始条件测试 (含新特性): phi1_0={self.phi1_initial}, phi2_0={self.phi2_initial}")
        self.run_simulation_and_assert(
            Pe1_nominal=0.5, Pe2_nominal=0.2, 
            actual_temp_K=DEFAULT_NOMINAL_TEMP_K + 10, # 略微升温
            evap_factor=1.2, # 略微增强蒸发
            test_name_suffix="CustomInit_NewFeatures"
        )
        
        self.phi1_initial = original_phi1
        self.phi2_initial = original_phi2

if __name__ == '__main__':
    unittest.main() 