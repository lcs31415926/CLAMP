import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # 用于 3D 绘图

# 从新的 matlab_pore_model_py 包中导入各个模块
from matlab_pore_model_py import alignment_module
from matlab_pore_model_py import catalyst_composition_module
from matlab_pore_model_py import pore_distribution_module

def run_and_plot_alignment_model():
    """运行排列子模型并绘图"""
    print("运行排列子模型 (alignment_submodel)...\n")
    # 使用默认参数调用，这些参数与 MATLAB 脚本中的参数一致
    x, delta_DA, theta, t = alignment_module.calculate_alignment_and_contact_angle()

    plt.figure(figsize=(8, 10))
    plt.suptitle('离子omer排列与接触角 (复现 alignment_submodel.m)', fontsize=14)

    plt.subplot(2, 1, 1)
    plt.plot(x, delta_DA, 'b-')
    plt.xlabel('位置 (x)')
    plt.ylabel('分子排列 (ΔDA)') # MATLAB y轴是 \Delta DA
    plt.title('离子omer的分子排列')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(x, theta, 'r-')
    plt.xlabel('位置 (x)')
    plt.ylabel('接触角 (°)')
    plt.title('接触角')
    plt.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # 为总标题留出空间
    plt.show()
    print("排列子模型绘图完成。\n")

def run_and_plot_catalyst_composition_model():
    """运行催化剂组成子模型并绘图"""
    print("运行催化剂组成子模型 (catalyst_layer_composition_model)...\n")
    # 使用默认参数调用，这些参数与 MATLAB 脚本中的参数一致
    x_sol, a_ion_sol, alpha_free_sol, t_membrane = catalyst_composition_module.solve_catalyst_composition()
    
    # 检查膜厚度的数量级，如果需要，进行调整以匹配nm
    # MATLAB的V_20=0.39 (cm3/gC), A_2=6.57 (m2/gC)
    # V_20/A_2 -> cm3/m2 = 1e-6 m。 MATLAB 图显示 nm，可能需要 *1e9。
    # 或者 V_20(cm3/gC) * 1e3 / A_2(m2/gC) -> 0.39e3 / 6.57 = 59 nm 左右
    # 实际的 t_membrane 值需要看一下才知道转换因子
    # 假设 MATLAB 脚本中的 V_20 和 A_2 单位输入后，直接计算的结果配合绘图是 nm 级别
    # 如果 t_membrane 均值在 1e-7 到 1e-9 之间，则乘以 1e9 得到 nm
    # 如果均值在 1e-4 到 1e-6 之间，则乘以 1e6 得到 um
    # 经过初步计算，t_membrane 的值在几十的量级，因此，可能 MATLAB 脚本的单位处理隐含了转换
    # 或者 V_20 和 A_2 的输入值本身就调整过了。我们暂时直接使用计算值。
    # 如果绘图显示的数量级不对，可以取消下面这行注释并调整因子：
    # t_membrane_plot = t_membrane * 1 # 示例：乘以转换因子（如果需要）
    t_membrane_plot = t_membrane # 暂时不加转换，观察输出

    plt.figure(figsize=(8, 10))
    plt.suptitle('催化剂层组成与膜厚度 (复现 catalyst_layer_composition_model.m)', fontsize=14)

    plt.subplot(2, 1, 1)
    plt.plot(x_sol, a_ion_sol, 'b-', label='离子omer覆盖率 (a_ion)')
    plt.plot(x_sol, alpha_free_sol, 'r--', label='自由孔隙体积分数 (α_free)')
    plt.legend()
    plt.xlabel('位置 (x)')
    plt.ylabel('分数')
    plt.title('离子omer覆盖率和自由孔隙体积比例随 x 的变化')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(x_sol, t_membrane_plot, 'g-')
    plt.xlabel('位置 (x)')
    plt.ylabel('膜厚度 (nm)') # 假设单位是 nm，基于 MATLAB 图
    plt.title('膜厚度 t 随 x 的变化')
    plt.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    print("催化剂组成子模型绘图完成。\n")

def run_and_plot_pore_distribution_model():
    """运行孔径和接触角的统计分布子模型并绘图"""
    print("运行孔径和接触角分布模型 (pore_angle_distribution_model)...\n")
    # 使用默认参数调用，这些参数与 MATLAB 脚本中的参数一致
    dist_results = pore_distribution_module.calculate_pore_and_angle_distributions()

    pore_radius = dist_results['pore_radius_array']
    pore_volume1 = dist_results['pore_volume1']
    pore_volume2 = dist_results['pore_volume2']
    contact_angle = dist_results['contact_angle_array']
    cad_normalized = dist_results['cad_normalized']
    R_mesh = dist_results['R_meshgrid']
    C_mesh = dist_results['C_meshgrid']
    Z_joint = dist_results['Z_joint_distribution']
    theta_eff_deg = dist_results['theta_pore_effective_deg']

    # 绘制孔径分布
    plt.figure(figsize=(8, 6))
    plt.plot(pore_radius, pore_volume1, 'b-', label='一次孔隙体积分数')
    plt.plot(pore_radius, pore_volume2, 'r--', label='二次孔隙体积分数')
    plt.legend()
    plt.xlabel('孔径 (nm)')
    plt.ylabel('孔隙体积分数 (归一化 PSD * φ)')
    plt.title('孔径分布 (复现 pore_angle_distribution_model.m)')
    plt.grid(True)
    plt.show()

    # 绘制接触角分布
    plt.figure(figsize=(8, 6))
    plt.plot(contact_angle, cad_normalized, 'g-')
    plt.xlabel('接触角 (°)')
    plt.ylabel('概率密度')
    plt.title(f'接触角分布 (均值 ≈ {theta_eff_deg:.1f}°，复现)')
    plt.grid(True)
    plt.show()

    # 绘制联合分布
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(R_mesh, C_mesh, Z_joint, cmap='viridis', edgecolor='none')
    # MATLAB 使用 shading interp，matplotlib 中 plot_surface 默认会插值，或者用 `ax.set_zlim` 等调整
    # 对于更平滑的表面，可以增加 meshgrid 的点数，或使用更高级的插值
    fig.colorbar(surf, shrink=0.5, aspect=10, label='联合概率密度')
    ax.set_xlabel('孔径 (nm)')
    ax.set_ylabel('接触角 (°)')
    ax.set_zlabel('联合概率密度')
    ax.set_title('孔径和接触角的联合分布 (复现)')
    plt.show()
    print("孔径和接触角分布模型绘图完成。\n")

if __name__ == "__main__":
    # 设置 Matplotlib 中文显示
    plt.rcParams['font.sans-serif'] = ['SimHei'] # 指定默认字体
    plt.rcParams['axes.unicode_minus'] = False   # 解决保存图像是负号'-'显示为方块的问题

    print("开始执行 MATLAB 孔隙模型 Python 复现脚本...\n")
    
    run_and_plot_alignment_model()
    run_and_plot_catalyst_composition_model()
    run_and_plot_pore_distribution_model()
    
    print("所有孔隙模型计算和绘图完成。") 