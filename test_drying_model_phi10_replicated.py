import numpy as np
import matplotlib.pyplot as plt
from drying_model_phi10_replicated import solve_phi10_model_identically

# 设置 Matplotlib 以正确显示中文字符和负号
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为 SimHei
plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像时负号'-'显示为方块的问题

def run_and_plot_phi10_simulation():
    """
    运行与 phi_10.py 的 grap_phi_1 函数行为一致的模拟并绘制结果。
    """
    # --- 参数设置 (与 phi_10.py 中的 grap_phi_1 基本一致) ---
    L = 1.0          # 系统长度 (单位任意，保持一致即可)
    Nx = 100         # 空间离散点数
    t_final = 0.1    # 总模拟时间
    # dt_original = 0.0001 # phi_10.py 中的 dt，这里用 t_eval 控制输出点
    num_time_points = 200 # 输出时间点的数量 (可调整)
    t_eval = np.linspace(0, t_final, num_time_points) # 需要求解器输出解的时间点

    D1 = 0.001       # 溶质1的扩散系数
    D2 = 0.0005      # 溶质2的扩散系数 (phi_10.py 中 D2 = D1/2)
    
    Pe1_val = 1.0      # 溶质1的佩克莱数
    Pe2_val = 0.5      # 溶质2的佩克莱数 (phi_10.py 中 Pe2 = Pe1/2)
    
    E = 1.0          # 蒸发速率或等效对流因子 (phi_10.py 中 E_input)

    # 初始浓度 (phi_10.py 中 phi1_0 和 phi2_0 都是均匀的 phi_0 / 2)
    phi_0_total_initial = 0.1 # 初始总浓度
    phi1_initial = np.full(Nx, phi_0_total_initial / 2.0)
    phi2_initial = np.full(Nx, phi_0_total_initial / 2.0)
    
    # 反应参数 (phi_10.py 中 K 和 Z 的值)
    # K_val 和 Z_val 将基于 phi1 和 phi2 的全局平均值在 phi_10.py 中计算，但这里我们直接使用其示例值
    # 如果严格复制 grap_phi_1, K 和 Z 会在求解前确定，而不是动态变化。
    # 从 phi_10.py 的 grap_phi_1 调用 pde_system 的地方看，K 和 Z 是作为标量传入的。
    # K = 10, Z = 1 是 grap_phi_1 中定义的示例值。
    K_val = 10.0
    Z_val = 1.0

    phi_max = 0.9    # 模拟终止时的最大总溶质浓度 (phi_total_max)

    print(f"开始模拟: L={L}, Nx={Nx}, t_final={t_final}")
    print(f"参数: D1={D1}, D2={D2}, Pe1={Pe1_val}, Pe2={Pe2_val}, E={E}")
    print(f"K={K_val}, Z={Z_val}, phi_max={phi_max}")
    print(f"初始总浓度: {phi_0_total_initial}")

    # --- 运行模拟 ---
    actual_t, x_coords, phi1_sol, phi2_sol = solve_phi10_model_identically(
        phi1_initial=phi1_initial,
        phi2_initial=phi2_initial,
        L=L,
        Nx=Nx,
        t_span=(0, t_final),
        t_eval=t_eval, 
        phi_max=phi_max,
        D1=D1, D2=D2,
        Pe1_val=Pe1_val, Pe2_val=Pe2_val,
        E=E,
        K_val=K_val, Z_val=Z_val
    )

    print(f"模拟完成。实际结束时间: {actual_t[-1]:.4f} (共 {len(actual_t)} 个时间点)")

    # --- 结果处理与绘图 ---
    # 获取最后一个时间点的数据进行绘图
    # 如果由于事件触发，actual_t 可能比 t_eval 短
    final_time_index = -1 # 默认最后一个时间点
    final_tau = actual_t[final_time_index]

    phi1_final = phi1_sol[:, final_time_index]
    phi2_final = phi2_sol[:, final_time_index]
    phi_total_final = phi1_final + phi2_final

    # 应用 y_plot = phi / (1 - E * tau) 变换
    # 确保 (1 - E * tau) 不为零或负数以避免问题
    scaling_denominator = 1 - E * final_tau
    if scaling_denominator <= 1e-9: # 接近零或为负，可能表示模型参数/时间选择不当
        print(f"警告: 分母 (1 - E * tau) = {scaling_denominator:.3e} 非常小或非正。将不进行缩放绘图。")
        y_phi1_plot = phi1_final
        y_phi2_plot = phi2_final
        y_phi_total_plot = phi_total_final
        plot_title_suffix = f" (原始值, t={final_tau:.3f})"
    else:
        y_phi1_plot = phi1_final / scaling_denominator
        y_phi2_plot = phi2_final / scaling_denominator
        y_phi_total_plot = phi_total_final / scaling_denominator
        plot_title_suffix = f" (缩放值 y/(1-E$\tau$), t={final_tau:.3f}, E$\tau$={E*final_tau:.3f})"

    plt.figure(figsize=(12, 8))
    
    plt.plot(x_coords, y_phi1_plot, label=r'$\phi_1$')
    plt.plot(x_coords, y_phi2_plot, label=r'$\phi_2$')
    plt.plot(x_coords, y_phi_total_plot, label=r'$\phi_{total}$', linestyle='--', color='k')

    plt.xlabel('空间坐标 x/L')
    plt.ylabel(r'浓度参数 $\Phi / (1 - E \tau)$' if scaling_denominator > 1e-9 else '浓度 $\Phi$')
    plt.title('干燥模型溶质分布 (phi_10.py 复现)' + plot_title_suffix)
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0) # 浓度不应为负
    plt.show()

    # 还可以绘制一些其他 phi_10.py 中可能绘制的图，例如总溶质随时间的变化等
    # ...

if __name__ == '__main__':
    run_and_plot_phi10_simulation() 