import numpy as np
from scipy.integrate import solve_ivp

# 内部核心事件逻辑：当 phi1_max + phi2_max >= phi_max_param 时触发
def _event_phi10_core_logic(t, y, phi_max_param, Nx_param_for_slice):
    """
    事件函数的核心逻辑，用于 solve_ivp。
    当两种溶质的最大浓度之和达到或超过 phi_max_param 时，模拟终止。

    参数:
    t (float): 当前时间 (未直接使用)。
    y (np.ndarray): 当前解向量 (phi1 和 phi2 的扁平化数组)。
    phi_max_param (float): 终止条件中的最大溶质总浓度。
    Nx_param_for_slice (int): 用于从 y 中正确切片出 phi1 和 phi2 的空间点数。

    返回:
    float: event_value，当 event_value <= 0 时事件触发。
    """
    phi1 = y[:Nx_param_for_slice]
    phi2 = y[Nx_param_for_slice : 2 * Nx_param_for_slice]
    return phi_max_param - (np.max(phi1) + np.max(phi2))

def pde_system_phi10_exact(t, y, D1, D2, Pe1_val, Pe2_val, K_val, Z_val, L, Nx, E):
    """
    定义了与 phi_10.py 中完全一致的偏微分方程组。
    描述了两种溶质 (phi1, phi2) 在一维空间中的对流、扩散和反应。

    dphi1/dt = D1 * d^2(phi1)/dx^2 - (Pe1_val*E) * d(phi1)/dx - K_val * phi1 * phi2 + Z_val * phi2
    dphi2/dt = D2 * d^2(phi2)/dx^2 - (Pe2_val*E) * d(phi2)/dx + K_val * phi1 * phi2 - Z_val * phi2

    边界条件: dphi1/dt = 0 和 dphi2/dt = 0 在 x=0 和 x=L 处。

    参数:
    t (float): 当前时间 (未直接使用，但 scipy.integrate.solve_ivp 需要)。
    y (np.ndarray): 当前解向量，包含 phi1 和 phi2 的扁平化值。
    D1, D2 (float): 溶质1和溶质2的扩散系数。
    Pe1_val, Pe2_val (float): 溶质1和溶质2的佩克莱数的基础值。
    K_val (float): 反应速率常数 (phi1 + phi2 -> 结合态)。
    Z_val (float): 反应速率常数 (结合态 -> phi1 + phi2)，此处简化为 phi2 -> phi1 + phi2 (phi_10.py 中 Z 是从 P 到 S2)。
                     在 phi_10.py 中 K 是 S1+S2 -> P, Z 是 P -> S2。
                     如果这里 K 是 phi1 * phi2, Z 是 phi2, 那么反应项是
                     -K*phi1*phi2 + Z*phi2  (for phi1)
                     +K*phi1*phi2 - Z*phi2  (for phi2, if Z represents dissociation of phi2 itself)
                     phi_10.py 的具体项是:
                     Reaction1 = - K_val * phi1 * phi2 + Z_val * phi2_P (phi2_P 代表产物浓度, 此处简化模型中没有显式产物)
                                 应为 K*phi1*phi2 (消耗) + Z*P (生成)
                                 phi_10.py: R1 = -K*S1*S2 + Z*P, R2 = -K*S1*S2 - Z*P (这里似乎有误,应该是 R2 = -K*S1*S2 + Z*P), RP = K*S1*S2 - Z*P
                                 参考 grap_phi_1 中的实际方程:
                                 dphi1_dt = ... - K_val * phi1 * phi2 + Z_val * phi2 (这里Z似乎直接作用于phi2)
                                 dphi2_dt = ... + K_val * phi1 * phi2 - Z_val * phi2 (这里Z似乎直接作用于phi2)
                                 这是 phi_10.py 中 pde_system 实现的方式。
    L (float): 系统长度。
    Nx (int): 空间离散点数。
    E (float): 蒸发速率或其他影响对流项的因子。

    返回:
    np.ndarray: phi1 和 phi2 的时间导数 (扁平化数组)。
    """
    phi1 = y[:Nx]
    phi2 = y[Nx : 2 * Nx]

    dx = L / (Nx - 1)

    # 使用 'edge' 填充模拟诺伊曼边界条件 (dphi/dx = 0) 用于计算导数
    phi1_padded = np.pad(phi1, 1, mode='edge')
    phi2_padded = np.pad(phi2, 1, mode='edge')

    # 计算一阶导数 (中心差分)
    grad_phi1_centered = (phi1_padded[2:] - phi1_padded[:-2]) / (2 * dx)
    grad_phi2_centered = (phi2_padded[2:] - phi2_padded[:-2]) / (2 * dx)

    # 计算二阶导数 (中心差分)
    # (phi_ip1 - 2*phi_i + phi_im1) / dx^2
    # 这等价于 ( (phi_padded[2:] - phi_padded[1:-1])/dx - (phi_padded[1:-1] - phi_padded[:-2])/dx ) / dx
    d2phi1_dx2 = (phi1_padded[2:] - 2 * phi1_padded[1:-1] + phi1_padded[:-2]) / (dx**2)
    d2phi2_dx2 = (phi2_padded[2:] - 2 * phi2_padded[1:-1] + phi2_padded[:-2]) / (dx**2)
    
    # 扩散项
    diff1 = D1 * d2phi1_dx2
    diff2 = D2 * d2phi2_dx2

    # 对流项 (v_eff = Pe_val * E)
    adv1 = -Pe1_val * E * grad_phi1_centered
    adv2 = -Pe2_val * E * grad_phi2_centered
    
    # 反应项 (与 phi_10.py 中的 pde_system 一致)
    reaction1 = -K_val * phi1 * phi2 + Z_val * phi2
    reaction2 = +K_val * phi1 * phi2 - Z_val * phi2

    dphi1_dt = diff1 + adv1 + reaction1
    dphi2_dt = diff2 + adv2 + reaction2

    # 应用边界条件: 两端浓度不随时间变化 (dphi/dt = 0)
    dphi1_dt[0] = 0
    dphi1_dt[-1] = 0
    dphi2_dt[0] = 0
    dphi2_dt[-1] = 0
    
    return np.concatenate((dphi1_dt, dphi2_dt))

def solve_phi10_model_identically(
    phi1_initial, phi2_initial, L, Nx, t_span, t_eval,
    phi_max, D1, D2, Pe1_val, Pe2_val, E, K_val, Z_val
):
    """
    求解与 phi_10.py 的 grap_phi_1 函数行为一致的干燥模型。

    参数:
    phi1_initial (np.ndarray): 溶质1的初始浓度分布。
    phi2_initial (np.ndarray): 溶质2的初始浓度分布。
    L (float): 系统长度。
    Nx (int): 空间离散点数。
    t_span (tuple): 时间积分的起始和结束时间 (t_start, t_end)。
    t_eval (np.ndarray): 需要输出解的特定时间点。
    phi_max (float): 溶质总浓度的上限，达到此值时模拟终止。
    D1, D2 (float): 溶质1和溶质2的扩散系数。
    Pe1_val, Pe2_val (float): 溶质1和溶质2的佩克莱数基础值。
    E (float): 蒸发速率或影响对流的因子。
    K_val (float): 反应速率常数 K。
    Z_val (float): 反应速率常数 Z。

    返回:
    tuple: (sol_t, x_coords, phi1_solution, phi2_solution)
        sol_t (np.ndarray): 求解器实际输出解的时间点。
        x_coords (np.ndarray): 空间坐标点。
        phi1_solution (np.ndarray): phi1 在 x_coords 和 sol_t 上的解，形状 (Nx, len(sol_t))。
        phi2_solution (np.ndarray): phi2 在 x_coords 和 sol_t 上的解，形状 (Nx, len(sol_t))。
    """
    if len(phi1_initial) != Nx or len(phi2_initial) != Nx:
        raise ValueError("初始条件数组的长度必须等于 Nx。")

    y0 = np.concatenate((phi1_initial, phi2_initial))
    x_coords = np.linspace(0, L, Nx)

    # solve_ivp 需要的参数元组 (传递给 pde_system_phi10_exact 和事件函数)
    pde_args = (D1, D2, Pe1_val, Pe2_val, K_val, Z_val, L, Nx, E)

    # 创建事件函数 lambda，它捕获 solve_phi10_model_identically 作用域中的 phi_max 和 Nx
    # 并将它们传递给 _event_phi10_core_logic。
    # *event_pde_args 对应于 pde_args。
    event_wrapper = lambda t_event, y_event, *event_pde_args: _event_phi10_core_logic(
        t_event, y_event, phi_max, event_pde_args[7] # event_pde_args[7] is Nx
    )
    event_wrapper.terminal = True  # 事件发生时终止积分
    event_wrapper.direction = 0    # 当事件函数从正到负或负到正穿过零时触发 (即 value = 0)

    # 使用 'BDF' 方法，与 phi_10.py 中的设置一致，适用于可能刚性的系统
    # dense_output=True 也是为了与 phi_10.py 保持一致，允许在 t_eval 之外进行插值
    sol = solve_ivp(
        pde_system_phi10_exact,
        t_span,
        y0,
        method='BDF',
        t_eval=t_eval,
        args=pde_args,
        events=event_wrapper,
        dense_output=True 
    )

    # 从解中提取 phi1 和 phi2
    # sol.y 的形状是 (2*Nx, n_time_points)
    phi1_solution = sol.y[:Nx, :]
    phi2_solution = sol.y[Nx:, :]
    
    return sol.t, x_coords, phi1_solution, phi2_solution

if __name__ == '__main__':
    # 此部分为模块的简单测试/演示代码
    print("drying_model_phi10_replicated.py 模块已加载。")
    print("它包含 solve_phi10_model_identically 函数用于模拟。")
    print("要运行具体模拟，请在测试脚本中导入并使用此函数。")

    # 可以在这里添加一个非常简单的示例运行，如果需要的话
    # 例如:
    # L_test = 1.0
    # Nx_test = 50
    # x_test = np.linspace(0, L_test, Nx_test)
    # phi1_0_test = 0.1 * np.exp(-((x_test - L_test/4)**2) / (2 * (L_test/10)**2))
    # phi2_0_test = 0.1 * np.exp(-((x_test - 3*L_test/4)**2) / (2 * (L_test/10)**2))
    # t_final_test = 0.1 # 短时间测试
    # t_points_test = np.linspace(0, t_final_test, 20)
    
    # D1_test, D2_test = 0.001, 0.0005
    # Pe1_val_test, Pe2_val_test = 1.0, 0.5
    # E_test = 1.0
    # K_val_test, Z_val_test = 10.0, 1.0
    # phi_max_test = 0.9
    
    # print(f"运行一个快速演示，t_final={t_final_test}...")
    # sol_t, x, phi1_sol, phi2_sol = solve_phi10_model_identically(
    #     phi1_initial=phi1_0_test,
    #     phi2_initial=phi2_0_test,
    #     L=L_test,
    #     Nx=Nx_test,
    #     t_span=(0, t_final_test),
    #     t_eval=t_points_test,
    #     phi_max=phi_max_test,
    #     D1=D1_test, D2=D2_test,
    #     Pe1_val=Pe1_val_test, Pe2_val=Pe2_val_test,
    #     E=E_test,
    #     K_val=K_val_test, Z_val=Z_val_test
    # )
    # print(f"模拟完成。 时间点数量: {len(sol_t)}")
    # print(f"phi1 在最后一个时间点的最大值: {np.max(phi1_sol[:, -1]):.4f}")
    # print(f"phi2 在最后一个时间点的最大值: {np.max(phi2_sol[:, -1]):.4f}")

    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(10, 6))
    # plt.plot(x, phi1_sol[:, -1], label=f'phi1 at t={sol_t[-1]:.2f}')
    # plt.plot(x, phi2_sol[:, -1], label=f'phi2 at t={sol_t[-1]:.2f}')
    # plt.xlabel('空间坐标 x')
    # plt.ylabel('浓度 phi')
    # plt.title('简单演示模拟结果')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    pass # 保持 if __name__ == '__main__': 块为空或仅包含打印语句以避免在导入时执行复杂操作 