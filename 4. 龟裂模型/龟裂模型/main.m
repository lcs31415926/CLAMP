%% 计算毛细压力
% 定义参数
gamma = 0.072; % 液体的表面张力，单位为 N/m
theta = 30; % 接触角，单位为度
r_p = 8e-6; % 孔隙半径，单位为 m

% 调用函数计算毛细压力
p_cap = capillary_pressure(gamma, theta, r_p);

% 输出结果
fprintf('毛细压力 p_cap = %.2f Pa\n', p_cap);


%% 计算临界断裂应力
% 定义参数
E = 4e9; % 杨氏模量，单位为 Pa
gamma = 1.0; % 表面能，单位为 N/m
a = 0.1e-3; % 裂纹长度，单位为 m

% 调用函数计算临界断裂应力
sigma_c = griffith_criteria(E, gamma, a);

% 输出结果
fprintf('临界断裂应力 sigma_c = %.2f MPa\n', sigma_c / 1e6);


%% 计算临界干燥应力
% 定义参数
r = 8e-6; % 颗粒半径，单位为 m
h_layer = 10e-6; % 涂层厚度，单位为 m
G = 0.6e9; % 剪切模量，单位为 Pa
M = 6; % 无量纲配位数
phi_rcp = 0.64; % 随机密堆积的无量纲体积分数
L = 0.072; % 液体的表面张力，单位为 N/m

% 调用函数计算临界干燥应力
sigma_crit = tirumkudulu_russel_model(r, h_layer, G, M, phi_rcp, L);

% 输出结果
fprintf('临界干燥应力 sigma_crit = %.2f Pa\n', sigma_crit);


%% 计算临界涂层厚度
% 定义参数
G = 0.6; % 剪切模量，单位为 GPa
M = 6; % 无量纲配位数
phi_rcp = 0.64; % 随机密堆积的无量纲体积分数
L = 0.072; % 液体的表面张力，单位为 N/m
p_cap_max = -304; % 最大毛细压力，单位为 Pa

% 调用函数计算临界涂层厚度
CCT = critical_coating_thickness(G, M, phi_rcp, L, p_cap_max);

% 输出结果
fprintf('临界涂层厚度 CCT = %.2f μm\n', CCT * 1e6);

%% 应力计算
% 定义参数
d_t = 0.01; % 时间相关的距离或长度，单位为 m
delta_m_t = 0.001; % 时间相关的质量变化，单位为 kg
g = 9.81; % 重力加速度，单位为 m/s^2
L = 0.1; % 特征长度，单位为 m
E_s = 70e9; % 杨氏模量，单位为 Pa
b = 0.05; % 宽度或线性尺寸，单位为 m
h_s = 0.001; % 某一层的厚度，单位为 m
h_c_t = 0.002; % 随时间变化的另一层厚度，单位为 m
v_s = 0.3; % 泊松比

% 调用函数计算应力
sigma_t = stress_function(d_t, delta_m_t, g, L, E_s, b, h_s, h_c_t, v_s);

% 输出结果
fprintf('随时间变化的应力 sigma(t) = %.2f Pa\n', sigma_t);