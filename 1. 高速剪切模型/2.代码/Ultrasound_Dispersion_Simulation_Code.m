clear; clc; close all;

%% ==============================
% 输入参数定义与验证
% ==============================

fprintf('\n=============================\n');
fprintf('        输入参数定义\n');
fprintf('=============================\n');

% 用户输入超声波功率级别
original_power_levels = [20, 40, 60]; % 原始超声波功率级别，单位W

% 用户输入功率列表
input_powers = input('请输入超声波功率列表（W），例如 [20, 30, 40, 50, 60]: ');
if isempty(input_powers) || ~isnumeric(input_powers) || any(input_powers <= 0)
    error('请输入有效的正数功率列表。');
end

% 用户输入搅拌时间
stirring_time_min = input('请输入搅拌时间（分钟）: ');
if isempty(stirring_time_min) || ~isnumeric(stirring_time_min) || stirring_time_min <= 0
    error('请输入有效的正数搅拌时间。');
end
stirring_time_sec = stirring_time_min * 60; % 转换为秒

% 试验数据
% 超声波功率 (W), 时间 (min), 剪切速率 (s^-1), 团簇平均尺寸 (μm), 剪切粘度 (Pa·s), 温度变化 (℃)
data = [
    20, 0, 6, 6.5, 0.08, 22.1;
    20, 5, 6, 4.6, 0.073, 23;
    20, 10, 6, 3.5, 0.065, 24;
    20, 15, 6, 3.0, 0.062, 25;
    40, 0, 12, 6.5, 0.081, 22.5;
    40, 5, 12, 4, 0.071, 24;
    40, 10, 12, 3.2, 0.063, 26;
    40, 15, 12, 2.8, 0.060, 27.5;
    60, 0, 19, 6.5, 0.083, 23;
    60, 5, 19, 3.7, 0.070, 25;
    60, 10, 19, 2.8, 0.062, 27;
    60, 15, 19, 2.5, 0.058, 29;
];

homogenizer_diameter = 12e-3; % 超声波均质器直径，单位m

insertion_depth = 15e-3; % 插入深度，单位m
beaker_volume = 100e-3; % 烧杯体积，单位L
beaker_volume_m3 = beaker_volume / 1000; % 转换为m^3
suspension_depth = 36e-3; % 悬浮液深度，单位m

% 粒子与悬浮液参数
d0 = 0.24e-6; % 一次颗粒直径，单位m（2.1 μm）
desired_N0 = 1.20e16; % 初始颗粒数量密度，单位个/m^3
rho = 800; % 悬浮液密度，单位kg/m^3
phi_p_max = 0.63; % 最大填充率
phi = desired_N0 * (4/3)*pi*(d0/2)^3; % 固体体积分数
fprintf('调整后的固体体积分数 phi = %.2e\n', phi);

% 模型参数
alpha_b = 0.60; % 布朗凝聚系数
alpha_s = 0.58; % 剪切凝聚系数
alpha_i = 0.60; % 离聚物凝集系数
F0 = 2.10e-04; % 颗粒间键合能，单位J
k_b = 1.38e-23; % 玻尔兹曼常数，单位J/K
beta_m = 46; % 根据描述
beta_0 = 8.85e-12; % 真空介电常数，单位F/m
z = 1; % 电荷数
e_charge = 1.6e-19; % 元电荷，单位C
phi_s = 0.03; % 固体体积分数
N_A = 6.022e23; % 阿伏伽德罗常数

% 验证输入参数
if any(original_power_levels <= 0)
    warning('超声波功率级别必须为正值。');
end
if homogenizer_diameter <= 0
    warning('超声波均质器直径必须为正值。');
end
if insertion_depth <= 0 || insertion_depth > beaker_volume_m3^(1/3)
    warning('插入深度必须为正值且不超过烧杯的尺寸。');
end
if beaker_volume_m3 <= 0
    warning('烧杯体积必须为正值。');
end
if suspension_depth <= 0 || suspension_depth > beaker_volume_m3^(1/3)
    warning('悬浮液深度必须为正值且不超过烧杯的尺寸。');
end

fprintf('\n=============================\n');
fprintf('        模型参数\n');
fprintf('=============================\n');
fprintf('布朗凝聚系数 alpha_b = %.2f\n', alpha_b);
fprintf('剪切凝聚系数 alpha_s = %.2f\n', alpha_s);
fprintf('离聚物凝集系数 alpha_i = %.2f\n', alpha_i);
fprintf('结合能 F0 = %.2e J\n', F0);
fprintf('初始颗粒总数 N_0 = %.2e 个/m^3\n', desired_N0);
fprintf('粒子直径 d0 = %.2e m\n', d0);
fprintf('悬浮液密度 rho = %.2f kg/m^3\n', rho);

%% ==============================
% 计算初始质量
% ==============================
initial_mass = phi * rho;
fprintf('\n初始总质量 = %.2e kg/m^3\n', initial_mass);

%% ==============================
% 计算 Cu_values
% ==============================

% 提取剪切速率的稳态值（假设稳态在搅拌时间）
% 原始数据中稳态是最后一个时间点
shear_rates_original = zeros(length(original_power_levels),1);
for p = 1:length(original_power_levels)
    data_power_orig = data(data(:,1) == original_power_levels(p), :);
    shear_rates_original(p) = data_power_orig(end,3); % 稳态剪切速率
end

% 基于参考文献计算 Cu_values，公式：C_u = 4.7e-7 * shear_rate + 8.9e-6
Cu_values_original = 4.7e-7 * shear_rates_original + 8.9e-6;
fprintf('\n=============================\n');
fprintf('        计算 Cu_values\n');
fprintf('=============================\n');
for p = 1:length(original_power_levels)
    fprintf('功率 %d W, 剪切速率 %.2f s^{-1}, C_u = %.2e\n', original_power_levels(p), shear_rates_original(p), Cu_values_original(p));
end

%% ==============================
% 模拟计算过程
% ==============================

% 预定义不同功率对应的标记类型
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};

% 4. 将搅拌时间作为程序输入
t_total = stirring_time_sec; % 总模拟时间，单位秒

% 1. 使用插值方法允许任意功率输入
% 提取独立的实验功率和对应的剪切速率、Cu_values
unique_powers = original_power_levels;
shear_rates = shear_rates_original;
Cu_values = Cu_values_original;

% 创建插值函数
Cu_interp = @(P) interp1(unique_powers, Cu_values, P, 'linear', 'extrap');
shear_rate_interp = @(P) interp1(unique_powers, shear_rates, P, 'linear', 'extrap');

% 预分配颜色
% 所有绘图为黑色，使用不同的标记区分功率
colors = repmat('k', length(marker_styles),1); % 所有线条为黑色

% 初始化结果存储
results = struct();

% 系统体积
V = beaker_volume_m3; % 单位m^3

% 循环处理每个功率级别
for idx = 1:length(input_powers)
    power = input_powers(idx);
    fprintf('\n=============================\n');
    fprintf('  正在计算功率为 %d W 的情况\n', power);
    fprintf('=============================\n');
    
    % 使用插值方法获取数据
    unique_times_min = unique(data(:,2));
    interpolated_shear_rate = zeros(length(unique_times_min),1);
    interpolated_d_c_exp = zeros(length(unique_times_min),1);
    interpolated_viscosity = zeros(length(unique_times_min),1);
    interpolated_T_data = zeros(length(unique_times_min),1);
    
    for t = 1:length(unique_times_min)
        current_time = unique_times_min(t);
        data_at_time = data(data(:,2) == current_time, :);
        
        interpolated_shear_rate(t) = interp1(original_power_levels, data_at_time(:,3), power, 'linear', 'extrap');
        interpolated_d_c_exp(t) = interp1(original_power_levels, data_at_time(:,4), power, 'linear', 'extrap');
        interpolated_viscosity(t) = interp1(original_power_levels, data_at_time(:,5), power, 'linear', 'extrap');
        interpolated_T_data(t) = interp1(original_power_levels, data_at_time(:,6), power, 'linear', 'extrap');
    end
    
    % 构建插值后的 data_power
    data_power = [power * ones(length(unique_times_min),1), unique_times_min, interpolated_shear_rate, interpolated_d_c_exp, interpolated_viscosity, interpolated_T_data];
    
    times_min = data_power(:,2); % 时间，单位min
    times = times_min * 60; % 转换为秒
    shear_dot_data = data_power(:,3); % 剪切速率，单位s^-1
    T_data = data_power(:,6) + 273.15; % 转换为K
    d_c_exp = data_power(:,4) * 1e-6; % 转换为m
    
    fprintf('提取实验数据完成。\n');
    fprintf('时间范围：0 至 %.2f 秒\n', times(end));
    
    %% 定义时间步长和时间范围
    dt = 1e-1; % 时间步长，单位秒，调整为0.1秒以提高效率
    t_total_sim = times(end);
    t_steps = 0:dt:t_total_sim;
    num_steps = length(t_steps);
    
    fprintf('时间步长 dt = %.4f 秒，总步数 = %d\n', dt, num_steps);
    
    %% 获取稳态簇大小（15分钟）
    steady_state_time = 15 * 60; % 15分钟转换为秒
    steady_state_idx = find(t_steps >= steady_state_time, 1, 'first');
    if isempty(steady_state_idx)
        error('稳态时间超出模拟时间范围！');
    end
    d_c_steady = d_c_exp(end); % 实验稳态簇大小，单位m
    
    %% 计算 n_steady （簇内平均粒子数）
    n_steady = (d_c_steady / d0)^3; % 近似估计簇大小k
    
    fprintf('稳态簇大小 d_c_steady = %.2e m\n', d_c_steady);
    fprintf('稳态平均颗粒数 n_steady = %.2e\n', n_steady);
    
    %% 计算 epsilon_max 和 epsilon
    epsilon_max = 1 - phi / phi_p_max; % 公式 (5)
    epsilon = epsilon_max * (1 - n_steady^-0.4); % 公式 (4)
    
    fprintf('计算得到的 epsilon_max = %.2f\n', epsilon_max);
    fprintf('计算得到的 epsilon = %.2f\n', epsilon);
    
    %% 计算 phi_a
    phi_a = phi / (1 - epsilon); % 公式 (10)
    fprintf('计算得到的 phi_a = %.2e\n', phi_a);
    
    %% 计算 f
    f = (8 * phi_p_max)^(1/3); % 公式 (15)
    fprintf('计算得到的 f = %.2f\n', f);
    
    %% 计算 gamma 和 lambda_gamma
    if f * (1 - (phi_a^(1/3) / f)) <= 0
        gamma_val = 0; % 避免负数或零
        warning('计算gamma时遇到无效值，设为0。');
    else
        gamma_val = (phi_a / (f * (1 - (phi_a^(1/3) / f))))^(1/3); % 公式 (14)
    end
    fprintf('计算得到的 gamma = %.2f\n', gamma_val);
    
    numerator = 4 * (1 - gamma_val^7);
    denominator = 4 * (1 - gamma_val^10) - 25 * gamma_val^3 * (1 + gamma_val^4) + 42 * gamma_val^5;
    if denominator == 0
        lambda_gamma = 0; % 避免除零
        warning('计算lambda_gamma时遇到除零，设为0。');
    else
        lambda_gamma = numerator / denominator; % 公式 (13)
    end
    fprintf('计算得到的 lambda_gamma = %.2f\n', lambda_gamma);
    
    %% 计算相对粘度 eta_r 和悬浮液粘度 eta
    eta_r = 1 + 2.5 * lambda_gamma * phi; % 公式 (12)
    fprintf('计算得到的相对粘度 eta_r = %.2f\n', eta_r);
    
    % 假设温度在稳态为 T_steady
    T_steady = mean(T_data);
    eta0_steady = (2.414e-5) * 10^(247.8 / (T_steady - 140)); % Pa·s
    eta_steady = eta_r * eta0_steady; % 公式 (12)
    fprintf('计算得到的悬浮液粘度 eta = %.2e Pa·s\n', eta_steady);
    
    %% 计算 dc 和 N_b
    k_max = 30; % 最大k值，允许更大的簇
    % 初始化n_k，所有初始为k=1
    n_k = zeros(k_max, num_steps);
    n_k(1,1) = desired_N0;
    fprintf('初始条件设置完成：n_k(1,1) = %.2e\n', desired_N0);
    
    %% 计算初始的 k_avg
    total_clusters_initial = sum(n_k(:,1));
    if total_clusters_initial == 0
        k_avg = 1; % 避免除以零，设定为1
        warning('初始时刻总簇数为零，k_avg 被设定为1。');
    else
        k_avg = sum((1:k_max)' .* n_k(:,1)) / total_clusters_initial;
    end
    fprintf('初始平均簇大小 k_avg = %.2f\n', k_avg);
    
    %% 计算 dc 和 N_b
    dc = (k_avg / (1 - epsilon))^(1/3) * d0; % 公式 dc = [k / (1 - epsilon)]^(1/3) * d0
    Nb = (k_avg * d0) / (2 * dc); % 公式 (11)
    
    fprintf('计算得到的 dc = %.2e m\n', dc);
    fprintf('计算得到的 N_b = %.2e\n', Nb);
    
    %% 初始化凝聚核函数矩阵和破碎核函数向量
    b = zeros(k_max, k_max); % 凝聚核函数矩阵
    G_shear = zeros(k_max, 1); % 剪切破碎核函数向量
    G_ultrasound = zeros(k_max, 1); % 超声波破碎核函数向量
    
    % 预计算凝聚核函数
    for i = 1:k_max
        for j = 1:k_max
            b(i,j) = coagulation_kernel(i, j, alpha_b, alpha_i, alpha_s, eta_steady, T_steady, gamma_val, d0, k_b);
        end
    end
    
    % 计算剪切破碎核函数 G_shear(k)
    for i = 1:k_max
        G_shear(i) = breakup_kernel(i, power, shear_dot_data(end), F0, Cu_interp(power));
    end
    
    % 计算超声波破碎核函数 G_ultrasound(k)
    G_ultrasound = Cu_interp(power) * G_shear;
    
    %% 预计算破碎分布函数矩阵（F_shear 和 F_ultrasound）
    F_shear = zeros(k_max, k_max); % 剪切破碎分布函数矩阵
    F_ultrasound = zeros(k_max, k_max); % 超声波破碎分布函数矩阵
    
    for k_cluster = 1:k_max
        % 剪切破碎
        if G_shear(k_cluster) > 0
            if k_cluster == 1
                % k=1 不能破碎，保持不变
                F_shear(k_cluster, k_cluster) = 1;
            elseif k_cluster == 2
                % k=2 分裂成两个 k'=1 的粒子
                F_shear(k_cluster,1) = 2;
            elseif k_cluster == 3
                % k=3 分裂成三个 k'=1 的粒子
                F_shear(k_cluster,1) = 3;
            elseif k_cluster == 4
                % k=4 分裂成两个 k'=2 的粒子
                F_shear(k_cluster,2) = 2;
            elseif k_cluster == 5
                % k=5 分裂成一个 k'=1 和一个 k'=4 的粒子
                F_shear(k_cluster,1) = 1;
                F_shear(k_cluster,4) = 1;
            elseif k_cluster == 6
                % k=6 分裂成三个 k'=2 的粒子
                F_shear(k_cluster,2) = 3;
            elseif k_cluster == 7
                % k=7 分裂成一个 k'=1 和一个 k'=6 的粒子
                F_shear(k_cluster,1) = 1;
                F_shear(k_cluster,6) = 1;
            elseif k_cluster == 8
                % k=8 分裂成两个 k'=4 的粒子
                F_shear(k_cluster,4) = 2;
            elseif k_cluster == 9
                % k=9 分裂成三个 k'=3 的粒子
                F_shear(k_cluster,3) = 3;
            elseif k_cluster == 10
                % k=10 分裂成一个 k'=2 和一个 k'=8 的粒子
                F_shear(k_cluster,2) = 1;
                F_shear(k_cluster,8) = 1;
            else
                % 对于更大的 k，按需要定义破碎模式
                % 例如，全部破碎成 k'=1
                F_shear(k_cluster,1) = k_cluster;
            end
        else
            % 如果剪切破碎率为0，则保持簇不变
            F_shear(k_cluster, k_cluster) = 1;
        end
        
        % 超声波破碎
        if G_ultrasound(k_cluster) > 0
            if k_cluster == 1
                % k=1 不能破碎，保持不变
                F_ultrasound(k_cluster, k_cluster) = 1;
            elseif mod(k_cluster, 2) == 0 && k_cluster/2 <= k_max
                % 偶数k，分裂为两个 k'=k/2的簇
                F_ultrasound(k_cluster, k_cluster/2) = 2;
            elseif k_cluster > 1 && floor(k_cluster/2) >= 1 && ceil(k_cluster/2) <= k_max
                % 奇数k，分裂为一个 k'=floor(k/2) 和一个 k'=ceil(k/2) 的簇
                F_ultrasound(k_cluster, floor(k_cluster/2)) = 1;
                F_ultrasound(k_cluster, ceil(k_cluster/2)) = 1;
            else
                % 其他情况，全部分裂成k'=1
                F_ultrasound(k_cluster,1) = k_cluster;
            end
        else
            % 如果超声波破碎率为0，则保持簇不变
            F_ultrasound(k_cluster, k_cluster) = 1;
        end
    end

    %% 开始时间步进模拟
    fprintf('开始时间步进模拟...\n');
    tic; % 记录时间步进时间
    
    % 为了节省内存，只存储每store_interval步的结果
    store_interval = floor(num_steps / 1000); % 每千步存储一次
    if store_interval < 1
        store_interval = 1;
    end
    store_steps = ceil(num_steps / store_interval);
    
    % 初始化存储结构
    stored_time = [];
    stored_d_c = [];
    stored_epsilon = []; % 存储平均孔隙大小
    stored_mass_error = []; % 存储质量误差
    stored_eta = [];
    stored_eta_r = [];
    stored_T_sim = []; % 存储模拟温度
    
    for t_idx = 2:num_steps
        % 当前时间
        current_time = t_steps(t_idx);
        
        % 根据实验数据插值获取当前剪切速率和温度
        gamma_dot = interp1(times, shear_dot_data, current_time, 'linear', 'extrap');
        T = interp1(times, T_data, current_time, 'linear', 'extrap');
        
        % 更新介质粘度
        eta_0 = (2.414e-5) * 10^(247.8 / (T - 140)); % Pa·s
        eta_r = 1 + 2.5 * lambda_gamma * phi; % 公式 (12)
        eta = eta_r * eta_0; % 公式 (12)
        
        % 重新计算凝聚核函数和破碎核函数
        for i = 1:k_max
            for j = 1:k_max
                b(i,j) = coagulation_kernel(i, j, alpha_b, alpha_i, alpha_s, eta_0, T, gamma_val, d0, k_b);
            end
            G_shear(i) = breakup_kernel(i, power, gamma_dot, F0, Cu_interp(power));
        end
        
        % 更新超声波破碎核函数
        G_ultrasound = Cu_interp(power) * G_shear;
        
        % 计算平均簇大小 k
        total_clusters = sum(n_k(:,t_idx-1));
        if total_clusters == 0
            k_avg = 1; % 避免除零，设定为1
        else
            k_avg = sum((1:k_max)' .* n_k(:,t_idx-1)) / total_clusters;
        end
        
        % 计算 epsilon_max 和 epsilon
        epsilon_max = 1 - phi / phi_p_max; % 公式 (5)
        if k_avg <= 0
            epsilon = epsilon_max; % 避免负数
        else
            epsilon = epsilon_max * (1 - k_avg^-0.4); % 公式 (4)
        end
        
        % 计算 phi_a
        phi_a = phi / (1 - epsilon); % 公式 (10)
        
        % 计算 f
        f = (8 * phi_p_max)^(1/3); % 公式 (15)
        
        % 计算 gamma 和 lambda_gamma
        if f * (1 - (phi_a^(1/3) / f)) <= 0
            gamma_val = 0; % 避免负数或零
            warning('计算gamma时遇到无效值，设为0。');
        else
            gamma_val = (phi_a / (f * (1 - (phi_a^(1/3) / f))))^(1/3); % 公式 (14)
        end
        
        numerator = 4 * (1 - gamma_val^7);
        denominator = 4 * (1 - gamma_val^10) - 25 * gamma_val^3 * (1 + gamma_val^4) + 42 * gamma_val^5;
        if denominator == 0
            lambda_gamma = 0; % 避免除零
            warning('计算lambda_gamma时遇到除零，设为0。');
        else
            lambda_gamma = numerator / denominator; % 公式 (13)
        end
        
        % 更新相对粘度和悬浮液粘度
        eta_r = 1 + 2.5 * lambda_gamma * phi; % 公式 (12)
        eta = eta_r * eta_0; % 公式 (12)
        
        % 计算 dc 和 N_b
        dc = (k_avg / (1 - epsilon))^(1/3) * d0; % 公式 dc = [k / (1 - epsilon)]^(1/3) * d0
        Nb = (k_avg * d0) / (2 * dc); % 公式 (11)
        
        % 防止N_b为零或负数
        if Nb <= 0
            Nb = 1e-10;
            warning('计算N_b时遇到无效值，设为1e-10。');
        end
        
        % 初始化 dn_k_dt
        dn_k_dt = zeros(k_max,1);
        
        %% 凝聚项
        for k_cluster = 1:k_max
            % 凝聚增益：从i和j聚合形成k_cluster
            for i = 1:k_cluster-1
                j = k_cluster - i;
                if j > k_max
                    continue; % 超出范围
                end
                dn_k_dt(k_cluster) = dn_k_dt(k_cluster) + 0.5 * b(i,j) * n_k(i,t_idx-1) * n_k(j,t_idx-1);
            end
            % 凝聚损失：k_cluster与所有其他聚集形成更大的聚集
            dn_k_dt(k_cluster) = dn_k_dt(k_cluster) - (b(k_cluster, :) * n_k(:, t_idx-1)) * n_k(k_cluster,t_idx-1);
        end
        
        %% 剪切破碎项（剪切破坏）
        for k_cluster = 1:k_max
            if G_shear(k_cluster) > 0
                % 破碎损失
                dn_k_dt(k_cluster) = dn_k_dt(k_cluster) - G_shear(k_cluster) * n_k(k_cluster,t_idx-1);
                
                % 破碎增益
                for k_prime = 1:k_max
                    dn_k_dt(k_prime) = dn_k_dt(k_prime) + F_shear(k_cluster, k_prime) * G_shear(k_cluster) * n_k(k_cluster,t_idx-1);
                end
            end
        end
        
        %% 超声波破碎项（超声波破坏）
        for k_cluster = 1:k_max
            if G_ultrasound(k_cluster) > 0
                % 破碎损失
                dn_k_dt(k_cluster) = dn_k_dt(k_cluster) - G_ultrasound(k_cluster) * n_k(k_cluster,t_idx-1);
                
                % 破碎增益
                for k_prime = 1:k_max
                    dn_k_dt(k_prime) = dn_k_dt(k_prime) + F_ultrasound(k_cluster, k_prime) * G_ultrasound(k_cluster) * n_k(k_cluster,t_idx-1);
                end
            end
        end
        
        %% 更新 n_k
        n_k(:,t_idx) = n_k(:,t_idx-1) + dn_k_dt * dt;
        % 确保 n_k 非负
        n_k(:,t_idx) = max(n_k(:,t_idx), 0);
        
        %% 计算总质量
        mass_total = sum(n_k(:,t_idx) .* ((4/3) * pi * ((1:k_max)' * d0 / 2).^3 )) * rho;
        
        %% 计算质量误差
        mass_error = mass_total - initial_mass;
        
        %% 每 store_interval 步或最后一步存储一次结果
        if mod(t_idx, store_interval) == 0 || t_idx == num_steps
            stored_time(end+1) = t_steps(t_idx) / 60; % 转换为分钟
            
            % 计算当前时间步的平均粒径
            total_size = sum((1:k_max)' .* n_k(:,t_idx) .* d0);
            total_particles = sum(n_k(:,t_idx));
            if total_particles == 0
                d_c_t = 0;
            else
                d_c_t = (total_size / total_particles) * 1e6; % 转换为 μm
            end
            stored_d_c(end+1) = d_c_t;
            
            % 计算并存储 epsilon
            stored_epsilon(end+1) = epsilon;
            
            % 存储粘度和相对粘度
            stored_eta(end+1) = eta;
            stored_eta_r(end+1) = eta_r;
            
            % 存储模拟温度
            stored_T_sim(end+1) = T - 273.15; % 转换为 ℃
            
            stored_mass_error(end+1) = mass_error;
            
            % 打印进度
            fprintf('功率 %d W，时间步 %d/%d (时间: %.2f 分钟)\n', ...
                power, t_idx, num_steps, t_steps(t_idx)/60);
            
            % 打印特定 k 值的颗粒数量
            for kk = 1:5 % 打印 k=1 到 k=5
                fprintf('  颗粒数量 k=%d: %.2e\n', kk, n_k(kk, t_idx));
            end
            
            % 打印总质量和质量误差
            fprintf('  当前总质量：%.2e kg/m^3\n', mass_total);
            fprintf('  质量误差：%.2e kg/m^3\n', mass_error);
            fprintf('  平均簇大小 k_avg：%.2f\n', k_avg); % 新增打印 k_avg
            fprintf('  平均孔隙大小 epsilon：%.2f\n', epsilon); % 新增打印 epsilon
        end
    end
    
    toc; % 显示时间步进时间
    fprintf('时间步进模拟完成。\n');
    
    %% 存储结果
    results(idx).power = power;
    results(idx).time = stored_time; % 转换为分钟
    results(idx).d_c = stored_d_c;
    results(idx).epsilon = stored_epsilon; % 存储平均孔隙大小
    results(idx).eta = stored_eta;
    results(idx).eta_r = stored_eta_r;
    results(idx).T_sim = stored_T_sim; % 存储模拟温度
    results(idx).d_c_exp = d_c_exp * 1e6; % 转换为 μm
    results(idx).time_exp = unique_times_min;
    results(idx).gamma_dot = shear_dot_data;
    results(idx).mass_error = stored_mass_error;
    
    fprintf('功率 %d W 的模拟结果已存储。\n', power);
end
toc; % 显示总时间步进时间
fprintf('时间步进模拟完成。\n');

%% ==============================
% 2. 平均粒径应该是 dc，不是 epsilon，应该单独放一张图
% 以及 质量误差随时间变化和温度随时间变化的子图
% ==============================

figure('Name', '平均粒径 dc、质量误差和温度随时间变化', 'NumberTitle', 'off');

% 子图1: 平均粒径 dc 随时间变化
subplot(3,1,1);
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).d_c(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('平均粒径 d_c (μm)', 'FontSize', 14, 'FontWeight', 'bold');
title('平均粒径 d_c 随时间变化', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','best');
grid on;

% 子图2: 质量误差随时间变化
subplot(3,1,2);
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).mass_error(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('质量误差 Δm (kg/m^3)', 'FontSize', 14, 'FontWeight', 'bold');
title('质量误差 Δm 随时间变化', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','best');
grid on;

% 子图3: 温度随时间变化
subplot(3,1,3);
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).T_sim(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('温度 T (℃)', 'FontSize', 14, 'FontWeight', 'bold');
title('温度 T 随时间变化', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','best');
grid on;

fprintf('平均粒径、质量误差和温度随时间变化 图已生成。\n');

%% ==============================
% 5. 体现各种参量对粘度的影响，用 subplot 绘图！
% ==============================

figure('Name', '各种参量对粘度的影响', 'NumberTitle', 'off');

% 子图1: 粘度随时间变化
subplot(2,1,1);
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).eta(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 12);
ylabel('粘度 η (Pa·s)', 'FontSize', 12);
title('粘度随时间变化', 'FontSize', 14);
legend('Location','best');
grid on;

% 子图2: 相对粘度随时间变化
subplot(2,1,2);
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).eta_r(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 12);
ylabel('相对粘度 η_r', 'FontSize', 12);
title('相对粘度随时间变化', 'FontSize', 14);
legend('Location','best');
grid on;

fprintf('各种参量对粘度的影响 图已生成。\n');

% 子图3: 平均孔隙大小随时间变化
figure('Name', '平均孔隙大小随时间变化', 'NumberTitle', 'off');
hold on;
for idx = 1:length(input_powers)
    plot(results(idx).time(1:100:end), results(idx).epsilon(1:100:end), '-k', 'LineWidth', 2, ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('时间 (min)', 'FontSize', 12);
ylabel('平均孔隙 epsilon', 'FontSize', 12);
title('平均孔隙大小随时间变化', 'FontSize', 14);
legend('Location','best');
grid on;

fprintf('平均孔隙大小随时间变化 图 已生成。\n');

%% ==============================
% 新增绘图1: 超声波破碎比率 C_u vs 剪切速率
% ==============================

figure('Name', '超声波破碎比率 C_u vs 剪切速率', 'NumberTitle', 'off');
hold on;
for idx = 1:length(input_powers)
    % 仅绘制稳态时刻的剪切速率和 C_u
    steady_gamma_dot = results(idx).gamma_dot(end);
    steady_Cu = Cu_interp(results(idx).power);
    plot(steady_gamma_dot, steady_Cu, 'o', 'Color','k', 'MarkerFaceColor','k', ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('剪切速率 γ̇ (s^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('超声波破碎比率 C_u', 'FontSize', 14, 'FontWeight', 'bold');
title('超声波破碎比率 C_u vs 剪切速率', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','best');
grid on;
fprintf('超声波破碎比率 C_u vs 剪切速率 图已生成。\n');

%% ==============================
% 新增绘图2: 超声波破碎比率 C_u vs 超声波功率
% ==============================

figure('Name', '超声波破碎比率 C_u vs 超声波功率', 'NumberTitle', 'off');
hold on;
for idx = 1:length(input_powers)
    % 仅绘制一个点
    steady_Cu = Cu_interp(results(idx).power);
    plot(results(idx).power, steady_Cu, 'o', 'Color','k', 'MarkerFaceColor','k', ...
        'Marker', marker_styles{mod(idx-1, length(marker_styles))+1}, ...
        'DisplayName', sprintf('%d W', results(idx).power));
end
hold off;
xlabel('超声波功率 P (W)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('超声波破碎比率 C_u', 'FontSize', 14, 'FontWeight', 'bold');
title('超声波破碎比率 C_u vs 超声波功率', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location','best');
grid on;
fprintf('超声波破碎比率 C_u vs 超声波功率 图已生成。\n');

%% ==============================
% 最终结果输出
% ==============================

fprintf('\n=============================\n');
fprintf('        最终结果输出\n');
fprintf('=============================\n');
for idx = 1:length(input_powers)
    fprintf('功率 %d W:\n', results(idx).power);
    fprintf('  最终平均粒径：%.2f μm\n', results(idx).d_c(end));
    fprintf('  最终时间：%.2f min\n', results(idx).time(end));
    fprintf('  最终粘度：%.4e Pa·s\n', results(idx).eta(end));
    fprintf('  最终相对粘度：%.2f\n', results(idx).eta_r(end));
    fprintf('  最终温度：%.2f ℃\n', results(idx).T_sim(end));
end

%% ==============================
% 模拟完成
% ==============================
fprintf('=============================\n');
fprintf('         模拟完成\n');
fprintf('=============================\n');

%% ==============================
% 辅助函数
% ==============================

function b_ij = coagulation_kernel(i, j, alpha_b, alpha_i, alpha_s, eta_0, T, gamma_dot, d0, k_b)
    % 计算凝聚核函数 b_ij
    % 公式 (18) 和 (19)

    r_i = (i)^(1/3) * d0 / 2; % 凝聚体 i 的半径，单位m
    r_j = (j)^(1/3) * d0 / 2; % 凝聚体 j 的半径，单位m

    % 布朗凝聚
    b_brown = (2 * (alpha_b + alpha_i) * k_b * T) / (3 * eta_0) * (r_i + r_j)^2 / (r_i * r_j); % 公式 (18)

    % 剪切凝聚
    b_shear = (4 * alpha_s * gamma_dot) / 3 * (r_i + r_j)^3; % 公式 (19)

    % 总凝聚核函数
    b_ij = b_brown + b_shear;
end

function G_k = breakup_kernel(k, power, gamma_dot, F0, C_u)
    % 计算破碎核函数 G(k)
    % 根据公式 (12)

    if k == 1
        G_k = 0;
    else
        G_k = C_u * gamma_dot / (F0 * k); % G(k) = C_u * shear_rate / (F0 * k)
    end
end
