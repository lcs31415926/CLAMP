%% 应力变化
function sigma_t = stress_function(d_t, delta_m_t, g, L, E_s, b, h_s, h_c_t, v_s)
    % 计算随时间变化的应力 sigma(t)
    
    % 常数 g（重力加速度）
    g = 9.81; % m/s^2
    
    % 计算分子部分
    numerator = (d_t - (3 * delta_m_t * g * L^3 / (2 * E_s * b * h_s^3))) * E_s * h_s^3;
    
    % 计算分母部分
    denominator = 3 * h_c_t * L^2 * (h_s + h_c_t) * (1 - v_s);
    
    % 计算应力
    sigma_t = numerator / denominator;
end

