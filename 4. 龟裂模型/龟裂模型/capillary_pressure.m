%% 毛细管压力
function p_cap = capillary_pressure(L, theta, r_p)
    % 将接触角从度转换为弧度
    theta_rad = deg2rad(theta);
    
    % 计算毛细压力
    p_cap = (2 * L * cos(theta_rad)) / r_p;
end