%% 裂纹生成的临界应力（Griffith 模型）
function sigma_c = griffith_criteria(E, gamma, a)
    % 计算临界断裂应力
    sigma_c = sqrt((2 * E * gamma) / (pi * a));
end