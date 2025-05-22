%% 临界涂层厚度

function CCT = critical_coating_thickness(G, M, phi_rcp, L, p_cap_max)
    % 计算临界涂层厚度
    term1 = (G * M * phi_rcp / (2 * L))^(1/2);
    term2 = (2 * L / - p_cap_max)^(3/2);
    CCT = 0.64 * term1 * term2;
end