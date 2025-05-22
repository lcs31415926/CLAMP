%% 裂纹生成的临界干燥应力
function sigma_crit = tirumkudulu_russel_model(r, h_layer, G, M, phi_rcp, L)
    % 计算临界干燥应力
    term1 = (2 * r / h_layer)^(2/3);
    term2 = (G * M * phi_rcp * r / (2 * L))^(1/3);
    term3 = (2 * L) / r; 
    sigma_crit = 0.1877 * term1 * term2 * term3;
end