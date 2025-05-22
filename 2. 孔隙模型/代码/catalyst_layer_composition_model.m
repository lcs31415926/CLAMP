%% 催化剂组成子模型
function catalyst_layer_composition_model
    k_A = 4.5; % 离子omer分散参数
    k_V = 0.6; % 离子omer体积减少因子
    V_20 = 0.39; % 二次孔隙体积 (cm3/gC)
    A_2 = 6.57; % 二次孔隙表面积 (m2/gC)
    x_max = 1.5; % x 的最大值

    % 初始条件
    y0 = [0; 1];

    % 求解微分方程
    [x, y] = ode45(@(x, y) ode_model(x, y, k_A, k_V), [0 x_max], y0);

    a_ion = y(:, 1);
    alpha_free = y(:, 2);

    % 膜厚度
    t = (V_20 .* (1 - alpha_free)) ./ (A_2 .* a_ion + eps); % 避免除以零

    figure;
    subplot(2, 1, 1);
    plot(x, a_ion, 'b-', x, alpha_free, 'r--');
    legend('a_ion', '\alpha_{free}');
    xlabel('x');
    ylabel('Fraction');
    title('离子omer覆盖率和自由孔隙体积比例随 x 的变化');

    subplot(2, 1, 2);
    plot(x, t, 'g-');
    xlabel('x');
    ylabel('Membrane Thickness (nm)');
    title('膜厚度 t 随 x 的变化');
end

% 微分方程函数
function dydx = ode_model(x, y, k_A, k_V)
    a_ion = y(1);
    alpha_free = y(2);

    da_ion_dx = (1 - a_ion) * k_A;
    d_alpha_free_dx = -alpha_free * k_V;

    dydx = [da_ion_dx; d_alpha_free_dx];
end