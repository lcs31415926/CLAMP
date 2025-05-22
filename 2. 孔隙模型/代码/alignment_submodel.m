%% 排列子模型
function alignment_submodel
    delta_DA0 = 0.25; % 初始排列值
    delta_DA_inf = 0.37; % 无限厚膜的排列值
    k_DA = 0.1; % 衰减常数
    t0 = 6.5; % 初始膜厚度 (nm)

    % 膜厚度变化
    x = 0:0.1:2; 
    t = t0 * exp(-k_DA * x);

    % 离子omer分子排列
    delta_DA = delta_DA_inf + (delta_DA0 - delta_DA_inf) .* exp(-(k_DA * (t - t0)));

    % 接触角参数
    theta_HI = 30;
    theta_HO = 120;

    % 计算接触角
    theta = acosd(delta_DA .* (cosd(theta_HI) - cosd(theta_HO)) + cosd(theta_HO));

    figure;
    subplot(2, 1, 1);
    plot(x, delta_DA, 'b-');
    xlabel('位置');
    ylabel('\Delta DA');
    title('离子omer的分子排列');

    subplot(2, 1, 2);
    plot(x, theta, 'r-');
    xlabel('位置');
    ylabel('接触角 (^\circ)');
    title('接触角');
end