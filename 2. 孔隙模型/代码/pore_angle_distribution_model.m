%% 孔径和接触角的统计分布子模型
function pore_angle_distribution_model
    % 参数设置
    r_mu1 = 3; % 一次孔径的几何平均值 (nm)
    sigma_r1 = 0.8; % 一次孔径分布的标准差
    r_mu2 = 80; % 二次孔径的几何平均值 (nm)
    sigma_r2 = 1.2; % 二次孔径分布的标准差
    phi1 = 0.21; % 一次孔隙体积分数
    phi2 = 0.79; % 二次孔隙体积分数
    a_Pt = 0.4; % Pt 纳米颗粒的表面覆盖率
    a_ion = 0.6; % 离子omer的表面覆盖率
    theta_Pt = 10; % Pt 的接触角 (度)
    theta_ion = 90; % 离子omer的接触角 (度)
    r = 20; % 当前孔径 (nm)

    pore_radius = 0.1:0.1:200; % 孔径范围 (nm)

    % 一次孔和二次孔的孔径分布
    psd1 = lognpdf(pore_radius, log(r_mu1), sigma_r1);
    psd2 = lognpdf(pore_radius, log(r_mu2), sigma_r2);

    % 归一化孔径分布
    psd1 = psd1 / sum(psd1);
    psd2 = psd2 / sum(psd2);

    % 计算一次孔和二次孔的孔隙体积分数
    pore_volume1 = phi1 * psd1;
    pore_volume2 = phi2 * psd2;

    % 接触角范围
    contact_angle = 0:1:180; % 接触角范围 (度)

    % 接触角分布
    rad_Pt = r - 1.5; % Pt 纳米颗粒的半径
    rad_ion = r - 5; % 离子omer的半径

    if rad_Pt < 0
        rad_Pt = 0;
    end
    if rad_ion < 0
        rad_ion = 0;
    end

    % Pt 纳米颗粒和离子omer的覆盖率
    a_Pt = a_Pt * (rad_Pt / r);
    a_ion = a_ion * (rad_ion / r);

    % 有效接触角
    cos_theta_pore = a_Pt * cosd(theta_Pt) + a_ion * cosd(theta_ion);
    theta_pore = acosd(cos_theta_pore);

    % 接触角分布
    cad = normpdf(contact_angle, theta_pore, 5); % 假设标准差为 5
    cad = cad / sum(cad);

    % 绘制孔径分布
    figure;
    plot(pore_radius, pore_volume1, 'b-', pore_radius, pore_volume2, 'r--');
    legend('一次孔', '二次孔');
    xlabel('孔径 (nm)');
    ylabel('孔隙体积分数');
    title('孔径分布');

    % 绘制接触角分布
    figure;
    plot(contact_angle, cad, 'g-');
    xlabel('接触角 (^\circ)');
    ylabel('概率分布');
    title('接触角分布');

    % 绘制联合分布
    [R, C] = meshgrid(pore_radius, contact_angle);
    Z = zeros(length(contact_angle), length(pore_radius));

    for i = 1:length(pore_radius)
        for j = 1:length(contact_angle)
            Z(j, i) = pore_volume1(i) * cad(j) + pore_volume2(i) * cad(j);
        end
    end

    figure;
    surf(R, C, Z);
    shading interp;
    xlabel('孔径 (nm)');
    ylabel('接触角 (^\circ)');
    zlabel('联合概率分布');
    title('孔径和接触角的联合分布');
end