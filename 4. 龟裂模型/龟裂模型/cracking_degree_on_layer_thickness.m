%% 不同聚合物含量的CB1悬浮液的开裂程度对层厚度的依赖
clear;clc;close all;
CB1_IC0_layer_thickness = [5, 10.2, 15]; % 层厚度 (um)
CB1_IC0_cracking_area = [0, 1, 4.7]; % 龟裂覆盖率 (%)

CB1_IC02_layer_thickness = [5, 11.4, 19 ,23];
CB1_IC02_cracking_area = [0, 1, 8.4, 9.7];

% CB1_IC035数据
CB1_IC035_layer_thickness = [5, 10, 14, 17, 21]; % 层厚度 (um)
CB1_IC035_cracking_area = [0, 0, 1, 4, 9]; % 龟裂覆盖率 (%)

% CB1_IC05数据
CB1_IC05_layer_thickness = [5, 10, 15, 19.7, 21]; % 层厚度 (um)
CB1_IC05_cracking_area = [0, 0, 0, 1, 1.8]; % 龟裂覆盖率 (%)

% CB1_IC075数据
CB1_IC075_layer_thickness = [5, 10, 15, 19.3]; % 层厚度 (um)
CB1_IC075_cracking_area = [0, 0, 0, 1]; % 龟裂覆盖率 (%)

% CB1_IC1数据
CB1_IC1_layer_thickness = [5, 10, 13.7, 17]; % 层厚度 (um)
CB1_IC1_cracking_area = [0, 0, 1, 5.4]; % 龟裂覆盖率 (%)


% 数据拟合 (指数模型)
CB1_IC0_fit_model = fit(CB1_IC0_layer_thickness', CB1_IC0_cracking_area', 'poly2');
CB1_IC02_fit_model = fit(CB1_IC02_layer_thickness', CB1_IC02_cracking_area', 'poly3');
CB1_IC035_fit_model = fit(CB1_IC035_layer_thickness', CB1_IC035_cracking_area', 'poly3');
CB1_IC05_fit_model = fit(CB1_IC05_layer_thickness', CB1_IC05_cracking_area', 'poly3');
CB1_IC075_fit_model = fit(CB1_IC075_layer_thickness', CB1_IC075_cracking_area', 'poly1');
CB1_IC1_fit_model = fit(CB1_IC1_layer_thickness', CB1_IC1_cracking_area', 'poly3');

% 计算拟合值

CB1_IC0_fitted_cracking = feval(CB1_IC0_fit_model, CB1_IC0_layer_thickness);
CB1_IC02_fitted_cracking = feval(CB1_IC02_fit_model, CB1_IC02_layer_thickness);
CB1_IC035_fitted_cracking = feval(CB1_IC035_fit_model, CB1_IC035_layer_thickness);
CB1_IC05_fitted_cracking = feval(CB1_IC05_fit_model, CB1_IC05_layer_thickness);
CB1_IC075_fitted_cracking = feval(CB1_IC075_fit_model, CB1_IC075_layer_thickness);
CB1_IC1_fitted_cracking = feval(CB1_IC1_fit_model, CB1_IC1_layer_thickness);

% 绘图
figure;
%plot(IC0_layer_thickness, CB1_IC0_cracking_area, 'ro', 'LineWidth', 1.5); hold on;
plot(CB1_IC0_layer_thickness, CB1_IC0_fitted_cracking, 'b-', 'LineWidth', 1.5);hold on;

%plot(CB1_IC02_layer_thickness, CB1_IC02_cracking_area, 'ro', 'LineWidth', 1.5); hold on;
plot(CB1_IC02_layer_thickness, CB1_IC02_fitted_cracking, 'b--', 'LineWidth', 1.5);hold on;

%plot(CB1_IC035_layer_thickness, CB1_IC035_cracking_area, 'ro', 'LineWidth', 1.5); hold on;
plot(CB1_IC035_layer_thickness, CB1_IC035_fitted_cracking, 'r-', 'LineWidth', 1.5);hold on;

%plot(CB1_IC05_layer_thickness, CB1_IC05_cracking_area, 'go', 'LineWidth', 1.5); hold on;
plot(CB1_IC05_layer_thickness, CB1_IC05_fitted_cracking, 'm-', 'LineWidth', 1.5);hold on;

%plot(CB1_IC075_layer_thickness, CB1_IC075_cracking_area, 'co', 'LineWidth', 1.5); hold on;
plot(CB1_IC075_layer_thickness, CB1_IC075_fitted_cracking, 'k-', 'LineWidth', 1.5);hold on;

%plot(CB1_IC1_layer_thickness, CB1_IC1_cracking_area, 'yo', 'LineWidth', 1.5); hold on;
plot(CB1_IC1_layer_thickness, CB1_IC1_fitted_cracking, 'g-', 'LineWidth', 1.5);
xlabel('层厚度 (um)');
ylabel('龟裂覆盖率 (%)');
legend('CB1_IC0 拟合曲线', 'CB1_IC02 拟合曲线', ...
       'CB1_IC035 拟合曲线', 'CB1_IC05 拟合曲线', ...
       'CB1_IC075 拟合曲线', 'CB1_IC1 拟合曲线', 'Location', 'NorthWest');
title('层厚与龟裂覆盖率关系');
grid on;

