import numpy as np
import matplotlib.pyplot as plt

def get_cracking_degree_data_and_fits():
    """
    定义不同聚合物含量的CB1悬浮液的开裂程度对层厚度的依赖数据，
    并进行多项式拟合。

    返回:
    dict: 包含原始数据和拟合函数/拟合值的字典。
          键的例子: 'CB1_IC0_layer_thickness', 'CB1_IC0_cracking_area',
                     'CB1_IC0_fit_coeffs', 'CB1_IC0_fitted_cracking_points'
    """
    datasets = {}

    # CB1_IC0 数据
    CB1_IC0_layer_thickness = np.array([5, 10.2, 15]) # 层厚度 (um)
    CB1_IC0_cracking_area = np.array([0, 1, 4.7])    # 龟裂覆盖率 (%)
    datasets['CB1_IC0'] = {'thickness': CB1_IC0_layer_thickness, 'area': CB1_IC0_cracking_area, 'degree': 2}

    # CB1_IC02 数据
    CB1_IC02_layer_thickness = np.array([5, 11.4, 19, 23])
    CB1_IC02_cracking_area = np.array([0, 1, 8.4, 9.7])
    datasets['CB1_IC02'] = {'thickness': CB1_IC02_layer_thickness, 'area': CB1_IC02_cracking_area, 'degree': 3}

    # CB1_IC035 数据
    CB1_IC035_layer_thickness = np.array([5, 10, 14, 17, 21])
    CB1_IC035_cracking_area = np.array([0, 0, 1, 4, 9])
    datasets['CB1_IC035'] = {'thickness': CB1_IC035_layer_thickness, 'area': CB1_IC035_cracking_area, 'degree': 3}

    # CB1_IC05 数据
    CB1_IC05_layer_thickness = np.array([5, 10, 15, 19.7, 21])
    CB1_IC05_cracking_area = np.array([0, 0, 0, 1, 1.8])
    datasets['CB1_IC05'] = {'thickness': CB1_IC05_layer_thickness, 'area': CB1_IC05_cracking_area, 'degree': 3}

    # CB1_IC075 数据
    CB1_IC075_layer_thickness = np.array([5, 10, 15, 19.3])
    CB1_IC075_cracking_area = np.array([0, 0, 0, 1])
    datasets['CB1_IC075'] = {'thickness': CB1_IC075_layer_thickness, 'area': CB1_IC075_cracking_area, 'degree': 1}

    # CB1_IC1 数据
    CB1_IC1_layer_thickness = np.array([5, 10, 13.7, 17])
    CB1_IC1_cracking_area = np.array([0, 0, 1, 5.4])
    datasets['CB1_IC1'] = {'thickness': CB1_IC1_layer_thickness, 'area': CB1_IC1_cracking_area, 'degree': 3}

    results = {}
    for name, data in datasets.items():
        results[f'{name}_layer_thickness'] = data['thickness']
        results[f'{name}_cracking_area'] = data['area']
        results[f'{name}_degree'] = data['degree']
        
        # 数据拟合
        # MATLAB的 'polyN' 中的 N 是阶数，numpy.polyfit 的第三个参数也是阶数
        coeffs = np.polyfit(data['thickness'], data['area'], data['degree'])
        fit_function = np.poly1d(coeffs)
        
        # 为了绘图，生成更平滑的拟合曲线点
        smooth_thickness = np.linspace(data['thickness'].min(), data['thickness'].max(), 100)
        fitted_cracking_points = fit_function(smooth_thickness)
        
        results[f'{name}_fit_coeffs'] = coeffs
        results[f'{name}_fit_function'] = fit_function # 可以直接调用此函数
        results[f'{name}_smooth_thickness_for_plot'] = smooth_thickness
        results[f'{name}_fitted_cracking_for_plot'] = fitted_cracking_points
        
    return results

# 如果直接运行此文件，可以添加一个简单的绘图示例
if __name__ == '__main__':
    # 设置中文显示
    plt.rcParams['font.sans-serif'] = ['SimHei'] # 指定默认字体
    plt.rcParams['axes.unicode_minus'] = False   # 解决保存图像是负号'-'显示为方块的问题

    plot_data = get_cracking_degree_data_and_fits()
    
    plt.figure(figsize=(10, 6))
    
    colors = ['b', 'g', 'r', 'c', 'm', 'k']
    line_styles = ['-', '--', '-.', ':', '-', '--']
    dataset_keys = ['CB1_IC0', 'CB1_IC02', 'CB1_IC035', 'CB1_IC05', 'CB1_IC075', 'CB1_IC1']
    legend_labels = []

    for i, key in enumerate(dataset_keys):
        # 绘制原始数据点 (可选)
        # plt.plot(plot_data[f'{key}_layer_thickness'], plot_data[f'{key}_cracking_area'], 
        #          'o', color=colors[i % len(colors)], markersize=5)
        
        # 绘制拟合曲线
        plt.plot(plot_data[f'{key}_smooth_thickness_for_plot'], 
                 plot_data[f'{key}_fitted_cracking_for_plot'], 
                 linestyle=line_styles[i % len(line_styles)], 
                 color=colors[i % len(colors)], 
                 linewidth=1.5)
        legend_labels.append(f'{key} 拟合曲线 (poly{plot_data[f"{key}_degree"]})')

    plt.xlabel('层厚度 (um)')
    plt.ylabel('龟裂覆盖率 (%)')
    plt.title('层厚与龟裂覆盖率关系 (MATLAB 脚本复现)')
    plt.legend(legend_labels, loc='upper left')
    plt.grid(True)
    plt.show() 