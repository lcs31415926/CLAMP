import numpy as np
import matplotlib.pyplot as plt
from .cracking_degree_module import get_cracking_degree_data_and_fits # 使用相对导入

def analyze_coating_cracks():
    """
    分析涂层开裂数据的主函数。
    该函数将:
    1. 从 cracking_degree_module 获取处理后的数据和拟合结果。
    2. (未来可以添加) 基于这些数据进行更深入的分析。
    3. (未来可以添加) 生成特定的图表或报告。
    """
    # 设置中文显示
    plt.rcParams['font.sans-serif'] = ['SimHei'] # 指定默认字体
    plt.rcParams['axes.unicode_minus'] = False   # 解决保存图像是负号'-'显示为方块的问题

    print("正在从 cracking_degree_module 获取数据...")
    crack_data_and_fits = get_cracking_degree_data_and_fits()
    print("数据获取完毕。")

    # 示例：打印出 CB1_IC0 的拟合系数
    print("\n示例数据处理:")
    if 'CB1_IC0_fit_coeffs' in crack_data_and_fits:
        print("CB1_IC0 拟合系数:", crack_data_and_fits['CB1_IC0_fit_coeffs'])
    if 'CB1_IC0_fit_function' in crack_data_and_fits:
        # 示例: 使用拟合函数预测特定厚度的开裂面积
        example_thickness = 12 # um
        predicted_area = crack_data_and_fits['CB1_IC0_fit_function'](example_thickness)
        print(f"CB1_IC0 在厚度为 {example_thickness}um 时，预测的龟裂覆盖率为: {predicted_area:.2f}%")

    # --- 绘图示例 ---
    # 这里可以根据需要绘制特定的图表
    # 例如，绘制所有数据集的原始数据点和拟合曲线，类似于 cracking_degree_module 中的示例
    
    plt.figure(figsize=(12, 7))
    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'orange', 'purple']
    line_styles = ['-', '--', '-.', ':', '-', '--', '-.', ':']
    # 这些键需要与 cracking_degree_module.py 中定义的 datasets 的键一致
    dataset_keys = ['CB1_IC0', 'CB1_IC02', 'CB1_IC035', 'CB1_IC05', 'CB1_IC075', 'CB1_IC1'] 
    
    print("\n开始生成图表...")
    for i, key in enumerate(dataset_keys):
        # 检查绘图所需的所有键是否存在于字典中
        thickness_key = f'{key}_layer_thickness'
        area_key = f'{key}_cracking_area'
        smooth_thickness_key = f'{key}_smooth_thickness_for_plot'
        fitted_area_key = f'{key}_fitted_cracking_for_plot'
        degree_key = f'{key}_degree'

        if all(k in crack_data_and_fits for k in [thickness_key, area_key, smooth_thickness_key, fitted_area_key, degree_key]):
            # 绘制原始数据点
            plt.plot(crack_data_and_fits[thickness_key], 
                     crack_data_and_fits[area_key], 
                     'o', 
                     color=colors[i % len(colors)], 
                     markersize=6,
                     label=f'{key} 原始数据')
            
            # 绘制拟合曲线
            plt.plot(crack_data_and_fits[smooth_thickness_key], 
                     crack_data_and_fits[fitted_area_key], 
                     linestyle=line_styles[i % len(line_styles)], 
                     color=colors[i % len(colors)], 
                     linewidth=2,
                     label=f'{key} 拟合 (poly{crack_data_and_fits[degree_key]})')
        else:
            print(f"警告: 数据集 {key} 的绘图所需数据不完整，已跳过。")
            # 打印出哪些键缺失，方便调试
            missing_keys = [k for k in [thickness_key, area_key, smooth_thickness_key, fitted_area_key, degree_key] if k not in crack_data_and_fits]
            print(f"数据集 {key} 缺失的键: {missing_keys}")


    plt.xlabel('层厚度 (μm)', fontsize=14)
    plt.ylabel('龟裂覆盖率 (%)', fontsize=14)
    plt.title('不同样品层厚度与龟裂覆盖率的关系分析', fontsize=16)
    
    # 仅当图例不为空时显示图例
    handles, labels = plt.gca().get_legend_handles_labels()
    if handles:
        plt.legend(loc='upper left', fontsize=10)
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout() # 自动调整子图参数，使之填充整个图像区域
    
    # 保存图像（可选）
    # output_figure_path = "coating_crack_analysis_plot.png"
    # plt.savefig(output_figure_path)
    # print(f"图表已保存至 {output_figure_path}")
    
    plt.show()
    print("图表显示完毕。")

    # --- 其他分析 ---
    # 在这里可以添加更多的分析代码
    # 例如，计算不同材料在特定厚度下的开裂程度差异，
    # 或者找出开裂程度对厚度最敏感的材料等。

    print("\ncoating_crack_analysis 执行完毕。")


if __name__ == '__main__':
    analyze_coating_cracks() 