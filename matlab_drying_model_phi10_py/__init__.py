# 初始化 matlab_drying_model_phi10_py 包
# 从子模块中导出主要功能，方便外部调用

from .core_solver import (
    solve_pde_phi10_based,
    PHI_MAX_PACKING,
    DEFAULT_NOMINAL_TEMP_K,
    K_original, 
    Z_original,
    pde_system_phi10_formulation,
    event_reached_phi_max_original,
    calculate_stratification
)

__all__ = [
    'solve_pde_phi10_based',
    'PHI_MAX_PACKING',
    'DEFAULT_NOMINAL_TEMP_K',
    'K_original',
    'Z_original',
    'pde_system_phi10_formulation',
    'event_reached_phi_max_original',
    'calculate_stratification'
] 