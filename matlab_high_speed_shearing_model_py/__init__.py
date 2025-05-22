# matlab_high_speed_shearing_model_py/__init__.py
# 从核心求解器模块导出主要功能

from .core_hss_solver import run_high_speed_shearing_model_from_document
 
__all__ = [
    'run_high_speed_shearing_model_from_document'
] 