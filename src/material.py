from enum import Enum

import numpy as np
from scipy.interpolate import interp1d

__all__ = ['Material_In_Use']
## CONSTANTs
MIU0 = 4 * np.pi * 1e-7


class MaterialType(Enum):
    Steel = 'steel'
    Copper = 'copper'
    Magnet = 'magnet'
    Other = 'other'


## Materials
# TODO external Material Database
material_db_dict = {
    'steel': {
        # 'demo_steel_miu': {
        'demo_steel': {
            'bh': None,  # if bh is None, seek miu_r
            'miu_r': 1000,
        },
        'demo_steel_bh': {
            'bh': np.array(
                [
                    [0, 0],
                    [0.2003, 238.7],
                    [0.3204, 318.3],
                    [0.40045, 358.1],
                    [0.50055, 437.7],
                    [0.5606, 477.5],
                    [0.7908, 636.6],
                    [0.931, 795.8],
                    [1.1014, 1114.1],
                    [1.2016, 1273.2],
                    [1.302, 1591.5],
                    [1.4028, 2228.2],
                    [1.524, 3183.1],
                    [1.626, 4774.6],
                    [1.698, 6366.2],
                    [1.73, 7957.7],
                    [1.87, 15915.5],
                    [1.99, 31831],
                    [2.04, 47746.5],
                    [2.07, 63662],
                    [2.095, 79577.5],
                    [2.2, 159155],
                    [2.4, 318310],
                ]
            )
        },
    },
    'copper': {
        'cu': {
            'miu_r': 1,
        },
    },
    'magnet': {},
    'other': {
        'air': {
            'miu_r': 1,
        },
    },
}


class Material:
    def __init__(self, material_type: str, miu_r: float, bh=None) -> None:
        self.material_type = material_type
        # B = miu * H
        # reluctivity = 1 / permeability(miu)
        if bh is not None:
            B_init = 1 / 1000 / MIU0
            interp_func = interp1d(bh[:, 0], bh[:, 1])
            self.reluctivity_func = lambda B: B_init if B == 0 else interp_func(B) / B
        else:
            reluctivity = 1 / miu_r / MIU0
            self.reluctivity_func = lambda _: reluctivity


class Material_In_Use:
    def __init__(self, material_in_use_dict: dict) -> None:
        self.material_in_use_dict = {}
        for material_name, value_dict in material_in_use_dict.items():
            material_type = value_dict['material_type']
            maetrial_dict = material_db_dict[material_type][material_name]
            bh = maetrial_dict.get('bh')
            self.material_in_use_dict[material_name] = Material(
                material_type=material_type, bh=bh, miu_r=maetrial_dict.get('miu_r')
            )

    def get_reluctivity(self, material_name: str, B: float | None = None):
        return self.material_in_use_dict[material_name].reluctivity_func(B)
