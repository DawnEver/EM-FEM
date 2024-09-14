from enum import Enum

import numpy as np

__all__ = ['Material_In_Use']
## CONSTANTs
MIU0 = 4 * np.pi * 1e-7


class Material(Enum):
    Steel = 'steel'
    Copper = 'copper'
    Magnet = 'magnet'
    Other = 'other'


## Materials
# TODO external Material Database
material_db_dict = {
    'steel': {
        'demo_steel': {
            'bh': None,  # if bh is None, seek miu_r
            'miu_r': 1000,
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


class Material_In_Use:
    def __init__(self, material_in_use_dict: dict) -> None:
        self.material_in_use_dict = {}
        for material_name, value_dict in material_in_use_dict.items():
            material_type = value_dict['material_type']
            maetrial_dict = material_db_dict[material_type][material_name]
            bh = maetrial_dict.get('bh')
            if bh is not None:
                ...
            else:
                miu = maetrial_dict['miu_r'] * MIU0

            self.material_in_use_dict[material_name] = {
                'material_type': material_type,
                'miu': miu,
            }

    def get_miu(self, material_name: str, B: float | None = None):
        return self.material_in_use_dict[material_name]['miu']
