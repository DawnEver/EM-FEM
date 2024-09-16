import numpy as np


def read_lipo_csv():
    trigInfoPath = '/Users/linxu/Files/Workspace/EM-FEM/data/lipo/trigInfo.csv'
    vertInfoPath = '/Users/linxu/Files/Workspace/EM-FEM/data/lipo/vertInfo.csv'

    with open(trigInfoPath, encoding='utf-8') as f:
        trigInfoMat = np.loadtxt(f, delimiter=',', dtype=np.int64)
    with open(vertInfoPath, encoding='utf-8') as f:
        vertInfoMat = np.loadtxt(f, delimiter=',')
    trigGroupMat = trigInfoMat[:, 3]
    trigInfoMat = trigInfoMat[:, :3] - 1  # matlab start from 1

    material_in_use_dict = {
        'cu': {'material_type': 'copper'},
        'demo_steel': {'material_type': 'steel'},
        'air': {'material_type': 'other'},
    }
    group_list = []
    max_group_id = trigGroupMat.max()
    # % Material 0: Air 1: Iron 2: Copper
    group_list = []
    group_cu = 2
    group_steel = 1
    group_air = 0
    for group_i in range(max_group_id + 1):
        if group_i == group_cu:
            material_type = 'copper'
            material_name = 'cu'
        elif group_i == group_steel:
            material_type = 'steel'
            material_name = 'demo_steel'
        elif group_i == group_air:
            material_type = 'other'
            material_name = 'air'
        else:
            material_type = 'other'
            material_name = 'air'

        trig_ids = np.where(trigGroupMat == group_i)  # index of target triangles
        group_list.append(
            {
                'group_name': group_i,
                'material_type': material_type,
                'material_name': material_name,
                'trig_ids': trig_ids,
            }
        )
    xs, ys = vertInfoMat[:, 0], vertInfoMat[:, 1]
    vertex_ids = np.concatenate(
        (np.where(xs == 0)[0], np.where(xs == xs.max())[0], np.where(ys == 0)[0], np.where(ys == ys.max())[0]), axis=0
    )
    boundary_dirichlet_0 = {0: vertex_ids}
    boundary_sym_odd = ('odd', None, None)
    boundary_sym_even = ('even', None, None)
    boundary_dict = {
        'dirichlet': boundary_dirichlet_0,
        'sym': (boundary_sym_odd, boundary_sym_even),
    }
    return trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict
