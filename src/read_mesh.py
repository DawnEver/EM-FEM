import numpy as np
import meshio

__all__ = ['read_gmsh', 'read_nastran', 'read_lipo_csv']


def read_gmsh(mesh_path: str):
    mesh = meshio.gmsh.read(mesh_path)
    vertInfoMat = mesh.points[:, 0:2]  # 2d fea
    trigInfoMat = mesh.cells_dict['triangle']
    n_trig = len(trigInfoMat)

    ## boundary condition & material in use & group
    boundary_dirichlet = {}
    boundary_sym_odd_1 = None
    boundary_sym_odd_2 = None
    boundary_sym_even_1 = None
    boundary_sym_even_2 = None

    group_list = []
    material_in_use_dict = {}
    n_tag_in_need = 2
    for key, value_dict in mesh.cell_sets_dict.items():
        tags = key.split('@')
        if len(tags) < n_tag_in_need:
            continue
        cell_type = tags[0]

        if cell_type == 'boundary':
            boundary_type = tags[1]
            boundary_id = tags[2]
            line_ids = value_dict['line']
            if boundary_type == 'dirichlet':
                lineInfoMat = mesh.cells_dict['line']
                vert_ids = np.unique(lineInfoMat[line_ids])
                boundary_dirichlet[float(boundary_id)] = vert_ids
            elif boundary_type == 'symmetry':
                if boundary_id == 'odd1':
                    boundary_sym_odd_1 = line_ids
                elif boundary_id == 'odd2':
                    boundary_sym_odd_2 = line_ids
                if boundary_id == 'even1':
                    boundary_sym_even_1 = line_ids
                elif boundary_id == 'even2':
                    boundary_sym_even_2 = line_ids
                else:
                    ...
            else:
                ...
        elif cell_type == 'group':
            group_name = tags[1]
            material_type = tags[2]
            material_name = tags[3]
            trig_ids = value_dict['triangle']
            group_list.append(
                {
                    'group_name': group_name,
                    'material_type': material_type,
                    'material_name': material_name,
                    'trig_ids': trig_ids,
                }
            )
            material_in_use_dict[material_name] = {
                'material_type': material_type,
            }
        else:
            ...

    boundary_sym_odd = ('odd', boundary_sym_odd_1, boundary_sym_odd_2)
    boundary_sym_even = ('even', boundary_sym_even_1, boundary_sym_even_2)

    boundary_dict = {
        'dirichlet': boundary_dirichlet,
        'sym': (boundary_sym_odd, boundary_sym_even),
    }

    # group id mat

    trigGroupMat = np.zeros((n_trig, 1), dtype=int)
    for i_group, group_dict in enumerate(group_list):
        trig_ids = group_dict['trig_ids']
        trigGroupMat[trig_ids] = i_group

    return trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict


def read_nastran(mesh_path: str):
    mesh = meshio.nastran.read(mesh_path)
    vertInfoMat = mesh.points[:, 0:2]  # 2d fea
    trigGroupMat = np.tile(mesh.cell_data_dict['nastran:ref']['triangle6'], 2)
    trigInfoMat = np.vstack((mesh.cells_dict['triangle6'][:, :3], mesh.cells_dict['triangle6'][:, 3:]))

    material_in_use_dict = {
        'cu': {'material_type': 'copper'},
        'demo_steel': {'material_type': 'steel'},
        'air': {'material_type': 'other'},
    }
    max_group_id = trigGroupMat.max()

    group_cu = [5, 9, 12, 14, 16, 18]
    group_steel = [1, 4, 8, 11, 13, 15, 17]
    group_air = [2, 3, 6, 7, 10]

    group_list = []
    for group_i in range(max_group_id + 1):
        if group_i in group_cu:
            material_type = 'copper'
            material_name = 'cu'
        elif group_i in group_steel:
            material_type = 'steel'
            material_name = 'demo_steel'
        elif group_i in group_air:
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
    min_error = 1e-6
    ro = 0.075
    ri = 0.0125
    xs, ys = vertInfoMat[:, 0], vertInfoMat[:, 1]
    vertex_ids = np.concatenate(
        (
            np.where(np.abs(xs**2 + ys**2 - ri**2) < min_error)[0],
            np.where(np.abs(xs**2 + ys**2 - ro**2) < min_error)[0],
        )
    )
    boundary_dirichlet_0 = {0: vertex_ids}
    vertex_ids_0 = np.where(xs < min_error)
    vertex_ids_1 = np.where(ys < min_error)
    boundary_sym_odd = ('odd', vertex_ids_0, vertex_ids_1)
    boundary_sym_even = ('even', None, None)
    boundary_dict = {
        'dirichlet': boundary_dirichlet_0,
        'sym': (boundary_sym_odd, boundary_sym_even),
    }
    return trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict


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
