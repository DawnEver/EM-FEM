import numpy as np
import matplotlib.pyplot as plt
import meshio

from utils import sum_area
from logger import log

## CONSTANTs
MIU0 = 4 * np.pi * 1e-7

## Materials
# TODO external Material Database
material_db_dict = {
    'steel': {
        'iron': {
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

## Read Mesh
geo_path = '/Users/linxu/Files/Workspace/EM-FEM/data/test.msh'
mesh = meshio.gmsh.read(geo_path)

vertInfoMat = mesh.points[:, 0:2]  # 2d fea
trigInfoMat_raw = mesh.cells_dict['triangle']
n_vert = len(vertInfoMat)
n_trig = len(trigInfoMat_raw)

## boundary condition & material in use & group
boundary_dirichlet = {}
boundary_sym_odd_1 = None
boundary_sym_odd_2 = None
boundary_sym_even_1 = None
boundary_sym_even_2 = None

group_list = []

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

    else:
        ...

n_group = len(group_list)  # num of groups

boundary_sym_odd = (boundary_sym_odd_1, boundary_sym_odd_2)
boundary_sym_even = (boundary_sym_even_1, boundary_sym_even_2)

# add group id to trigInfoMat

group_id_mat = np.zeros((n_trig, 1), dtype=int)


for i_group, group_dict in enumerate(group_list):
    trig_ids = group_dict['trig_ids']
    group_id_mat[trig_ids] = i_group

trigInfoMat = np.hstack((group_id_mat, trigInfoMat_raw), dtype=int)


## current -> current density
group_current_dict = {
    0: 289.7777,  # group_id : I (current,A)
}
group_current_density_list = np.zeros(n_group)  # group_id : J (current_density,A/m^2)

for group_id, current in group_current_dict.items():
    trig_ids = group_list[group_id]['trig_ids']  # group_dict['trig_ids']
    vertex_mat = vertInfoMat[trigInfoMat_raw[trig_ids]]
    a_copper = sum_area(vertex_mat)  # sum the area of all triangles
    # a_copper = 1.8098e-4
    group_current_density_list[group_id] = current / a_copper


### S_mat * A_mat = T_mat: prepare S_mat & T_mat

S_mat = np.zeros((n_vert, n_vert))
A_mat = np.zeros(n_vert)
T_mat = np.zeros(n_vert)
last_B_mat = np.ones(n_trig)
B_mat = np.zeros(n_trig)
last_Energy_mat = np.ones(n_trig)
Energy_mat = np.zeros(n_trig)

rtol_B = 1e-2  # relative error tolerance of B
rtol_Energy = 1e-2
max_num_iter = 10
miu = 0

is_convergence = False

for i_iter in range(max_num_iter + 1):
    for i_trig in range(n_trig):
        # get x,y
        # [?]TODO calculate these coordinates once and save as matrix (Space exchange time)
        trig = trigInfoMat[i_trig]
        vertex_ids = trig[1:4]
        [[xi, yi], [xj, yj], [xk, yk]] = vertInfoMat[vertex_ids]

        ai = xj * yk - xk * yj
        bi = yj - yk
        ci = xk - xj
        aj = xk * yi - xi * yk
        bj = yk - yi
        cj = xi - xk
        ak = xi * yj - xj * yi
        bk = yi - yj
        ck = xj - xi
        delta = ((xj * yk - xk * yj) + xi * (yj - yk) + yi * (xk - xj)) / 2  # area

        # get miu of material
        group_id = trig[0]
        group_dict = group_list[group_id]
        material_type = group_dict['material_type']
        material_name = group_dict['material_name']
        material_dict = material_db_dict[material_type][material_name]
        if material_type == 'copper':
            # current
            current = group_current_density_list[group_id] * delta / 3
            T_mat[vertex_ids] += current

            miu = material_dict.get('miu_r') * MIU0
        elif material_type == 'steel':
            # calculate miu
            miu = material_dict.get('miu_r') * MIU0
            # last_B = last_B_mat[i_trig]
            # miu = f(B)
        elif material_type == 'magnet':
            miu = material_dict.get('miu_r') * MIU0
        else:
            miu = material_dict.get('miu_r') * MIU0

        # calculate B, Energy
        [Ai, Aj, Ak] = A_mat[vertex_ids]
        current_B = (
            np.sqrt((bi * Ai + bj * Aj + bk * Ak) ** 2 + (ci * Ai + cj * Aj + ck * Ak) ** 2) / 2 / delta
        )  # magnetic flux density in the
        B_mat[i_trig] = current_B
        Energy_mat[i_trig] = current_B**2 * delta * miu / 2

        Sii_raw = bi**2 + ci**2
        Sjj_raw = bj**2 + cj**2
        Skk_raw = bk**2 + ck**2
        Sij_raw = Sji_raw = bi * bj + ci * cj
        Sjk_raw = Skj_raw = bj * bk + cj * ck
        Sik_raw = Ski_raw = bi * bk + ci * ck

        S_mat[np.ix_(vertex_ids, vertex_ids)] = (
            np.array(
                [
                    [Sii_raw, Sij_raw, Sik_raw],
                    [Sji_raw, Sjj_raw, Sjk_raw],
                    [Ski_raw, Skj_raw, Skk_raw],
                ]
            )
            * miu
            / 4
            / delta
        )

    if is_convergence:
        # last loop to update B_mat & Energy_mat
        break

    # current error < relative error -> exit loop
    if np.allclose(Energy_mat, last_Energy_mat, rtol=rtol_Energy) and np.allclose(B_mat, last_B_mat, rtol=rtol_B):
        is_convergence = True

    last_Energy_mat = Energy_mat
    last_B_mat = B_mat

    # boundary condition
    # constant A
    for A_const, vertex_ids in boundary_dirichlet.items():
        n_vertex_ids = len(vertex_ids)
        if n_vertex_ids < 1:
            continue
        T_mat[vertex_ids] = 3 * A_const
        matrix = np.zeros((n_vertex_ids, n_vert))
        matrix[np.arange(n_vertex_ids), vertex_ids] = 1
        S_mat[vertex_ids] = matrix

    # symmetry boundary odd/even
    vertex_ids_1, vertex_ids_2 = boundary_sym_odd
    if vertex_ids_1 is not None and vertex_ids_2 is not None:
        n_vertex_ids_1 = len(vertex_ids_1)
        n_vertex_ids_2 = len(vertex_ids_2)
        if n_vertex_ids_1 > 1 and n_vertex_ids_1 != n_vertex_ids_2:
            T_mat[vertex_ids] = 0
            matrix = np.zeros((n_vertex_ids_1, n_vert))
            matrix[np.arange(n_vertex_ids_1), vertex_ids_1] = 1
            matrix[np.arange(n_vertex_ids_1), vertex_ids_2] = 1
            S_mat[vertex_ids_1] = matrix

    ### Solve
    A_mat = np.linalg.solve(S_mat, T_mat)  # magnetic vector potential

    log(f'Current iteration: {i_iter}', 'INFO')

### Post Process

total_Energy = sum(Energy_mat)
total_flux = max(A_mat)
msg = f'Total Energy: {total_Energy}; Total Flux: {total_flux}'
log(msg, 'INFO')


figure = plt.figure()
plt.title('Flux Density')
plt.scatter(vertInfoMat[:, 0], vertInfoMat[:, 1], c=A_mat)
plt.colorbar()
plt.show()
