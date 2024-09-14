import numpy as np
import matplotlib.pyplot as plt

from utils import sum_area
from logger import log
from solver import solve_magnetostatic
from mesh import read_nastran

## Read Mesh
mesh_path = """/Users/linxu/Files/Workspace/EM-FEM/data/
extraFine/Synchronous_Reluctance_Machine_MeshExport_ExtraFine_NoSlide.nas"""
mesh_path = """/Users/linxu/Files/Workspace/EM-FEM/data/
coarse/Synchronous_Reluctance_Machine_MeshExport_Coarse_NoSlide.nas"""
trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict = read_nastran(
    mesh_path=mesh_path
)

## current -> current density
n_turn = 100
group_current_dict = {  # group_id : I (current,A)
    5: 289.7777 * n_turn,
    9: 289.7777 * n_turn,
    12: -212.1320 * n_turn,
    14: -212.1320 * n_turn,
    16: -77.6457 * n_turn,
    18: -77.6457 * n_turn,
}
n_group = len(group_list)  # num of groups
group_current_density_list = np.zeros(n_group)  # group_id : J (current_density,A/m^2)

for group_id, current in group_current_dict.items():
    trig_ids = group_list[group_id]['trig_ids']  # group_dict['trig_ids']
    vertex_mat = vertInfoMat[trigInfoMat[trig_ids]]
    a_copper = sum_area(vertex_mat)  # sum the area of all triangles
    # a_copper = 1.8098e-4
    group_current_density_list[group_id] = current / a_copper

### Solve
S_mat, A_mat, T_mat, B_mat, Energy_mat = solve_magnetostatic(
    trigInfoMat=trigInfoMat,
    trigGroupMat=trigGroupMat,
    vertInfoMat=vertInfoMat,
    group_list=group_list,
    group_current_density_list=group_current_density_list,
    boundary_dict=boundary_dict,
    material_in_use_dict=material_in_use_dict,
)

### Post Process

total_Energy = sum(Energy_mat)
total_flux = max(A_mat)
msg = f'Total Energy: {total_Energy}; Total Flux: {total_flux}'
log(msg, 'INFO')

# Plot
figure = plt.figure()
plt.title('Flux Density')
plt.scatter(vertInfoMat[:, 0], vertInfoMat[:, 1], c=A_mat)
plt.colorbar()

plt.show()
