import numpy as np

from utils import sum_area, calc_centroid
from logger import log
from solver import solve_magnetostatic
from mesh import read_gmsh
from plot import plot_map, PlotMapType

## Read Mesh
mesh_path = '/Users/linxu/Files/Workspace/EM-FEM/data/test-DI_Chong.msh'
mesh_path = '/Users/linxu/Files/Workspace/EM-FEM/data/test-T_T_Lipo.msh'
trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict = read_gmsh(mesh_path=mesh_path)

## current -> current density
group_current_dict = {
    0: 1,  # group_id : I (current,A)
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
plot_map(
    title='Magnetic Vector Potential[Wb/m]',
    vertInfoMat=vertInfoMat,
    c_mat=A_mat,
    boundary=(-1000, 1000),
    plot_type=PlotMapType.Coutour,
)
vertex_mat = vertInfoMat[trigInfoMat]
centroid_mat = calc_centroid(vertex_mat)
plot_map(
    title='Flux Density[T]',
    vertInfoMat=centroid_mat,
    c_mat=B_mat,
    boundary=(-1000, 1000),
    plot_type=PlotMapType.Scatter,
)
plot_map(
    title='Energy[W]', vertInfoMat=centroid_mat, c_mat=Energy_mat, boundary=(-1000, 1000), plot_type=PlotMapType.Scatter
)
