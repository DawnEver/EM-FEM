import os
import re
import numpy as np

from utils import sum_area, calc_centroid
from logger import log
from solver import solve_magnetostatic
from read_mesh import read_gmsh, read_lipo_csv
from plot import plot_map, PlotMapType

script_path = os.path.abspath(__file__)
root_path = re.search(r'.*(EM-FEM)', script_path).group(0)
data_path = os.path.join(root_path, 'data')

## Read Mesh
FLAG_Lipo = False
if FLAG_Lipo:
    trigInfoPath = os.path.join(data_path, 'lipo', 'trigInfo.csv')
    vertInfoPath = os.path.join(data_path, 'lipo', 'vertInfo.csv')

    (trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict) = read_lipo_csv(
        trigInfoPath=trigInfoPath, vertInfoPath=vertInfoPath
    )
else:
    mesh_path = os.path.join(data_path, 'lipo', 'test-T_A_Lipo.msh')
    (trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict) = read_gmsh(
        mesh_path=mesh_path
    )

## current -> current density
group_current_dict = {
    2: 4000,  # group_id : I (current,A)
}
n_group = len(group_list)  # num of groups
group_current_density_list = np.zeros(n_group)  # group_id : J (current_density,A/m^2)

for group_id, current in group_current_dict.items():
    trig_ids = group_list[group_id]['trig_ids']  # group_dict['trig_ids']
    vertex_mat = vertInfoMat[trigInfoMat[trig_ids]]
    a_copper = sum_area(vertex_mat)  # sum the area of all triangles
    group_current_density_list[group_id] = current / a_copper

if 1:
    np.save(os.path.join(data_path, 'lipo', 'vertInfoMat_inductor.npy'), vertInfoMat)

### Solve
S_mat, A_mat, T_mat, B_mat, B_norm_mat, Energy_mat = solve_magnetostatic(
    trigInfoMat=trigInfoMat,
    trigGroupMat=trigGroupMat,
    vertInfoMat=vertInfoMat,
    group_list=group_list,
    group_current_density_list=group_current_density_list,
    boundary_dict=boundary_dict,
    material_in_use_dict=material_in_use_dict,
    depth=1,
)

### Post Process
total_Energy = np.sum(Energy_mat)
total_flux = np.max(np.abs(A_mat))
msg = f'Total Energy: {total_Energy}; Total Flux: {total_flux}'
log(msg, 'INFO')

### Plot
if 0:
    plot_map(
        title='Magnetic Vector Potential[Wb/m]',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat,
        plot_type=PlotMapType.Coutour,
    )
    plot_map(
        title='Matrix T[A]',
        vertInfoMat=vertInfoMat,
        c_mat=T_mat,
        plot_type=PlotMapType.Coutourf,
    )
    vertex_mat = vertInfoMat[trigInfoMat]
    centroid_mat = calc_centroid(vertex_mat)

    plot_map(
        title='Flux Density[T]',
        vertInfoMat=centroid_mat,
        c_mat=B_mat,
        plot_type=PlotMapType.Quiver,
    )
    plot_map(
        title='Flux Density[T]',
        vertInfoMat=centroid_mat,
        c_mat=B_norm_mat,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Energy[W]',
        vertInfoMat=centroid_mat,
        c_mat=Energy_mat,
        plot_type=PlotMapType.Coutourf,
    )
if 1 and FLAG_Lipo:
    file_path = os.path.join(root_path, 'mat-diff', 'data', 'A-matlab.csv')
    with open(file_path, encoding='utf-8') as f:
        A_mat_matlab = np.loadtxt(f)
    diff_A_mat = A_mat - A_mat_matlab
    plot_map(
        title="Error of Magnetic Vector Potential\nbetween Lipo's and My Python Code[Wb/m]",
        vertInfoMat=vertInfoMat,
        c_mat=diff_A_mat,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title="Magnetic Vector Potential[Wb/m](Lipo's Code)",
        vertInfoMat=vertInfoMat,
        c_mat=A_mat_matlab,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Magnetic Vector Potential[Wb/m](My Python Code)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat,
        plot_type=PlotMapType.Coutourf,
    )
elif 1 and not FLAG_Lipo:
    file_path = os.path.join(data_path, 'lipo', 'A_mat_inductor.npy')
    A_mat_femm = np.load(file_path)
    print(max(A_mat_femm))
    diff_A_mat = A_mat - A_mat_femm
    plot_map(
        title='Error of Magnetic Vector Potential\nbetween FEMM and My Pyhton Code[Wb/m]',
        vertInfoMat=vertInfoMat,
        c_mat=diff_A_mat,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Magnetic Vector Potential[Wb/m](FEMM)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat_femm,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Magnetic Vector Potential[Wb/m](My Python Code)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat,
        plot_type=PlotMapType.Coutourf,
    )
