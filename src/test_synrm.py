import os
import re

import numpy as np
from plot import plot_map, PlotMapType

from utils import sum_area, calc_centroid
from logger import log
from solver import solve_magnetostatic
from read_mesh import read_nastran, read_gmsh

script_path = os.path.abspath(__file__)
root_path = re.search(r'.*(EM-FEM)', script_path).group(0)
data_path = os.path.join(root_path, 'data')

flag_compare = 0
## Read Mesh
if flag_compare:
    mesh_path = os.path.join(data_path, 'synrm', 'synrm_circular.msh')
    trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict = read_gmsh(
        mesh_path=mesh_path
    )

else:
    if 0:
        mesh_path = os.path.join(data_path, 'extraFine', 'synrm_extrafine.nas')
    else:
        mesh_path = os.path.join(data_path, 'coarse', 'synrm_coarse.nas')

    trigInfoMat, trigGroupMat, vertInfoMat, group_list, boundary_dict, material_in_use_dict = read_nastran(
        mesh_path=mesh_path
    )


## current -> current density
n_turn = 100
if flag_compare:
    winding_flag = 3
    if winding_flag == 0:
        group_current_dict = {  # group_id : I (current,A)
            0: 289.7777 * n_turn,
            1: 289.7777 * n_turn,
            2: 289.7777 * n_turn,
            3: 289.7777 * n_turn,
            4: 289.7777 * n_turn,
            5: 289.7777 * n_turn,
        }
        A_mat_femm_filename = 'A_mat_synrm_same.npy'
    elif winding_flag == 1:
        group_current_dict = {  # group_id : I (current,A)
            2: -212.1320 * n_turn,
            3: -212.1320 * n_turn,
        }
        A_mat_femm_filename = 'A_mat_synrm_23.npy'
    elif winding_flag == 2:
        group_current_dict = {  # group_id : I (current,A)
            4: 289.7777 * n_turn,
            5: 289.7777 * n_turn,
        }
        A_mat_femm_filename = 'A_mat_synrm_45.npy'
    elif winding_flag == 3:
        group_current_dict = {  # group_id : I (current,A)
            0: -77.6457 * n_turn,
            1: -77.6457 * n_turn,
            2: -212.1320 * n_turn,
            3: -212.1320 * n_turn,
            4: 289.7777 * n_turn,
            5: 289.7777 * n_turn,
        }
        A_mat_femm_filename = 'A_mat_synrm.npy'
else:
    group_current_dict = {  # group_id : I (current,A)
        5: 289.7777 * n_turn,
        9: 289.7777 * n_turn,
        12: -212.1320 * n_turn,
        14: -212.1320 * n_turn,
        16: -77.6457 * n_turn,
        18: -77.6457 * n_turn,
    }
    A_mat_femm_filename = 'A_mat_synrm_nastran.npy'

n_group = len(group_list)  # num of groups
group_current_density_list = np.zeros(n_group)  # group_id : J (current_density,A/m^2)

for group_id, current in group_current_dict.items():
    trig_ids = group_list[group_id]['trig_ids']  # group_dict['trig_ids']
    vertex_mat = vertInfoMat[trigInfoMat[trig_ids]]
    a_copper = sum_area(vertex_mat)  # sum the area of all triangles
    # a_copper = 1.8098e-4
    group_current_density_list[group_id] = current / a_copper

### Solve
S_mat, A_mat, T_mat, B_mat, B_norm_mat, Energy_mat = solve_magnetostatic(
    trigInfoMat=trigInfoMat,
    trigGroupMat=trigGroupMat,
    vertInfoMat=vertInfoMat,
    group_list=group_list,
    group_current_density_list=group_current_density_list,
    boundary_dict=boundary_dict,
    material_in_use_dict=material_in_use_dict,
    depth=0.0375,
)

### Post Process
total_Energy = np.sum(Energy_mat)
total_flux = np.max(np.abs(A_mat))
msg = f'Total Energy: {total_Energy}; Total Flux: {total_flux}'
log(msg, 'INFO')

### Plot

vertex_mat = vertInfoMat[trigInfoMat]
centroid_mat = calc_centroid(vertex_mat)

if 1 and flag_compare:
    A_mat_femm = np.load(os.path.join(data_path, 'synrm', A_mat_femm_filename))

    plot_map(
        title='Magnetic Vector Potential[Wb/m](FEMM)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat_femm,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Magnetic Vector Potential[Wb/m](Python)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat,
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Magnetic Vector Potential[Wb/m](Difference)',
        vertInfoMat=vertInfoMat,
        c_mat=A_mat - A_mat_femm,
        plot_type=PlotMapType.Coutourf,
    )

if 1:
    if not flag_compare:
        plot_map(
            title='Magnetic Vector Potential[Wb/m]',
            vertInfoMat=vertInfoMat,
            c_mat=A_mat,
            plot_type=PlotMapType.Coutourf,
        )
    boundary = 1e5
    plot_map(
        title='Flux Density[T]',
        vertInfoMat=centroid_mat,
        c_mat=B_mat,
        boundary=(-boundary, boundary),
        plot_type=PlotMapType.Quiver,
    )
    plot_map(
        title='Flux Density[T]',
        vertInfoMat=centroid_mat,
        c_mat=B_norm_mat,
        boundary=(-boundary, boundary),
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Energy[W]',
        vertInfoMat=centroid_mat,
        c_mat=Energy_mat,
        boundary=(-boundary, boundary),
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Matrix T[A]',
        vertInfoMat=vertInfoMat,
        c_mat=T_mat,
        boundary=(-boundary, boundary),
        plot_type=PlotMapType.Coutourf,
    )
    plot_map(
        title='Group',
        vertInfoMat=centroid_mat,
        c_mat=trigGroupMat,
        boundary=(-1, 100),
        plot_type=PlotMapType.Scatter,
    )
