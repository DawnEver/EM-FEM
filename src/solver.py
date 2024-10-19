import numpy as np
from scipy.sparse import coo_matrix, linalg as sp_linalg

from logger import log
from material import Material_In_Use

__all__ = ['solve_magnetostatic']


def solve_magnetostatic(
    trigInfoMat: np.ndarray,
    trigGroupMat: np.ndarray,
    vertInfoMat: np.ndarray,
    group_list: list,
    group_current_density_list: list,
    material_in_use_dict: dict,
    boundary_dict: dict,
    depth: float = 0.005,
):
    n_vert = len(vertInfoMat)
    n_trig = len(trigInfoMat)
    boundary_dirichlet = boundary_dict['dirichlet']
    boundary_sym = boundary_dict['sym']

    ### Selected material

    material_in_use = Material_In_Use(material_in_use_dict)

    ### loop
    A_mat = np.zeros(n_vert)
    B_mat = np.ones((n_trig, 2))
    B_norm_mat = np.ones(n_trig)
    Energy_mat = np.zeros(n_trig)
    last_A_mat = np.ones(n_vert)
    last_B_norm_mat = np.zeros(n_trig)

    rtol_B = 1e-2  # relative error tolerance of B
    rtol_A = 1e-2
    max_B_norm = 1e5  # 2.4  # consider saturation
    max_num_iter = 3
    reluctivity = 0

    is_convergence = False
    msg = f'S({n_vert, n_vert}) * A({n_vert}) = T({n_vert})'
    log(msg, 'INFO')
    for i_iter in range(max_num_iter):
        n_val = n_trig * 9
        T_mat = np.zeros(n_vert)
        S_row_mat = np.zeros(n_val, dtype=int)
        S_col_mat = np.zeros(n_val, dtype=int)
        S_val_mat = np.zeros(n_val)
        ### S_mat * A_mat = T_mat: prepare S_mat & T_mat
        for i_trig in range(n_trig):
            print(f'\rPrepare equations triangles: {i_trig+1}/{n_trig}', end='')
            if i_trig == n_trig - 1:
                print('')

            # get x,y
            # [?]TODO calculate these coordinates once and save as matrix (Space exchange time)
            vertex_ids = trigInfoMat[i_trig]
            [[xi, yi], [xj, yj], [xk, yk]] = vertInfoMat[vertex_ids]

            # ai = xj * yk - xk * yj
            bi = yj - yk
            ci = xk - xj
            # aj = xk * yi - xi * yk
            bj = yk - yi
            cj = xi - xk
            # ak = xi * yj - xj * yi
            bk = yi - yj
            ck = xj - xi
            delta = np.abs((xj * yk - xk * yj) + xi * (yj - yk) + yi * (xk - xj)) / 2  # area

            # calculate B
            [Ai, Aj, Ak] = A_mat[vertex_ids]
            current_B = np.array([bi * Ai + bj * Aj + bk * Ak, -ci * Ai - cj * Aj - ck * Ak]) / 2 / delta
            # magnetic flux density

            current_B_norm = np.linalg.norm(current_B)
            current_B_norm = min(current_B_norm, max_B_norm)
            if current_B_norm > max_B_norm:
                # scaler to max_B_norm
                scale = max_B_norm / current_B_norm
                current_B *= scale
            B_mat[i_trig] = current_B
            B_norm_mat[i_trig] = current_B_norm

            # get reluctivity of material
            group_id = int(trigGroupMat[i_trig])
            group_dict = group_list[group_id]
            material_type = group_dict['material_type']
            material_name = group_dict['material_name']
            reluctivity = material_in_use.get_reluctivity(material_name, B=current_B_norm)

            # calculate Energy
            Energy_mat[i_trig] = current_B_norm**2 * delta * depth * reluctivity / 2

            ## current
            if material_type == 'copper':
                # current
                current = group_current_density_list[group_id] * delta / 3
                T_mat[vertex_ids] += current

            Sii_raw = bi**2 + ci**2
            Sjj_raw = bj**2 + cj**2
            Skk_raw = bk**2 + ck**2
            Sij_raw = Sji_raw = bi * bj + ci * cj
            Sjk_raw = Skj_raw = bj * bk + cj * ck
            Sik_raw = Ski_raw = bi * bk + ci * ck

            start_id = 9 * i_trig
            end_id = start_id + 9
            S_val_mat[start_id:end_id] = (
                reluctivity
                / 4
                / delta
                * np.array([Sii_raw, Sij_raw, Sik_raw, Sji_raw, Sjj_raw, Sjk_raw, Ski_raw, Skj_raw, Skk_raw])
            )
            S_row_mat[start_id:end_id] = np.repeat(vertex_ids, 3)
            S_col_mat[start_id:end_id] = np.tile(vertex_ids, 3)

        if is_convergence:
            # last loop to update B_mat & A_mat
            break
        # current error < relative error -> exit loop
        if np.allclose(A_mat, last_A_mat, rtol=rtol_A) and np.allclose(B_norm_mat, last_B_norm_mat, rtol=rtol_B):
            is_convergence = True

        last_A_mat = A_mat
        last_B_norm_mat = B_norm_mat

        ### boundary condition
        # constant A
        for A_const, _vertex_ids in boundary_dirichlet.items():
            vertex_ids = np.unique(_vertex_ids)
            n_vertex_ids = len(vertex_ids)
            if n_vertex_ids < 1:
                continue
            T_mat[vertex_ids] = 3 * A_const

            # set 0
            S_row_ids = np.where(np.isin(S_row_mat, vertex_ids))
            S_val_mat[S_row_ids] = 0
            # set 1
            S_row_mat = np.concatenate((S_row_mat, vertex_ids))
            S_col_mat = np.concatenate((S_col_mat, vertex_ids))
            S_val_mat = np.concatenate((S_val_mat, np.ones(n_vertex_ids)))

        # symmetry boundary odd/even
        for boundary_sym_i in boundary_sym:
            tag, vertex_ids_1, vertex_ids_2 = boundary_sym_i
            if vertex_ids_1 is not None and vertex_ids_2 is not None:
                vertex_ids_1 = np.unique(vertex_ids_1)
                vertex_ids_2 = np.unique(vertex_ids_2)
                n_vertex_ids_1 = len(vertex_ids_1)
                n_vertex_ids_2 = len(vertex_ids_2)
                if n_vertex_ids_1 > 1 and n_vertex_ids_1 == n_vertex_ids_2:
                    vertex_ids = np.concatenate((vertex_ids_1, vertex_ids_2))
                    T_mat[vertex_ids] = 0
                    # set 0
                    S_row_ids = np.where(np.isin(S_row_mat, vertex_ids))
                    S_val_mat[S_row_ids] = 0
                    S_row_mat = np.concatenate((S_row_mat, vertex_ids))
                    S_col_mat = np.concatenate((S_col_mat, vertex_ids))
                    if tag == 'odd':
                        value = 1
                    else:
                        value = -1
                    S_val_mat = np.concatenate((S_val_mat, np.ones(n_vertex_ids_1), np.full(n_vertex_ids_2, value)))
                else:
                    msg = f'{tag} symmetry boundary conditions.\
                    len(vertex_ids_1)={n_vertex_ids_1} len(vertex_ids_2)={n_vertex_ids_2}'
                    log(msg=msg, level='ERROR')

        ### Solve A_mat: magnetic vector potential
        # https://www.osgeo.cn/scipy/reference/sparse.linalg.html#

        S_sp_mat = coo_matrix((S_val_mat, (S_row_mat, S_col_mat)), shape=(n_vert, n_vert))
        A_mat, _ = sp_linalg.cg(S_sp_mat, T_mat)

        log(f'Iteration: {i_iter+1}/{max_num_iter}', 'INFO')

    S_mat = S_sp_mat.toarray()
    return S_mat, A_mat, T_mat, B_mat, B_norm_mat, Energy_mat
