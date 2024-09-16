import numpy as np

from logger import log
from material import Material_In_Use
from utils import gauss_seidel

__all__ = ['solve_magnetostatic']


def solve_magnetostatic(
    trigInfoMat: np.ndarray,
    trigGroupMat: np.ndarray,
    vertInfoMat: np.ndarray,
    group_list: list,
    group_current_density_list: list,
    material_in_use_dict: dict,
    boundary_dict: dict,
):
    n_vert = len(vertInfoMat)
    n_trig = len(trigInfoMat)
    boundary_dirichlet = boundary_dict['dirichlet']
    boundary_sym = boundary_dict['sym']

    ### Selected material

    material_in_use = Material_In_Use(material_in_use_dict)

    ### loop
    S_mat = np.zeros((n_vert, n_vert))
    A_mat = np.zeros(n_vert)
    T_mat = np.zeros(n_vert)
    B_mat = np.ones(n_trig)
    Energy_mat = np.zeros(n_trig)
    last_A_mat = np.ones(n_vert)
    last_B_mat = np.zeros(n_trig)

    rtol_B = 1e-2  # relative error tolerance of B
    rtol_A = 1e-2
    max_B = 2.4  # consider saturation
    max_num_iter = 5
    reluctivity = 0

    is_convergence = False
    msg = f'S({n_vert, n_vert}) * A({n_vert}) = T({n_vert})'
    log(msg, 'INFO')
    for i_iter in range(max_num_iter):
        ### S_mat * A_mat = T_mat: prepare S_mat & T_mat
        for i_trig in range(n_trig):
            print(f'\rPrepare equations triangles: {i_trig}/{n_trig}', end='')
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
            current_B = (
                np.sqrt((bi * Ai + bj * Aj + bk * Ak) ** 2 + (ci * Ai + cj * Aj + ck * Ak) ** 2) / 2 / delta
            )  # magnetic flux density in the
            current_B = min(current_B, max_B)
            B_mat[i_trig] = current_B

            # get reluctivity of material
            group_id = int(trigGroupMat[i_trig])
            group_dict = group_list[group_id]
            material_type = group_dict['material_type']
            material_name = group_dict['material_name']
            reluctivity = material_in_use.get_reluctivity(material_name, B=current_B)

            # calculate Energy
            Energy_mat[i_trig] = current_B**2 * delta * reluctivity / 2

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

            S_mat[np.ix_(vertex_ids, vertex_ids)] += (
                np.array(
                    [
                        [Sii_raw, Sij_raw, Sik_raw],
                        [Sji_raw, Sjj_raw, Sjk_raw],
                        [Ski_raw, Skj_raw, Skk_raw],
                    ]
                )
                * reluctivity
                / 4
                / delta
            )

        if is_convergence:
            # last loop to update B_mat & A_mat
            break
        # current error < relative error -> exit loop
        if np.allclose(A_mat, last_A_mat, rtol=rtol_A) and np.allclose(B_mat, last_B_mat, rtol=rtol_B):
            is_convergence = True

        last_A_mat = A_mat
        last_B_mat = B_mat

        ### boundary condition
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
        for boundary_sym_i in boundary_sym:
            tag, vertex_ids_1, vertex_ids_2 = boundary_sym_i
            if vertex_ids_1 is not None and vertex_ids_2 is not None:
                n_vertex_ids_1 = len(vertex_ids_1)
                n_vertex_ids_2 = len(vertex_ids_2)
                if n_vertex_ids_1 > 1 and n_vertex_ids_1 != n_vertex_ids_2:
                    T_mat[vertex_ids] = 0
                    matrix = np.zeros((n_vertex_ids_1, n_vert))
                    matrix[np.arange(n_vertex_ids_1), vertex_ids_1] = 1
                    if tag == 'odd':
                        value = 1
                    else:
                        value = -1
                    matrix[np.arange(n_vertex_ids_1), vertex_ids_2] = value
                    S_mat[vertex_ids_1] = matrix

        ### Solve A_mat: magnetic vector potential
        # LU decomposition/Gaussian Elimination
        # TODO Gauss-Seidel Iteration
        flag = 1
        match flag:
            case 0:
                # A_mat = np.linalg.solve(S_mat, T_mat)
                # A_mat = np.linalg.tensorsolve(S_mat, T_mat)
                A_mat = np.linalg.inv(S_mat).dot(T_mat)
            case _:
                A_mat = gauss_seidel(S_mat, T_mat, A_mat)

        log(f'Iteration: {i_iter+1}/{max_num_iter}', 'INFO')
    return S_mat, A_mat, T_mat, B_mat, Energy_mat
