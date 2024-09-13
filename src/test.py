import numpy as np
import pygmsh

from logger import log

### Mesh
with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        mesh_size=0.5,
    )
    mesh = geom.generate_mesh()


vertInfoMat = mesh.points[:, 0:2]  # 2d fea

trigInfoMatRaw = mesh.cells_dict['triangle']
n_vert = len(vertInfoMat)
n_trig = len(trigInfoMatRaw)
group_ids = np.ones((n_trig, 1), dtype=np.uint64)
# mock
group_ids[[0, 1, 2]] = 2
group_ids[[-1, -2]] = 5

trigInfoMat = np.hstack((group_ids, trigInfoMatRaw))

### definitions

# TODO Material Database

miu0 = 4*np.pi*1e-7
miu_copper = 1 * miu0
miu_iron = 1000 * miu0
miu_air = 1 * miu0

r_stator_o = 0.075
r_stator_i = 0.04
r_rotor_o = 0.039
r_rotor_i = 0.0125
len_stack = 0.0375
a_copper = 1.8098e-4
n_turn = 100

domain_current_density_dict = {
    5: 289.7777 / a_copper/3,
    9: 289.7777 / a_copper/3,
    12: -212.1320 / a_copper/3,
    14: -212.1320 / a_copper/3,
    16: -77.6457 / a_copper/3,
    18: -77.6457 / a_copper/3,
} # J/3: divide 3 here

domain_copper = list(domain_current_density_dict.keys())
domain_iron = [1, 4, 8, 11, 13, 15, 17]
domain_air = [2, 3, 6, 7, 10]

A_constant_boundaries =[(0,[])] # [(A_const,[vertex_ids])]
sym_boundaries = [('odd',[],[])] # [('odd',[vertex_ids_1],[vertex_ids_2])]

### S_mat * A_mat = T_mat: prepare S_mat & T_mat

S_mat = np.zeros((n_vert, n_vert))
A_mat = np.zeros(n_vert)
T_mat = np.zeros(n_vert)
last_B_mat = np.ones(n_trig)
B_mat = np.zeros(n_trig)
last_Energy_mat = np.ones(n_trig)
Energy_mat = np.zeros(n_trig)

rtol_B = 1e-2 # relative error tolerance of B
rtol_Energy = 1e-2
max_num_iter = 10
miu = 0

is_convergence = False

for i_iter in range(max_num_iter+1):
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
        if group_id in domain_copper:
            # current 
            current = domain_current_density_dict[group_id] * delta
            T_mat[vertex_ids] += current
            
            miu = miu_copper
        elif group_id in domain_iron:
            # calculate miu
            miu = miu_iron
            # last_B = last_B_mat[i_trig]
            # miu = f(B)
        else:
            miu = miu_air
        
        # calculate B, Energy
        [Ai,Aj,Ak] = A_mat[vertex_ids]
        current_B = np.sqrt((bi*Ai + bj*Aj + bk*Ak)**2 + (ci*Ai + cj*Aj + ck*Ak)**2) /2 /delta  # magnetic flux density in the
        B_mat[i_trig] = current_B
        Energy_mat[i_trig] = current_B**2*delta*miu/2
        
        Sii_raw = bi ** 2 + ci ** 2
        Sjj_raw = bj ** 2 + cj ** 2
        Skk_raw = bk ** 2 + ck ** 2
        Sij_raw = Sji_raw = bi * bj + ci * cj
        Sjk_raw = Skj_raw = bj * bk + cj * ck
        Sik_raw = Ski_raw = bi * bk + ci * ck
        
        S_mat[np.ix_(vertex_ids, vertex_ids)] = np.array([
            [Sii_raw, Sij_raw, Sik_raw],
            [Sji_raw, Sjj_raw, Sjk_raw],
            [Ski_raw, Skj_raw, Skk_raw],]
        )*miu/4/delta

    if is_convergence:
        # last loop to update B_mat & Energy_mat
        break

    # current error < relative error -> exit loop
    if np.allclose(Energy_mat,last_Energy_mat,rtol = rtol_Energy) and np.allclose(B_mat,last_B_mat,rtol = rtol_B):
        is_convergence = True

    last_Energy_mat = Energy_mat
    last_B_mat = B_mat
    
    # boundary condition
    # constant A
    for boundary in A_constant_boundaries:
        A_const,vertex_ids = boundary
        n_vertex_ids = len(vertex_ids)
        if n_vertex_ids < 1:
            continue
        T_mat[vertex_ids] = 3 * A_const
        matrix = np.zeros((n_vertex_ids, n_vert))
        matrix[np.arange(n_vertex_ids), vertex_ids] = 1
        S_mat[vertex_ids] = matrix

    # symmetry boundary odd/even
    for boundary in sym_boundaries:
        symbol, vertex_ids_1, vertex_ids_2 = boundary
        n_vertex_ids_1 = len(vertex_ids_1)
        n_vertex_ids_2 = len(vertex_ids_2)
        if n_vertex_ids_1 < 1 or n_vertex_ids_1 != n_vertex_ids_2:
            continue
        if symbol == 'odd':
            sym = 1
        else: # even
            sym = -1
        T_mat[vertex_ids] = 0
        matrix = np.zeros((n_vertex_ids_1, n_vert))
        matrix[np.arange(n_vertex_ids_1), vertex_ids_1] = sym
        matrix[np.arange(n_vertex_ids_1), vertex_ids_2] = sym
        S_mat[vertex_ids_1] = matrix
        

    ### Solve
    A_mat = np.linalg.solve(S_mat,T_mat) # magnetic vector potential

    log(f'Current iteration: {i_iter}','INFO')

### Post Process
...

