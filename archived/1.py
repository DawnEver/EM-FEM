import numpy as np
import pygmsh

### Model
model = [
    {
        'shape': 'polygon',  #
        'points': [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        'material_type': 'copper',
        'material_name': 'Cu',
        'mesh_size': 0.5,
    },
    {
        'shape': 'polygon',
        'points': [
            [1.1, 0.0],
            [1.1, 1.1],
            [0.0, 1.1],
            [0.0, 2.0],
            [2.0, 2.0],
            [2.0, 0.0],
        ],
        'material_type': 'steel',
        'material_name': 'iron',
        'mesh_size': 0.5,
    },
    {
        'shape': 'polygon',
        'points': [
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.1],
            [1.1, 1.1],
            [1.1, 0.0],
        ],
        'material_type': 'other',
        'material_name': 'air',
        'mesh_size': 0.05,
    },
]


### Mesh
floor_id = 0
with pygmsh.geo.Geometry() as geom:
    for _, section in enumerate(model):
        shape = section['shape']
        points = section['points']
        mesh_size = section['mesh_size']
        if shape == 'polygon':
            loopPath = geom.add_polygon(points, mesh_size=mesh_size)
        else:
            ...

    mesh = geom.generate_mesh()

    # this_vert_mat = mesh.points[:, 0:2]  # 2d fea
    # len_this_vert_mat = len(this_vert_mat)
    # this_trig_mat = mesh.cells_dict['triangle'] + floor_id
    # len_this_trig_mat = len(this_trig_mat)
    # group_ids = np.full((len_this_trig_mat, 1), i_section, dtype=np.uint64)

    # floor_id += len_this_vert_mat
    # # update the vertex id floor
    # # 0 ~ n --> 0 + sum(len_this_vert_mat) ~ n + sum(len_this_vert_mat)


vertInfoMat = mesh.points[:, 0:2]  # 2d fea

trigInfoMatRaw = mesh.cells_dict['triangle']
n_vert = len(vertInfoMat)
n_trig = len(trigInfoMatRaw)

i_section = 0
group_ids = np.full((n_trig, 1), i_section, dtype=np.uint64)
trigInfoMat = np.hstack((group_ids, trigInfoMatRaw))


lines = mesh.cells_dict['line']
