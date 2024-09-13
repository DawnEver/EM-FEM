import numpy as np
import matplotlib.pyplot as plt

### Read NAS

### Read CSV
if 1:  # coarse mesh
    trigInfoPath = './coarse/trigInfo.csv'
    vertInfoPath = './coarse/vertInfo.csv'
else:  # extra fine mesh
    trigInfoPath = './extraFine/trigInfo.csv'
    vertInfoPath = './extraFine/vertInfo.csv'

with open(trigInfoPath, encoding='utf-8') as f:
    trigInfoMat = np.loadtxt(f, delimiter=',', dtype=np.int16)
with open(vertInfoPath, encoding='utf-8') as f:
    vertInfoMat = np.loadtxt(f, delimiter=',', dtype=np.float16)


### definitions
permeability_copper = 1
permeability_iron = 1000
permeability_air = 1

domain_copper = [5, 9, 12, 14, 16, 18]
domain_iron = [1, 4, 8, 11, 13, 15, 17]
domain_air = [2, 3, 6, 7, 10]

r_stator_o = 0.075
r_stator_i = 0.04
r_rotor_o = 0.039
r_rotor_i = 0.0125
len_stack = 0.0375
a_copper = 1.8098e-4
n_turn = 100


### Plot Mesh
if 1:
    colors = ['r', 'g', 'b', 'y', 'black', 'brown', 'grey']
    num_colors = len(colors)

    vertInfoMat = vertInfoMat * 1e4
    # trigInfoMat = trigInfoMat[:20]

    plt.figure()
    for row in trigInfoMat:
        id_group = row[1]
        # if id_group not in [3]:
        #     continue
        ids_vert = row[2:8]
        color = colors[id_group % num_colors]
        points = vertInfoMat[ids_vert]
        plt.plot(points[:, 0], points[:, 1], color=color)
        plt.scatter(points[:, 0], points[:, 1], color=color)
        # break
        plt.show()

    # plt.scatter(vertInfoMat[:,0],vertInfoMat[:,1])
    # plt.show()

### Simulate
