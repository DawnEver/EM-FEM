import argparse
import os
import re

import numpy as np
import matplotlib.pyplot as plt

script_path = os.path.abspath(__file__)
root_path = re.search(r'.*(mat-diff)', script_path).group(0)
data_path = os.path.join(root_path, 'data')


def read_file(file_path: str):
    # TODO if not csv
    with open(file_path, encoding='utf-8') as f:
        matrix = np.loadtxt(f, delimiter=',')
    return matrix


def show_diff(mat_diff: np.ndarray, mat_1: np.ndarray, mat_2: np.ndarray):
    # TODO show diff in excel or html
    err_tol = 1e-5
    mat_is_diff = np.isclose(mat_1, mat_2)
    count = np.bincount(mat_is_diff.flatten())
    num_F = count[0]
    num_T = count[1]
    msg = f'Max Error: {mat_diff.max()}/{mat_diff.min()}, Error Tolerance: {err_tol},\nUnequal: {num_F}/{num_F+num_T},'
    print(msg)
    if len(mat_diff.shape) > 1:
        if 0:
            fig, axes = plt.subplots(1, 2)
            ax0 = axes[0].imshow(mat_1)
            fig.colorbar(ax0)
            ax1 = axes[1].imshow(mat_diff)
            fig.colorbar(ax1)
        else:
            plt.imshow(mat_diff)
            plt.colorbar()
        plt.show()


def mat_diff(file_path_1: str, file_path_2: str):
    mat_1 = read_file(file_path_1)
    mat_2 = read_file(file_path_2)
    mat_diff = mat_1 - mat_2
    show_diff(mat_diff, mat_1, mat_2)


if __name__ == '__main__':
    # get args
    if 0:
        parser = argparse.ArgumentParser()
        parser.add_argument('file_path_1')
        parser.add_argument('file_path_2')
        args = parser.parse_args()
        file_path_1 = args.file_path_1
        file_path_2 = args.file_path_2
    else:
        flag = 2
        if flag == 0:
            f1 = 'T-matlab.csv'
            f2 = 'T-python.csv'
        elif flag == 1:
            f1 = 'A-matlab.csv'
            f2 = 'A-python.csv'
        else:
            f1 = 'S-python-v1.csv'
            # f2 = 'S-python-v0.csv'
            f2 = 'S-matlab.csv'
        file_path_1 = os.path.join(data_path, f1)
        file_path_2 = os.path.join(data_path, f2)
    mat_diff(file_path_1=file_path_1, file_path_2=file_path_2)
