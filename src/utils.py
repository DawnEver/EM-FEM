import numpy as np

__all__ = ['sum_area']


def sum_area(vertex_mat: np.ndarray) -> np.number:
    """
    # sum the area of a series of triangles

    >>> vertex_mat = np.array([[[1, 0], [0, 1], [0, 0]], [[1, 1], [0, 1], [0, 0]]])
    >>> np.isclose(sum_area(vertex_mat), 1.0)
    np.True_

    >>> vertex_mat = np.array([[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]])
    >>> np.isclose(sum_area(vertex_mat), 0.0)
    np.True_

    >>> vertex_mat = np.array([[[1, 1], [2, 2], [3, 3]], [[4, 4], [5, 5], [6, 6]]])
    >>> np.isclose(sum_area(vertex_mat), 1.0)
    np.False_

    >>> vertex_mat = np.array([[[0.,1.],[0.25,1.25],[0.5,1.]],[[1.,1.],[1.25,0.75],[1.,0.5]], \
     [[1.5,1.5],[1.25,1.25],[1.,1.5]],[[1.5,0.],[1.25,0.25],[1.5,0.5]], \
     [[0.5,1.5],[0.75,1.25],[0.5,1.]],[[0.5,1.],[0.25,1.25],[0.5,1.5]], \
     [[1.,1.],[1.25,1.25],[1.5,1.]],[[1.5,1.],[1.25,0.75],[1.,1.]], \
     [[1.,1.5],[0.75,1.25],[0.5,1.5]],[[1.25,0.75],[1.25,0.25],[1.,0.5]], \
     [[0.,1.5],[0.25,1.25],[0.,1.]],[[1.5,0.5],[1.25,0.25],[1.25,0.75]], \
     [[1.,0.],[1.25,0.25],[1.5,0.]],[[1.5,0.5],[1.25,0.75],[1.5,1.]], \
     [[1.5,1.],[1.25,1.25],[1.5,1.5]],[[1.,1.],[0.75,1.25],[1.25,1.25]], \
     [[1.,0.5],[1.25,0.25],[1.,0.]],[[0.5,1.5],[0.25,1.25],[0.,1.5]], \
     [[0.5,1.],[0.75,1.25],[1.,1.]],[[1.25,1.25],[0.75,1.25],[1.,1.5]]])
    >>> np.isclose(sum_area(vertex_mat), 1.25)
    np.True_
    """
    # extend each matrix to a square matrix (3x3)
    vertex_mat_extended = np.concatenate((vertex_mat, np.ones((vertex_mat.shape[0], vertex_mat.shape[1], 1))), axis=2)
    # S = sum( 1/2 * det([[xi,yi,1],[xj,yj,1],[xk,yk,1]]) )
    return np.sum(np.abs(np.linalg.det(vertex_mat_extended))) / 2


def _test():
    import doctest

    doctest.testmod()


if __name__ == '__main__':
    _test()
