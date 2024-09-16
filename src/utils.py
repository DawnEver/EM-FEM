import numpy as np

__all__ = ['sum_area', 'calc_centroid', 'gauss_seidel']


def calc_centroid(vertex_mat: np.ndarray) -> np.ndarray:
    """
    calculate the centroid coordinates of triangles
    >>> vertex_mat = np.array([[[1, 1], [2, 2], [3, 3]], [[4, 4], [5, 5], [6, 6]]])
    >>> np.all(np.isclose(calc_centroid(vertex_mat), np.array([[2., 2.],[5., 5.]])))
    np.True_
    """
    return np.sum(vertex_mat, axis=1) / 3


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


# def gauss_seidel(A, b, x0, tol=1e-4, max_iter=5):
#     """
#     Solve the system of linear equations using the Gauss-Seidel method.

#     Parameters:
#         A (numpy.ndarray): Coefficient matrix of shape (n, n).
#         b (numpy.ndarray): Right-hand side vector of shape (n,).
#         x0 (numpy.ndarray): Initial guess for the solution vector of shape (n,).
#         tol (float): Tolerance for convergence (default: 1e-4).
#         max_iter (int): Maximum number of iterations (default: 20).

#     Returns:
#         numpy.ndarray: Solution vector of shape (n,).
#     """
#     x = np.copy(x0)
#     L = np.tril(A)
#     U = A - L
#     for i_iter in range(max_iter+1):
#         print(f'\rGauss-Seidel itering: {i_iter}/{max_iter}', end='')
#         x_new = np.linalg.solve(L, b - np.dot(U, x))
#         error = np.linalg.norm(x_new - x)
#         if np.linalg.norm(x_new - x) < tol:
#             print('')
#             return x_new
#         x = x_new
#     print(f'\nNot converge. Error={error}/tol')
#     return x


def gauss_seidel(A, b, x0, tol=1e-6, max_iter=5):
    """
    Solve the system of linear equations using the Gauss-Seidel method.

    Parameters:
        A (numpy.ndarray): Coefficient matrix of shape (n, n).
        b (numpy.ndarray): Right-hand side vector of shape (n,).
        x0 (numpy.ndarray): Initial guess for the solution vector of shape (n,).
        tol (float): Tolerance for convergence (default: 1e-4).
        max_iter (int): Maximum number of iterations (default: 20).

    Returns:
        numpy.ndarray: Solution vector of shape (n,).
    """
    n = len(b)
    x = x0.copy()

    for i_iter in range(max_iter):
        print(f'\rGauss-Seidel itering: {i_iter+1}/{max_iter}', end='')
        x_new = np.zeros(n)

        for i in range(n):
            x_new[i] = (b[i] - np.dot(A[i, :i], x_new[:i]) - np.dot(A[i, i + 1 :], x[i + 1 :])) / A[i, i]
        error = np.linalg.norm(x - x_new)
        if error < tol:
            return x_new

        x = x_new.copy()
    print(f'\nNot converge: {error}/{tol}')
    return x


def _test():
    import doctest

    doctest.testmod()


if __name__ == '__main__':
    _test()
