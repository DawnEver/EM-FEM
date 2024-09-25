import numpy as np

__all__ = ['sum_area', 'calc_centroid', 'gauss_seidel_sor', 'jacobi_iteration']


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


def gauss_seidel_sor(A, b, x0=None, omega=1, tol=1e-6, max_iter=100, is_print: bool = True):
    """
    Solve the system of linear equations using the Gauss-Seidel method with Successive Over-Relaxation(SOR).

    Parameters:
        A (numpy.ndarray): Coefficient matrix of shape (n, n).
        b (numpy.ndarray): Right-hand side vector of shape (n,).
        x0 (numpy.ndarray): Initial guess for the solution vector of shape (n,).
        omega (float): Relaxation parameter (0 < omega < 2).
        tol (float): Tolerance for convergence (default: 1e-4).
        max_iter (int): Maximum number of iterations (default: 20).

    Returns:
        numpy.ndarray: Solution vector of shape (n,).

    >>> A = np.array([[4, -1, 0], [-1, 4, -1], [0, -1, 4]])
    >>> b = np.array([9, 5, 0])
    >>> np.allclose(gauss_seidel(A, b, max_iter=100,is_print=False), np.linalg.solve(A, b), atol=1e-6)
    True
    """
    n = len(b)
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.copy()

    for i_iter in range(max_iter):
        x_new = np.copy(x)
        for i in range(n):
            sum1 = np.dot(A[i, :i], x_new[:i])
            sum2 = np.dot(A[i, i + 1 :], x[i + 1 :])
            x_new[i] = (1 - omega) * x[i] + (omega / A[i, i]) * (b[i] - sum1 - sum2)

        error = np.linalg.norm(x_new - x)
        if is_print:
            print(f'\rGauss-Seidel SOL itering: {i_iter+1}/{max_iter} | error:{error}', end='')
        if error < tol:
            break
        x = x_new

    if is_print:
        print(f'\nNot converge: {error}/{tol}')
    return x


def jacobi_iteration(A, b, x0=None, tol=1e-6, max_iter=100, is_print: bool = True):
    """
    Solve the system of linear equations using the Jacobi iteration method.

    Parameters:
        A (numpy.ndarray): Coefficient matrix of shape (n, n).
        b (numpy.ndarray): Right-hand side vector of shape (n,).
        x0 (numpy.ndarray): Initial guess for the solution vector of shape (n,).
        tol (float): Tolerance for convergence (default: 1e-4).
        max_iter (int): Maximum number of iterations (default: 20).

    Returns:
        numpy.ndarray: Solution vector of shape (n,).

    >>> A = np.array([[4, -1, 0], [-1, 4, -1], [0, -1, 4]])
    >>> b = np.array([9, 5, 0])
    >>> np.allclose(gauss_seidel(A, b, max_iter=100,is_print=False), np.linalg.solve(A, b), atol=1e-6)
    True
    """
    n = len(b)
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.copy()
    A_diag_mat = np.diagonal(A)
    T = A - np.diag(A_diag_mat)

    for i_iter in range(max_iter):
        x_new = (b - np.dot(T, x)) / A_diag_mat
        error = np.linalg.norm(x_new - x)
        if is_print:
            print(f'\r Jacobi itering: {i_iter+1}/{max_iter} | error:{error}', end='')
        if error < tol:
            break
        x = x_new

    if is_print:
        print(f'\nNot converge: {error}/{tol}')
    return x


def _test():
    import doctest

    doctest.testmod()


if __name__ == '__main__':
    _test()
