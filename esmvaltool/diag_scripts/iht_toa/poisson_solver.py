"""Poisson solver for the full ocean-atmosphere column. The Poisson equation is
solved by numerically using the biconjugate gradient stabilized (BiCGSTAB)
method.

The solution is achieved when the difference between the input field (radiative
flux) and the Laplacian of the output field is less than the stated tolerance.
If the solver fails to converge, the tolerance can be increased.

Convergence is achieved faster by using a preconditioner on the output field.

The heat transport is calculated as the gradient of the scalar
p-field output of the Poisson solver.
"""

import numpy as np
from numba import jit


def set_metrics():
    # Define variables used in the solver
    hp = np.zeros(M)
    hv = np.zeros(M + 1)

    y = -0.5 * np.pi + 0.5 * dy
    hv[0] = 0.0
    for j in range(0, M):
        hp[j] = np.cos(y)
        hv[j + 1] = np.cos(y + 0.5 * dy)
        y += dy
    hv[-1] = 0.0

    return hp, hv


def set_matrix(hp, hv):
    # Storing the full matrix
    A_p = np.zeros([M, N])
    A_e = np.zeros([M, N])
    A_w = np.zeros([M, N])
    A_s = np.zeros([M, N])
    A_n = np.zeros([M, N])

    # ILU factors
    M_p = np.zeros([M + 1, N + 1])
    M_e = np.zeros([M + 1, N + 1])
    M_w = np.zeros([M + 1, N + 1])
    M_s = np.zeros([M + 1, N + 1])
    M_n = np.zeros([M + 1, N + 1])

    # Spherical Laplacian variables
    A = 1.0 / (dx**2.)
    B = 1.0 / (dy**2.)

    # First calculate the Poisson equations 5-point stencil
    # A_w is the contribution from i-1, A_e is from i+1,
    # A_s is j-1, A_n is j+1, and A_p is the diagonal
    for j in range(0, M):
        tx = A / hp[j]**2.0
        ty = B / hp[j]

        for i in range(0, N):
            A_w[j, i] = tx
            A_e[j, i] = tx
            A_s[j, i] = ty * hv[j]
            A_n[j, i] = ty * hv[j + 1]
            A_p[j, i] = -(A_w[j, i] + A_e[j, i] + A_s[j, i] + A_n[j, i])

    # ILU/SIP preconditioner factors: alf = 0.0 is ILU
    alf = 0.9
    M_p += 1.0

    for j in range(1, M + 1):
        for i in range(1, N + 1):
            M_s[j, i] = A_s[j - 1, i - 1] / (1.0 + alf * M_e[j - 1, i])
            M_w[j, i] = A_w[j - 1, i - 1] / (1.0 + alf * M_n[j, i - 1])

            M_p[j, i] = A_p[j-1, i-1] - \
                M_s[j, i]*(M_n[j-1, i] - alf*M_e[j-1, i]) -\
                M_w[j, i]*(M_e[j, i-1] - alf*M_n[j, i-1])
            M_p[j, i] = 1.0 / M_p[j, i]

            M_e[j, i] = (A_e[j - 1, i - 1] -
                         alf * M_s[j, i] * M_e[j - 1, i]) * M_p[j, i]
            M_n[j, i] = (A_n[j - 1, i - 1] -
                         alf * M_w[j, i] * M_n[j, i - 1]) * M_p[j, i]

    return A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s, M_n, M_p


def swap_bounds(fld):
    # Extends the array by one in all directions
    # As the array is periodic it allows for easier computations at boundaries
    wrap_pnt = int(N / 2 + 1)
    for i in range(1, N + 1):
        fld[0, i] = fld[1, wrap_pnt]
        fld[M + 1, i] = fld[M, wrap_pnt]
        wrap_pnt += 1
        if wrap_pnt > N:
            wrap_pnt = 1

    fld[:, 0] = fld[:, N]
    fld[:, N + 1] = fld[:, 1]

    return fld


def calc_Ax(x, A_e, A_w, A_s, A_n, A_p):
    # Laplacian equation
    Ax = np.zeros([M + 2, N + 2])

    x = swap_bounds(x)
    # Ax[j, i] = A_s[j-1, i-1]*x[j-1, i] + A_w[j-1, i-1]*x[j, i-1] + \
    #     A_e[j-1, i-1]*x[j, i+1] + A_n[j-1, i-1]*x[j+1, i] + \
    #     A_p[j-1, i-1]*x[j, i]
    Ax[1:M+1, 1:N+1] = A_s[0:M, 0:N] * x[0:M, 1:N+1] + A_w[0:M, 0:N] * \
        x[1:M+1, 0:N] + A_e[0:M, 0:N] * x[1:M+1, 2:N+2] + A_n[0:M, 0:N] * \
        x[2:M+2, 1:N+1] + A_p[0:M, 0:N] * x[1:M+1, 1:N+1]
    Ax = swap_bounds(Ax)
    return Ax


def dot_prod(x, y):
    # Calculate dot product of two matrices
    return (x[1:M + 1, 1:N + 1] * y[1:M + 1, 1:N + 1]).sum()


def precon(x, M_e, M_w, M_s, M_n, M_p):
    # Preconditioner
    Cx = np.zeros([M + 2, N + 2])
    precon_a(x, M_w, M_s, M_p, Cx)
    Cx = swap_bounds(Cx)
    precon_b(M_e, M_n, Cx)
    Cx = swap_bounds(Cx)
    return Cx


@jit
def precon_a(x, M_w, M_s, M_p, Cx):
    for j in range(1, M + 1):
        for i in range(1, N + 1):
            Cx[j, i] = M_p[j, i] * (x[j, i] - M_s[j, i] * Cx[j - 1, i] -
                                    M_w[j, i] * Cx[j, i - 1])


@jit
def precon_b(M_e, M_n, Cx):
    for j in range(M, 0, -1):
        for i in range(N, 0, -1):
            Cx[j,
               i] = Cx[j,
                       i] - M_e[j, i] * Cx[j, i + 1] - M_n[j, i] * Cx[j + 1, i]


def bicgstab(logger, x, b, A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s, M_n, M_p):
    # Bi-conjugate gradient stabilized numerical solver
    sc_err = dot_prod(b, b)

    Ax = calc_Ax(x, A_e, A_w, A_s, A_n, A_p)

    r = b - Ax
    cr = r

    alf = 1.0
    omg = 1.0
    nrm = 1.0
    imx = 1000

    p = np.zeros([M + 2, N + 2])
    v = np.zeros([M + 2, N + 2])

    it = 0
    while it < imx:
        rho = dot_prod(r, cr)

        bet = (rho / nrm) * (alf / omg)

        t = r - bet * omg * v

        s = precon(t, M_e, M_w, M_s, M_n, M_p)
        p = s + bet * p

        v = calc_Ax(p, A_e, A_w, A_s, A_n, A_p)
        nrm = dot_prod(cr, v)

        alf = rho / nrm
        s = r - alf * v

        cs = precon(s, M_e, M_w, M_s, M_n, M_p)
        t = calc_Ax(cs, A_e, A_w, A_s, A_n, A_p)

        tt = dot_prod(t, t)
        ts = dot_prod(t, s)
        omg = ts / tt

        x = x + alf * p + omg * cs
        r = s - omg * t

        nrm = rho

        if abs(omg) < 1.0e-16:
            logger.info('Terminating Poisson solver.')
            quit()

        err = np.sqrt(dot_prod(r, r) / sc_err)
        if err < tol:
            break

        it += 1

    return x


def calc_mht(sol):
    # MHT = Grad(P)
    y = np.arange(-0.5 * np.pi + 0.5 * dy, 0.5 * np.pi, dy)

    grad_phi = np.gradient(sol, dy, axis=0)
    grad_phi = grad_phi[1:-1, 1:-1]
    MHT = np.sum((grad_phi.T * np.cos(y) * dx).T, axis=1)

    return MHT


def spherical_poisson(logger, forcing, tolerance=None):
    # Define global variables
    logger.info("spherical_poisson: Setting global variables")

    global N, M, dx, dy, tol
    (M, N) = (np.shape(forcing))
    N = int(N)
    M = int(M)
    dx = 2.0 * np.pi / N
    dy = np.pi / M
    if tolerance is None:
        tol = 2.0e-4
    else:
        tol = tolerance

    logger.info("spherical_poisson: Calling set_metrics")
    hp, hv = set_metrics()

    logger.info("spherical_poisson: Calling set_matrix")
    A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s, M_n, M_p = set_matrix(hp, hv)

    # Solving the Poisson equation
    rhs = np.zeros([M + 2, N + 2])
    sol = np.zeros([M + 2, N + 2])
    logger.info("spherical_poisson: Solving poisson equation")
    rhs[1:-1, 1:-1] = forcing
    rhs = swap_bounds(rhs)
    sol = bicgstab(logger, sol, rhs, A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s,
                   M_n, M_p)

    # Calculating meridional heat transport
    logger.info("spherical_poisson: Calculating MHT")
    mht = calc_mht(sol)

    return sol, mht
