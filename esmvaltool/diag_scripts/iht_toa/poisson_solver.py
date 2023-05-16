"""(C) Crown Copyright 2023, the Met Office. Poisson solver for the full ocean-
atmosphere column. The Poisson equation is solved by numerically using the
biconjugate gradient stabilized (BiCGSTAB) method.

The solution is achieved when the difference between the input field (radiative
flux) and the Laplacian of the output field is less than the stated tolerance.
If the solver fails to converge, the tolerance can be increased.

Convergence is achieved faster by using a preconditioner on the output field.

The heat transport is calculated as the gradient of the energy flux potential,
the output of the Poisson solver.
"""

import numpy as np

# from numba import jit


def set_metrics():
    """Define variables used in the solver."""
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
    """Calculate the A-matrix (Eq.

    8) that defines the five-point stencil. The A_[s,n,w,s] are the
    values are the contributions from each of the four neighbouring
    cells, while A_p is the contribution from the given cell.
    """
    # Storing the full matrix
    a_p = np.zeros([M, N])
    a_e = np.zeros([M, N])
    a_w = np.zeros([M, N])
    a_s = np.zeros([M, N])
    a_n = np.zeros([M, N])

    # ILU factors
    m_p = np.zeros([M + 1, N + 1])
    m_e = np.zeros([M + 1, N + 1])
    m_w = np.zeros([M + 1, N + 1])
    m_s = np.zeros([M + 1, N + 1])
    m_n = np.zeros([M + 1, N + 1])

    # Spherical Laplacian variables
    a = 1.0 / (dx**2.)
    b = 1.0 / (dy**2.)

    # First calculate the Poisson equations 5-point stencil
    # A_w is the contribution from i-1, A_e is from i+1,
    # A_s is j-1, A_n is j+1, and A_p is the diagonal
    for j in range(0, M):
        tx = a / hp[j]**2.0
        ty = b / hp[j]

        for i in range(0, N):
            a_w[j, i] = tx
            a_e[j, i] = tx
            a_s[j, i] = ty * hv[j]
            a_n[j, i] = ty * hv[j + 1]
            a_p[j, i] = -(a_w[j, i] + a_e[j, i] + a_s[j, i] + a_n[j, i])

    # ILU/SIP preconditioner factors: alf = 0.0 is ILU
    alf = 0.9
    m_p += 1.0

    for j in range(1, M + 1):
        for i in range(1, N + 1):
            m_s[j, i] = a_s[j - 1, i - 1] / (1.0 + alf * m_e[j - 1, i])
            m_w[j, i] = a_w[j - 1, i - 1] / (1.0 + alf * m_n[j, i - 1])

            m_p[j, i] = a_p[j-1, i-1] - \
                m_s[j, i]*(m_n[j-1, i] - alf*m_e[j-1, i]) -\
                m_w[j, i]*(m_e[j, i-1] - alf*m_n[j, i-1])
            m_p[j, i] = 1.0 / m_p[j, i]

            m_e[j, i] = (a_e[j - 1, i - 1] -
                         alf * m_s[j, i] * m_e[j - 1, i]) * m_p[j, i]
            m_n[j, i] = (a_n[j - 1, i - 1] -
                         alf * m_w[j, i] * m_n[j, i - 1]) * m_p[j, i]

    return a_e, a_w, a_s, a_n, a_p, m_e, m_w, m_s, m_n, m_p


def swap_bounds(fld):
    """Extends the array by one in all directions.

    As the array is periodic it allows for easier computations at
    boundaries.
    """
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


def calc_ax(x, a_e, a_w, a_s, a_n, a_p):
    """Matrix calculation of the Laplacian equation, LHS of Eq.

    (9) in Pearce and Bodas-Salcedo (2023).
    """
    # Laplacian equation
    ax = np.zeros([M + 2, N + 2])

    x = swap_bounds(x)
    # Ax[j, i] = A_s[j-1, i-1]*x[j-1, i] + A_w[j-1, i-1]*x[j, i-1] + \
    #     A_e[j-1, i-1]*x[j, i+1] + A_n[j-1, i-1]*x[j+1, i] + \
    #     A_p[j-1, i-1]*x[j, i]
    ax[1:M+1, 1:N+1] = a_s[0:M, 0:N] * x[0:M, 1:N+1] + a_w[0:M, 0:N] * \
        x[1:M+1, 0:N] + a_e[0:M, 0:N] * x[1:M+1, 2:N+2] + a_n[0:M, 0:N] * \
        x[2:M+2, 1:N+1] + a_p[0:M, 0:N] * x[1:M+1, 1:N+1]
    ax = swap_bounds(ax)
    return ax


def dot_prod(x, y):
    """Calculate dot product of two matrices."""
    return (x[1:M + 1, 1:N + 1] * y[1:M + 1, 1:N + 1]).sum()


def precon(x, m_e, m_w, m_s, m_n, m_p):
    """Preconditioner.

    This is a wrapper to two steps that are optimised using jit.
    """
    cx = np.zeros([M + 2, N + 2])
    precon_a(x, m_w, m_s, m_p, cx)
    cx = swap_bounds(cx)
    precon_b(m_e, m_n, cx)
    cx = swap_bounds(cx)
    return cx


# @jit
def precon_a(x, m_w, m_s, m_p, cx):
    """First step of preconditioner."""
    for j in range(1, M + 1):
        for i in range(1, N + 1):
            cx[j, i] = m_p[j, i] * (x[j, i] - m_s[j, i] * cx[j - 1, i] -
                                    m_w[j, i] * cx[j, i - 1])


# @jit
def precon_b(m_e, m_n, cx):
    """Second step of preconditioner."""
    for j in range(M, 0, -1):
        for i in range(N, 0, -1):
            cx[j,
               i] = cx[j,
                       i] - m_e[j, i] * cx[j, i + 1] - m_n[j, i] * cx[j + 1, i]


def bicgstab(logger, x, b, a_e, a_w, a_s, a_n, a_p, m_e, m_w, m_s, m_n, m_p):
    """Bi-conjugate gradient stabilized numerical solver.

    van der Vorst, H. A., 1992: Bi-cgstab: A fast and smoothly
    converging variant of bi-cg for the solution of nonsymmetric linear
    systems. SIAM Journal on Scientific and Statistical Computing,
    https://doi.org/10.1137/0913035.
    """
    sc_err = dot_prod(b, b)

    ax = calc_ax(x, a_e, a_w, a_s, a_n, a_p)

    r = b - ax
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

        s = precon(t, m_e, m_w, m_s, m_n, m_p)
        p = s + bet * p

        v = calc_ax(p, a_e, a_w, a_s, a_n, a_p)
        nrm = dot_prod(cr, v)

        alf = rho / nrm
        s = r - alf * v

        cs = precon(s, m_e, m_w, m_s, m_n, m_p)
        t = calc_ax(cs, a_e, a_w, a_s, a_n, a_p)

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
    """Calculation of the meridional heat transport using the gradient of the
    energy flux potential.

    Equation (11) in Pearce and Bodas-Salcedo (2023)
    """
    y = np.arange(-0.5 * np.pi + 0.5 * dy, 0.5 * np.pi, dy)

    grad_phi = np.gradient(sol, dy, axis=0)
    grad_phi = grad_phi[1:-1, 1:-1]
    mht = np.sum((grad_phi.T * np.cos(y) * dx).T, axis=1)

    return mht


def spherical_poisson(logger, forcing, tolerance=None):
    """Wrapper function that solves the Poisson equation for a given source
    term (forcing)."""
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
    a_e, a_w, a_s, a_n, a_p, m_e, m_w, m_s, m_n, m_p = set_matrix(hp, hv)

    # Solving the Poisson equation
    rhs = np.zeros([M + 2, N + 2])
    sol = np.zeros([M + 2, N + 2])
    logger.info("spherical_poisson: Solving poisson equation")
    rhs[1:-1, 1:-1] = forcing
    rhs = swap_bounds(rhs)
    sol = bicgstab(logger, sol, rhs, a_e, a_w, a_s, a_n, a_p, m_e, m_w, m_s,
                   m_n, m_p)

    # Calculating meridional heat transport
    logger.info("spherical_poisson: Calculating MHT")
    mht = calc_mht(sol)

    return sol, mht
