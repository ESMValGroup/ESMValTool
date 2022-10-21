"""Poisson solver for over the ocean only. The Poisson equation is solved
numerically using the biconjugate gradient stabilized (BiCGSTAB) method. To
ensure a solution over the ocean only, Neumann boundary conditions are
implemented such that the normal derivative at the boundary is zero. The
simplest method of applying this boundary condition is to zero the land flux.

The solution is achieved when the difference between the net surface flux and
the Laplacian of the output field is less than the stated tolerance. If the
solver fails to converge, the tolerance can be increased.

Convergence is achieved faster by using a preconditioner on the output field.

The meridional heat transport is estimated as the gradient of the scalar
p-field output of the Poisson solver. This is done using a central difference
method with exceptions at the land-sea boundaries.
"""

import numpy as np
from scipy.ndimage import binary_fill_holes


def set_metrics():
    # Defines variables used in the solver
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


def set_matrix(hp, hv, mask):
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

        if mask[j + 1, i + 1] == 0:
            A_w[j, i] = 0.0
            A_e[j, i] = 0.0
            A_s[j, i] = 0.0
            A_n[j, i] = 0.0
            A_p[j, i] = 1.0

    # ILU/SIP preconditioner factors: alf = 0.0 is ILU
    alf = 0.9
    M_p += 1.0

    # Preconditioner only needs to act at sea points
    for j in range(1, M + 1):
        for i in range(1, N + 1):
            if (mask[j, i] == 1):
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
    # As the array is periodic is allows for easier comparison at boundaries
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


def calc_Ax(x, A_e, A_w, A_s, A_n, mask):
    # Laplcian equation
    # Only solved at ocean points (mask = 1)
    # Uses mask to remove contribution of land points at boundary
    Ax = np.zeros([M + 2, N + 2])

    x = swap_bounds(x)
    for j in range(1, M + 1):
        for i in range(1, N + 1):
            if mask[j, i] == 1:
                Ax[j, i] = mask[j-1, i]*A_s[j-1, i-1]*(x[j-1, i]-x[j, i]) + \
                    mask[j, i-1]*A_w[j-1, i-1]*(x[j, i-1]-x[j, i]) + \
                    mask[j, i+1]*A_e[j-1, i-1]*(x[j, i+1]-x[j, i]) + \
                    mask[j+1, i]*A_n[j-1, i-1]*(x[j+1, i]-x[j, i])

    Ax = swap_bounds(Ax)
    return Ax


def dot_prod(x, y):
    # Calculate dot product of two matrices
    dot_prod = 0

    for j in range(1, M + 1):
        for i in range(1, N + 1):
            dot_prod += x[j, i] * y[j, i]

    return dot_prod


def precon(x, M_e, M_w, M_s, M_n, M_p):
    # Preconditioner
    Cx = np.zeros([M + 2, N + 2])
    for j in range(1, M + 1):
        for i in range(1, N + 1):
            Cx[j, i] = M_p[j, i] * (x[j, i] - M_s[j, i] * Cx[j - 1, i] -
                                    M_w[j, i] * Cx[j, i - 1])

    Cx = swap_bounds(Cx)

    for j in range(M, 0, -1):
        for i in range(N, 0, -1):
            Cx[j,
               i] = Cx[j,
                       i] - M_e[j, i] * Cx[j, i + 1] - M_n[j, i] * Cx[j + 1, i]

    Cx = swap_bounds(Cx)
    return Cx


def bicgstab(logger, x, b, A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s, M_n, M_p,
             mask):
    # Bi-conjugate gradient stabilized numerical solver
    sc_err = dot_prod(b, b)

    Ax = calc_Ax(x, A_e, A_w, A_s, A_n, mask)

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

        v = calc_Ax(p, A_e, A_w, A_s, A_n, mask)
        nrm = dot_prod(cr, v)

        alf = rho / nrm
        s = r - alf * v

        cs = precon(s, M_e, M_w, M_s, M_n, M_p)
        t = calc_Ax(cs, A_e, A_w, A_s, A_n, mask)

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


def calc_mht(sol, mask, wrap_mask):
    # MHT = Grad(P)
    grad_phi = np.zeros([M, N])
    wrap_mask = wrap_mask[:, 1:-1]
    sol = sol[:, 1:-1]

    y = np.arange(-0.5 * np.pi + 0.5 * dy, 0.5 * np.pi, dy)

    for i in range(1, M + 1):
        for j in range(0, N):
            # Ocean points only
            if wrap_mask[i, j] == 1:
                if (wrap_mask[i + 1, j] == 0) and (wrap_mask[i - 1, j] == 0):
                    # If land N+S, gradient = 0
                    grad_phi[i - 1, j] = 0
                elif wrap_mask[i + 1, j] == 0:
                    # If land N, assume N = current ocean
                    grad_phi[i - 1, j] = (sol[i, j] - sol[i - 1, j]) / (2 * dy)
                elif wrap_mask[i - 1, j] == 0:
                    # If land S, assume S = current ocean
                    grad_phi[i - 1, j] = (sol[i + 1, j] - sol[i, j]) / (2 * dy)
                else:
                    # Surrounded by ocean
                    grad_phi[i - 1,
                             j] = (sol[i + 1, j] - sol[i - 1, j]) / (2 * dy)

    # Multiply by mask to get individual basin MHT
    grad_phi *= mask
    mht = np.sum((grad_phi.T * np.cos(y) * dx).T, axis=1)

    return mht


def calc_gradient(sol, wrap_mask):
    # MHT = Grad(P)
    # Used to calculate gradient in both x and y directions
    # Needed for quiver plot arrows showing transport direction/magnitude

    # Northwards transport
    grad_y = np.zeros([M, N])
    wrap_y = wrap_mask[:, 1:-1]
    sol_y = sol[:, 1:-1]

    for i in range(1, M + 1):
        for j in range(0, N):
            # Ocean points only
            if wrap_y[i, j] == 1:
                if (wrap_y[i + 1, j] == 0) and (wrap_y[i - 1, j] == 0):
                    grad_y[i - 1, j] = 0
                elif wrap_y[i + 1, j] == 0:
                    grad_y[i - 1,
                           j] = (sol_y[i, j] - sol_y[i - 1, j]) / (2 * dy)
                elif wrap_y[i - 1, j] == 0:
                    grad_y[i - 1,
                           j] = (sol_y[i + 1, j] - sol_y[i, j]) / (2 * dy)
                else:
                    grad_y[i - 1,
                           j] = (sol_y[i + 1, j] - sol_y[i - 1, j]) / (2 * dy)

    del wrap_y, sol_y

    # Eastwards transport
    grad_x = np.zeros([M, N])
    wrap_x = wrap_mask[1:-1, :]
    sol_x = sol[1:-1, :]

    for i in range(0, M):
        for j in range(1, N + 1):
            # Ocean points only
            if wrap_x[i, j] == 1:
                if (wrap_x[i, j + 1] == 0) and (wrap_x[i, j - 1] == 0):
                    grad_x[i, j - 1] = 0
                elif wrap_x[i, j + 1] == 0:
                    grad_x[i,
                           j - 1] = (sol_x[i, j] - sol_x[i, j - 1]) / (2 * dx)
                elif wrap_x[i, j - 1] == 0:
                    grad_x[i,
                           j - 1] = (sol_x[i, j + 1] - sol_x[i, j]) / (2 * dx)
                else:
                    grad_x[i, j -
                           1] = (sol_x[i, j + 1] - sol_x[i, j - 1]) / (2 * dx)

    del wrap_x, sol_x

    return grad_x, grad_y


def define_masks(mask):
    # BASINS
    # Drake Passage: 71 W = 360-71 = 289 E
    # Southern Ocean: 147 E
    # Southern Ocean: 20 E
    # Arctic: 80 W = 360-80 = 280 E
    # Artic: 100 E = 100

    xd = 360 / N
    yd = 180 / M

    # Pacific
    pacific = mask.copy()
    pacific = 1 - pacific
    pacific[:, int(np.ceil(289 / xd)):] = 1  # Drake Passage
    pacific[:int(np.ceil(81 / yd)), :int(np.floor(147 /
                                                  xd))] = 1  # Southern Ocean
    pacific[int(np.floor(97 / yd)):, int(np.floor(280 / xd)):] = 1  # Arctic
    pacific[int(np.floor(138 / yd)):, :int(np.floor(54 / xd))] = 1  # Arctic
    pacific[int(np.floor(75 / yd)):int(np.ceil(100 /
                                               yd)), :int(np.ceil(38 /
                                                                  xd))] = 1
    pacific = binary_fill_holes(pacific)
    # Binary fill holes used to remove inland seas
    # Also removes e.g. Med if the resolution is too low so the sea is cut off
    # from the ocean

    # Indian
    indian = mask.copy()
    indian = 1 - indian + 1 - pacific
    indian[:, :int(np.floor(20 / xd))] = 1
    indian[:, int(np.ceil(200 / xd)):] = 1
    indian[int(np.floor(150 / yd)):, :] = 1

    # Atlantic
    atlantic = mask.copy()
    atlantic = 1 - atlantic + 1 - pacific + 1 - indian

    masks = {
        'global': mask,
        'atlantic': 1 - atlantic,
        'pacific': 1 - pacific,
        'indian': 1 - indian
    }

    return masks


def spherical_ocean_poisson(logger, forcing, mask, tolerance=None):
    # Define global variables
    logger.info("spherical_ocean_poisson: Setting global variables")

    global N, M, dx, dy, tol
    (M, N) = (np.shape(mask))
    N = int(N)
    M = int(M)
    dx = 2.0 * np.pi / N
    dy = np.pi / M
    if tolerance is None:
        tol = 3.0e-2
    else:
        tol = tolerance

    logger.info("spherical_ocean_poisson: Calling set_metrics")
    hp, hv = set_metrics()

    # Defining masks
    wrap_mask = np.zeros([M + 2, N + 2])
    wrap_mask[1:-1, 1:-1] = mask
    wrap_mask = swap_bounds(wrap_mask)

    basins = define_masks(mask)

    # Solving the Poisson equation
    logger.info("spherical_ocean_poisson: Calling set_matrix")
    A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s, M_n, M_p = set_matrix(
        hp, hv, wrap_mask)

    rhs = np.zeros([M + 2, N + 2])
    rhs[1:-1, 1:-1] = forcing
    rhs = swap_bounds(rhs)

    logger.info("spherical_ocean_poisson: Solving poisson equation")
    sol = np.zeros([M + 2, N + 2])
    sol = bicgstab(logger, sol, rhs, A_e, A_w, A_s, A_n, A_p, M_e, M_w, M_s,
                   M_n, M_p, wrap_mask)

    # Calculating meridional heat transport
    logger.info("spherical_ocean_poisson: Calculating MHT")  #
    mht = {}
    for x in basins:
        mht[x] = calc_mht(sol, basins[x], wrap_mask)

    x, y = calc_gradient(sol, wrap_mask)

    return sol, mht, x, y
