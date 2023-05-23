"""(C) Crown Copyright 2023, the Met Office.

Poisson solver for the full ocean-atmosphere column. The Poisson equation
is solved by numerically using the biconjugate gradient stabilized (BiCGSTAB)
method.

The solution is achieved when the difference between the input field (radiative
flux) and the Laplacian of the output field is less than the stated tolerance.
If the solver fails to converge, the tolerance can be increased.

Convergence is achieved faster by using a preconditioner on the output field.

The heat transport is calculated as the gradient of the energy flux potential,
the output of the Poisson solver.
"""

import sys
import numpy as np
from numba import jit


def swap_bounds(fld):
    """Extend the array by one in all directions.

    As the array is periodic it allows for easier computations at
    boundaries.
    """
    shp0, shp1 = np.array(fld.shape) - 2
    wrap_pnt = int(shp1 / 2 + 1)
    for i in range(1, shp1 + 1):
        fld[0, i] = fld[1, wrap_pnt]
        fld[shp0 + 1, i] = fld[shp0, wrap_pnt]
        wrap_pnt += 1
        if wrap_pnt > shp1:
            wrap_pnt = 1

    fld[:, 0] = fld[:, shp1]
    fld[:, shp1 + 1] = fld[:, 1]

    return fld


def dot_prod(xxx, yyy):
    """Calculate dot product of two matrices only over source term size."""
    shp0, shp1 = np.array(xxx.shape) - 2
    return (xxx[1:shp0 + 1, 1:shp1 + 1] * yyy[1:shp0 + 1, 1:shp1 + 1]).sum()


@jit
def precon(xxx, m_matrix):
    """Preconditioner.

    This is a wrapper to two steps that are optimised using jit.
    """
    cxxx = np.zeros(np.array(xxx.shape))
    precon_a(xxx, m_matrix[1], m_matrix[2], m_matrix[4], cxxx)
    cxxx = swap_bounds(cxxx)
    precon_b(m_matrix[0], m_matrix[3], cxxx)
    cxxx = swap_bounds(cxxx)
    return cxxx


@jit
def precon_a(xxx, m_w, m_s, m_p, cxxx):
    """First step of preconditioner."""
    shp0, shp1 = np.array(cxxx.shape) - 2
    for j in range(1, shp0 + 1):
        for i in range(1, shp1 + 1):
            cxxx[j, i] = m_p[j, i] * \
                (xxx[j, i] - m_s[j, i] * cxxx[j - 1, i] -
                 m_w[j, i] * cxxx[j, i - 1])


@jit
def precon_b(m_e, m_n, cxxx):
    """Second step of preconditioner."""
    shp0, shp1 = np.array(cxxx.shape) - 2
    for j in range(shp0, 0, -1):
        for i in range(shp1, 0, -1):
            cxxx[j, i] = cxxx[j, i] - m_e[j, i] * cxxx[j, i + 1] - \
                m_n[j, i] * cxxx[j + 1, i]


class SphericalPoisson:
    """Poisson solver over the sphere.

    Solve Poisson equation for a given source term (forcing) and
    calculate MHT.
    """

    def __init__(self, logger, source, tolerance=2.0e-4):
        """Initialise solver with source field, metrics and matrices."""
        self.logger = logger
        self.source = source
        self.tolerance = tolerance
        self.efp = None
        self.mht = None
        logger.info("Initialising Poisson solver.")
        self.set_matrices()

    def set_matrices(self):
        """Calculate A and M matrices.

        A is the matrix that defines the five-point stencil (Eq. 8). The
        A_matrix are the values are the contributions from each of the
        four neighbouring cells: e,w,s,n,p.
        """
        # Calculate metrics hpi and hvj
        src_shape = np.array(self.source.shape)
        hpi = np.zeros(src_shape[0])
        hvj = np.zeros(src_shape[0] + 1)
        deltay = np.pi / src_shape[0]
        yyy = -0.5 * np.pi + 0.5 * deltay
        hvj[0] = 0.0
        for j in range(0, src_shape[0]):
            hpi[j] = np.cos(yyy)
            hvj[j + 1] = np.cos(yyy + 0.5 * deltay)
            yyy += deltay
        hvj[-1] = 0.0

        # Storing the full matrix
        a_matrix = np.zeros((5, *src_shape))

        # ILU factors
        m_matrix = np.zeros((5, *(src_shape + 1)))

        # Spherical Laplacian variables
        aaa = 1.0 / ((2.0 * np.pi / src_shape[1])**2.)
        bbb = 1.0 / ((np.pi / src_shape[0])**2.)

        # First calculate the Poisson equations 5-point stencil
        # A_w is the contribution from i-1, A_e is from i+1,
        # A_s is j-1, A_n is j+1, and A_p is the diagonal
        for j in range(0, src_shape[0]):
            txa = aaa / hpi[j]**2.0
            tyb = bbb / hpi[j]

            for i in range(0, src_shape[1]):
                a_matrix[0, j, i] = txa
                a_matrix[1, j, i] = txa
                a_matrix[2, j, i] = tyb * hvj[j]
                a_matrix[3, j, i] = tyb * hvj[j + 1]
                a_matrix[4, j, i] = -a_matrix[0:4, j, i].sum()

        # ILU/SIP preconditioner factors: alf = 0.0 is ILU
        alf = 0.9
        m_matrix[4] += 1.0

        for j in range(1, src_shape[0] + 1):
            for i in range(1, src_shape[1] + 1):
                m_matrix[2, j, i] = a_matrix[2, j - 1, i - 1] / \
                    (1.0 + alf * m_matrix[0, j - 1, i])

                m_matrix[1, j, i] = a_matrix[1, j - 1, i - 1] / \
                    (1.0 + alf * m_matrix[3, j, i - 1])

                m_matrix[4, j, i] = a_matrix[4, j - 1, i - 1] - \
                    m_matrix[2, j, i] * (m_matrix[3, j - 1, i] -
                                         alf * m_matrix[0, j - 1, i]) - \
                    m_matrix[1, j, i] * (m_matrix[0, j, i - 1] -
                                         alf * m_matrix[3, j, i - 1])

                m_matrix[4, j, i] = 1.0 / m_matrix[4, j, i]

                m_matrix[0, j, i] = (a_matrix[0, j - 1, i - 1] -
                                     alf * m_matrix[2, j, i] *
                                     m_matrix[0, j - 1, i]) * \
                    m_matrix[4, j, i]

                m_matrix[3, j, i] = (a_matrix[3, j - 1, i - 1] -
                                     alf * m_matrix[1, j, i] *
                                     m_matrix[3, j, i - 1]) * \
                    m_matrix[4, j, i]

        self.a_matrix = a_matrix
        self.m_matrix = m_matrix

    def solve(self, max_iterations=1000):
        """Solve equation for the source term.

        Bi-conjugate gradient stabilized numerical solver: van der
        Vorst, H. A., 1992: Bi-cgstab: A fast and smoothly converging
        variant of bi-cg for the solution of nonsymmetric linear
        systems. SIAM Journal on Scientific and Statistical Computing,
        https://doi.org/10.1137/0913035.
        """
        bbb = np.zeros(np.array(self.source.shape) + 2)
        xxx = np.zeros(np.array(self.source.shape) + 2)
        bbb[1:-1, 1:-1] = self.source
        bbb = swap_bounds(bbb)

        sc_err = dot_prod(bbb, bbb)

        # Group some temporal variables
        stv = {
            'alf': 1.0,
            'omg': 1.0,
            'nrm': 1.0,
            'rrr': bbb - self.calc_ax(xxx)
        }
        stv['crrr'] = stv['rrr'].copy()

        ppp = np.zeros(np.array(self.source.shape) + 2)
        vvv = np.zeros(np.array(self.source.shape) + 2)

        iteration = 0
        while iteration < max_iterations:
            rho = dot_prod(stv['rrr'], stv['crrr'])

            bet = (rho / stv['nrm']) * (stv['alf'] / stv['omg'])

            ttt = stv['rrr'] - bet * stv['omg'] * vvv

            sss = precon(ttt, self.m_matrix)
            ppp = sss + bet * ppp

            vvv = self.calc_ax(ppp)
            stv['nrm'] = dot_prod(stv['crrr'], vvv)

            stv['alf'] = rho / stv['nrm']
            sss = stv['rrr'] - stv['alf'] * vvv

            csss = precon(sss, self.m_matrix)
            ttt = self.calc_ax(csss)

            stv['omg'] = dot_prod(ttt, sss) / dot_prod(ttt, ttt)

            xxx = xxx + stv['alf'] * ppp + stv['omg'] * csss
            stv['rrr'] = sss - stv['omg'] * ttt

            stv['nrm'] = rho

            if abs(stv['omg']) < 1.0e-16:
                self.logger.info('Terminating Poisson solver.')
                break

            err = np.sqrt(dot_prod(stv['rrr'], stv['rrr']) / sc_err)
            if err < self.tolerance:
                self.logger.info('Poisson solver has converged.')
                break

            iteration += 1

            if iteration == max_iterations:
                self.logger.error('Poisson solver has not converged.')
                sys.exit(1)

        self.efp = xxx

    def calc_mht(self):
        """Meridional heat transport of energy flux potential.

        Calculate of the meridional heat transport using the gradient of
        the energy flux potential. Equation (11) in Pearce and Bodas-
        Salcedo (2023).
        """
        deltax = 2.0 * np.pi / self.source.shape[1]
        deltay = np.pi / self.source.shape[0]
        yyy = np.arange(-0.5 * np.pi + 0.5 * deltay, 0.5 * np.pi, deltay)
        grad_phi = np.gradient(self.efp, deltay, axis=0)
        grad_phi = grad_phi[1:-1, 1:-1]
        self.mht = np.sum((grad_phi.T * np.cos(yyy) * deltax).T, axis=1)

    def calc_ax(self, xxx):
        """Matrix calculation of the Laplacian equation, LHS of Eq.

        (9) in Pearce and Bodas-Salcedo (2023).
        """
        # Laplacian equation
        src_shape = np.array(self.source.shape)
        axxx = np.zeros(src_shape + 2)
        xxx = swap_bounds(xxx)
        shp0, shp1 = src_shape
        axxx[1:shp0 + 1, 1:shp1 + 1] = \
            self.a_matrix[2, 0:shp0, 0:shp1] * xxx[0:shp0, 1:shp1 + 1] + \
            self.a_matrix[1, 0:shp0, 0:shp1] * xxx[1:shp0 + 1, 0:shp1] + \
            self.a_matrix[0, 0:shp0, 0:shp1] * xxx[1:shp0 + 1, 2:shp1 + 2] + \
            self.a_matrix[3, 0:shp0, 0:shp1] * xxx[2:shp0 + 2, 1:shp1 + 1] + \
            self.a_matrix[4, 0:shp0, 0:shp1] * xxx[1:shp0 + 1, 1:shp1 + 1]
        axxx = swap_bounds(axxx)
        return axxx
