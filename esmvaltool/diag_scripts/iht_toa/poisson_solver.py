# (C) Crown Copyright 2023, the Met Office.
"""Poisson solver for the full ocean-atmosphere column.

The Poisson equation is solved by numerically using the bi-conjugate
gradient stabilized (BiCGSTAB) method.

The solution is achieved when the difference between the input field (radiative
flux) and the Laplacian of the output field is less than the stated tolerance.
If the solver fails to converge, the tolerance can be increased.

Convergence is achieved faster by using a preconditioner on the output field.

The heat transport is calculated as the gradient of the energy flux potential,
the output of the Poisson solver.
"""

import numpy as np
from numba import jit


def swap_bounds(array):
    """Extend the array by one in all directions.

    As the array is periodic it allows for easier computations at
    boundaries.
    """
    shape0, shape1 = np.array(array.shape) - 2
    wrap_point = int(shape1 / 2 + 1)
    for i in range(1, shape1 + 1):
        array[0, i] = array[1, wrap_point]
        array[shape0 + 1, i] = array[shape0, wrap_point]
        wrap_point += 1
        if wrap_point > shape1:
            wrap_point = 1

    array[:, 0] = array[:, shape1]
    array[:, shape1 + 1] = array[:, 1]

    return array


def dot_prod(a_matrix, b_matrix):
    """Calculate dot product of two matrices only over source term size."""
    shape0, shape1 = np.array(a_matrix.shape) - 2
    return (
        a_matrix[1 : shape0 + 1, 1 : shape1 + 1]
        * b_matrix[1 : shape0 + 1, 1 : shape1 + 1]
    ).sum()


def precon(x_matrix, m_matrix):
    """Preconditioner.

    This is a wrapper to two steps that are optimised using jit.
    It implements the preconditioning step of van der Vorst, H. A., 1992.
    https://doi.org/10.1137/0913035.
    """
    cx_matrix = np.zeros(np.array(x_matrix.shape))
    precon_a(x_matrix, m_matrix[1], m_matrix[2], m_matrix[4], cx_matrix)
    cx_matrix = swap_bounds(cx_matrix)
    precon_b(m_matrix[0], m_matrix[3], cx_matrix)
    cx_matrix = swap_bounds(cx_matrix)
    return cx_matrix


@jit(nopython=True)
def precon_a(x_matrix, m_w, m_s, m_p, cx_matrix):
    """First step of preconditioner."""
    shape0, shape1 = np.array(cx_matrix.shape) - 2
    for j in range(1, shape0 + 1):
        for i in range(1, shape1 + 1):
            cx_matrix[j, i] = m_p[j, i] * (
                x_matrix[j, i]
                - m_s[j, i] * cx_matrix[j - 1, i]
                - m_w[j, i] * cx_matrix[j, i - 1]
            )


@jit(nopython=True)
def precon_b(m_e, m_n, cx_matrix):
    """Second step of preconditioner."""
    shape0, shape1 = np.array(cx_matrix.shape) - 2
    for j in range(shape0, 0, -1):
        for i in range(shape1, 0, -1):
            cx_matrix[j, i] = (
                cx_matrix[j, i]
                - m_e[j, i] * cx_matrix[j, i + 1]
                - m_n[j, i] * cx_matrix[j + 1, i]
            )


class SphericalPoisson:
    """Poisson solver over the sphere.

    Solve Poisson equation for a given source term (forcing) and
    calculate meridional heat transport (MHT).
    """

    def __init__(self, logger, source, tolerance=2.0e-4):
        """Initialise solver with source field, metrics and matrices."""
        self.logger = logger
        self.source = source
        self.tolerance = tolerance
        self.energy_flux_potential = None
        self.meridional_heat_transport = None
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
        aaa = 1.0 / ((2.0 * np.pi / src_shape[1]) ** 2.0)
        bbb = 1.0 / ((np.pi / src_shape[0]) ** 2.0)

        # First calculate the Poisson equations 5-point stencil
        # A_w is the contribution from i-1, A_e is from i+1,
        # A_s is j-1, A_n is j+1, and A_p is the diagonal
        for j in range(0, src_shape[0]):
            txa = aaa / hpi[j] ** 2.0
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
                m_matrix[2, j, i] = a_matrix[2, j - 1, i - 1] / (
                    1.0 + alf * m_matrix[0, j - 1, i]
                )

                m_matrix[1, j, i] = a_matrix[1, j - 1, i - 1] / (
                    1.0 + alf * m_matrix[3, j, i - 1]
                )

                m_matrix[4, j, i] = (
                    a_matrix[4, j - 1, i - 1]
                    - m_matrix[2, j, i]
                    * (m_matrix[3, j - 1, i] - alf * m_matrix[0, j - 1, i])
                    - m_matrix[1, j, i]
                    * (m_matrix[0, j, i - 1] - alf * m_matrix[3, j, i - 1])
                )

                m_matrix[4, j, i] = 1.0 / m_matrix[4, j, i]

                m_matrix[0, j, i] = (
                    a_matrix[0, j - 1, i - 1]
                    - alf * m_matrix[2, j, i] * m_matrix[0, j - 1, i]
                ) * m_matrix[4, j, i]

                m_matrix[3, j, i] = (
                    a_matrix[3, j - 1, i - 1]
                    - alf * m_matrix[1, j, i] * m_matrix[3, j, i - 1]
                ) * m_matrix[4, j, i]

        self.a_matrix = a_matrix
        self.m_matrix = m_matrix

    def solve(self, max_iterations=1000):
        """Solve equation for the source term.

        Bi-conjugate gradient stabilized numerical solver: van der
        Vorst, H. A., 1992: Bi-cgstab: A fast and smoothly converging
        variant of bi-cg for the solution of nonsymmetric linear
        systems. SIAM Journal on Scientific and Statistical Computing,
        https://doi.org/10.1137/0913035.
        This solver implements the preconditioned Bi-CGSTAB algorithm,
        described in page 638 of that paper.
        """
        bbb = np.zeros(np.array(self.source.shape) + 2)
        xxx = np.zeros(np.array(self.source.shape) + 2)
        bbb[1:-1, 1:-1] = self.source
        bbb = swap_bounds(bbb)

        sc_err = dot_prod(bbb, bbb)

        # Group some temporal variables
        stv = {
            "alf": 1.0,
            "omg": 1.0,
            "nrm": 1.0,
            "rrr": bbb - self.calc_ax(xxx),
        }
        stv["crrr"] = stv["rrr"].copy()

        ppp = np.zeros(np.array(self.source.shape) + 2)
        vvv = np.zeros(np.array(self.source.shape) + 2)

        iteration = 0
        while iteration < max_iterations:
            rho = dot_prod(stv["rrr"], stv["crrr"])

            bet = (rho / stv["nrm"]) * (stv["alf"] / stv["omg"])

            ttt = stv["rrr"] - bet * stv["omg"] * vvv

            sss = precon(ttt, self.m_matrix)
            ppp = sss + bet * ppp

            vvv = self.calc_ax(ppp)
            stv["nrm"] = dot_prod(stv["crrr"], vvv)

            stv["alf"] = rho / stv["nrm"]
            sss = stv["rrr"] - stv["alf"] * vvv

            csss = precon(sss, self.m_matrix)
            ttt = self.calc_ax(csss)

            stv["omg"] = dot_prod(ttt, sss) / dot_prod(ttt, ttt)

            xxx = xxx + stv["alf"] * ppp + stv["omg"] * csss
            stv["rrr"] = sss - stv["omg"] * ttt

            stv["nrm"] = rho

            if abs(stv["omg"]) < 1.0e-16:
                self.logger.info("Terminating Poisson solver.")
                break

            err = np.sqrt(dot_prod(stv["rrr"], stv["rrr"]) / sc_err)
            if err < self.tolerance:
                self.logger.info("Poisson solver has converged.")
                break

            iteration += 1

            if iteration == max_iterations:
                raise RuntimeError("Poisson solver has not converged.")

        self.energy_flux_potential = xxx

    def calc_meridional_heat_transport(self):
        """Meridional heat transport of energy flux potential.

        Calculate of the meridional heat transport using the gradient of
        the energy flux potential. Equation (11) in Pearce and Bodas-
        Salcedo (2023).
        """
        deltax = 2.0 * np.pi / self.source.shape[1]
        deltay = np.pi / self.source.shape[0]
        yvalues = np.arange(-0.5 * np.pi + 0.5 * deltay, 0.5 * np.pi, deltay)
        grad_phi = np.gradient(self.energy_flux_potential, deltay, axis=0)
        grad_phi = grad_phi[1:-1, 1:-1]
        self.meridional_heat_transport = np.sum(
            (grad_phi.T * np.cos(yvalues) * deltax).T, axis=1
        )

    def calc_ax(self, x_matrix):
        """Matrix calculation of the Laplacian equation, LHS of Eq.

        (9) in Pearce and Bodas-Salcedo (2023).
        """
        # Laplacian equation
        src_shape = np.array(self.source.shape)
        ax_matrix = np.zeros(src_shape + 2)
        x_matrix = swap_bounds(x_matrix)
        shape0, shape1 = src_shape
        ax_matrix[1 : shape0 + 1, 1 : shape1 + 1] = (
            self.a_matrix[2, 0:shape0, 0:shape1]
            * x_matrix[0:shape0, 1 : shape1 + 1]
            + self.a_matrix[1, 0:shape0, 0:shape1]
            * x_matrix[1 : shape0 + 1, 0:shape1]
            + self.a_matrix[0, 0:shape0, 0:shape1]
            * x_matrix[1 : shape0 + 1, 2 : shape1 + 2]
            + self.a_matrix[3, 0:shape0, 0:shape1]
            * x_matrix[2 : shape0 + 2, 1 : shape1 + 1]
            + self.a_matrix[4, 0:shape0, 0:shape1]
            * x_matrix[1 : shape0 + 1, 1 : shape1 + 1]
        )
        ax_matrix = swap_bounds(ax_matrix)
        return ax_matrix
