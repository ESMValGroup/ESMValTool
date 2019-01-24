"""
rdp
from rdp 0.8 (not installable with conda)
~~~

Python implementation of the Ramer-Douglas-Peucker algorithm.

:copyright: 2014-2016 Fabian Hirschmann <fabian@hirschmann.email>
:license: MIT, see LICENSE_rdp.txt for more details.

"""
import sys
from functools import partial
import numpy as np

if sys.version_info[0] >= 3:
    xrange = range


def pldist(point, start, end):
    """
    Calculates the distance from ``point`` to the line given
    by the points ``start`` and ``end``.

    :param point: a point
    :type point: numpy array
    :param start: a point of the line
    :type start: numpy array
    :param end: another point of the line
    :type end: numpy array
    """
    if np.all(np.equal(start, end)):
        return np.linalg.norm(point - start)

    return np.divide(
        np.abs(np.linalg.norm(np.cross(end - start, start - point))),
        np.linalg.norm(end - start))


def rdp_rec(M, epsilon, dist=pldist):
    """
    Simplifies a given array of points.

    Recursive version.

    :param M: an array
    :type M: numpy array
    :param epsilon: epsilon in the rdp algorithm
    :type epsilon: float
    :param dist: distance function
    :type dist: function with signature ``f(point, start, end)``
      -- see :func:`rdp.pldist`
    """
    dmax = 0.0
    index = -1

    for i in xrange(1, M.shape[0]):
        d = dist(M[i], M[0], M[-1])

        if d > dmax:
            index = i
            dmax = d

    if dmax > epsilon:
        r1 = rdp_rec(M[:index + 1], epsilon, dist)
        r2 = rdp_rec(M[index:], epsilon, dist)

        return np.vstack((r1[:-1], r2))
    else:
        return np.vstack((M[0], M[-1]))


def _rdp_iter(M, start_index, last_index, epsilon, dist=pldist):
    stk = []
    stk.append([start_index, last_index])
    global_start_index = start_index
    indices = np.ones(last_index - start_index + 1, dtype=bool)

    while stk:
        start_index, last_index = stk.pop()

        dmax = 0.0
        index = start_index

        for i in xrange(index + 1, last_index):
            if indices[i - global_start_index]:
                d = dist(M[i], M[start_index], M[last_index])
                if d > dmax:
                    index = i
                    dmax = d

        if dmax > epsilon:
            stk.append([start_index, index])
            stk.append([index, last_index])
        else:
            for i in xrange(start_index + 1, last_index):
                indices[i - global_start_index] = False

    return indices


def rdp_iter(M, epsilon, dist=pldist, return_mask=False):
    """
    Simplifies a given array of points.

    Iterative version.

    :param M: an array
    :type M: numpy array
    :param epsilon: epsilon in the rdp algorithm
    :type epsilon: float
    :param dist: distance function
    :type dist: function with signature ``f(point, start, end)``
     -- see :func:`rdp.pldist`
    :param return_mask: return the mask of points to keep instead
    :type return_mask: bool
    """
    mask = _rdp_iter(M, 0, len(M) - 1, epsilon, dist)

    if return_mask:
        return mask

    return M[mask]


def rdp(M, epsilon=0, dist=pldist, algo="iter", return_mask=False):
    """
    Simplifies a given array of points using the Ramer-Douglas-Peucker
    algorithm.

    Example:

    >>> from rdp import rdp
    >>> rdp([[1, 1], [2, 2], [3, 3], [4, 4]])
    [[1, 1], [4, 4]]

    This is a convenience wrapper around both :func:`rdp.rdp_iter`
    and :func:`rdp.rdp_rec` that detects if the input is a numpy array
    in order to adapt the output accordingly. This means that
    when it is called using a Python list as argument, a Python
    list is returned, and in case of an invocation using a numpy
    array, a NumPy array is returned.

    The parameter ``return_mask=True`` can be used in conjunction
    with ``algo="iter"`` to return only the mask of points to keep. Example:

    >>> from rdp import rdp
    >>> import numpy as np
    >>> arr = np.array([1, 1, 2, 2, 3, 3, 4, 4]).reshape(4, 2)
    >>> arr
    array([[1, 1],
           [2, 2],
           [3, 3],
           [4, 4]])
    >>> mask = rdp(arr, algo="iter", return_mask=True)
    >>> mask
    array([ True, False, False,  True], dtype=bool)
    >>> arr[mask]
    array([[1, 1],
           [4, 4]])

    :param M: a series of points
    :type M: numpy array with shape ``(n,d)`` where ``n`` is the number of
     points and ``d`` their dimension
    :param epsilon: epsilon in the rdp algorithm
    :type epsilon: float
    :param dist: distance function
    :type dist: function with signature ``f(point, start, end)``
     -- see :func:`rdp.pldist`
    :param algo: either ``iter`` for an iterative algorithm or ``rec``
     for a recursive algorithm
    :type algo: string
    :param return_mask: return mask instead of simplified array
    :type return_mask: bool
    """

    if algo == "iter":
        algo = partial(rdp_iter, return_mask=return_mask)
    elif algo == "rec":
        if return_mask:
            raise NotImplementedError("return_mask=True not " +
                                      "supported with algo=\"rec\"")
        algo = rdp_rec

    if "numpy" in str(type(M)):
        return algo(M, epsilon, dist)

    return algo(np.array(M), epsilon, dist).tolist()
