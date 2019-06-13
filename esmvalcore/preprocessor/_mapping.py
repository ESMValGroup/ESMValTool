# -*- coding: utf-8 -*-
"""Provides mapping of a cube."""

import collections

import iris
import numpy as np
import six


def _is_single_item(testee):
    """
    Check if testee is a single item.

    Return whether this is a single item, rather than an iterable.
    We count string types as 'single', also.
    """
    return (isinstance(testee, six.string_types)
            or not isinstance(testee, collections.Iterable))


def _as_list_of_coords(cube, names_or_coords):
    """Convert a name, coord, or list of names/coords to a list of coords."""
    # If not iterable, convert to list of a single item
    if _is_single_item(names_or_coords):
        names_or_coords = [names_or_coords]
    coords = []
    for name_or_coord in names_or_coords:
        if isinstance(name_or_coord, (iris.coords.Coord, six.string_types)):
            coords.append(cube.coord(name_or_coord))
        else:
            # Don't know how to handle this type
            msg = ("Don't know how to handle coordinate of type %s. "
                   "Ensure all coordinates are of type six.string_types "
                   "or iris.coords.Coord.") % (type(name_or_coord), )
            raise TypeError(msg)
    return coords


def ref_to_dims_index_as_coordinate(cube, ref):
    """Get dims for coord ref."""
    coord = _as_list_of_coords(cube, ref)[0]
    dims = cube.coord_dims(coord)
    if not dims:
        msg = ('Requested an iterator over a coordinate ({}) '
               'which does not describe a dimension.')
        msg = msg.format(coord.name())
        raise ValueError(msg)
    return dims


def ref_to_dims_index_as_index(cube, ref):
    """Get dim for index ref."""
    try:
        dim = int(ref)
    except (ValueError, TypeError):
        raise ValueError('{} Incompatible type {} for '
                         'slicing'.format(ref, type(ref)))
    if dim < 0 or dim > cube.ndim:
        msg = ('Requested an iterator over a dimension ({}) '
               'which does not exist.'.format(dim))
        raise ValueError(msg)
    dims = [dim]
    return dims


def ref_to_dims_index(cube, ref_to_slice):
    """
    Map a list of :class:`iris.coords.DimCoord` to a tuple of indices.

    This method finds the indices of the dimensions in a cube that collectively
    correspond to the given list of :class:`iris.coords.DimCoord`.

    Parameters
    ----------
    cube: :class:`iris.cube.Cube`
        The cube to examine.
    ref_to_slice: iterable of or single :class:`iris.coords.DimCoord`
        Specification of the dimensions in terms of coordinates.

    Returns
    -------
    tuple:
        A tuple of indices corresponding to the given dimensions.
    """
    # Required to handle a mix between types
    if _is_single_item(ref_to_slice):
        ref_to_slice = [ref_to_slice]
    dim_to_slice = []
    dim_to_slice_set = set()
    for ref in ref_to_slice:
        try:
            dims = ref_to_dims_index_as_coordinate(cube, ref)
        except TypeError:
            dims = ref_to_dims_index_as_index(cube, ref)
        for dim in dims:
            if dim not in dim_to_slice_set:
                dim_to_slice.append(dim)
                dim_to_slice_set.add(dim)
    return dim_to_slice


def get_associated_coords(cube, dimensions):
    """
    Return all coords containing any of the given dimensions.

    Return all coords, dimensional and auxiliary, that contain any of
    the given dimensions.
    """
    dims = []
    dim_set = set()
    for dim in dimensions:
        if dim not in dim_set:
            dims.append(dim)
            dim_set.add(dim)
    dim_coords = []
    for i in dims:
        coords = cube.coords(contains_dimension=i, dim_coords=True)
        if coords:
            dim_coords.append(coords[0])
    aux_coords = []
    for i in dims:
        coords = cube.coords(contains_dimension=i, dim_coords=False)
        if coords:
            aux_coords.append(coords[0])
    return dim_coords, aux_coords


def get_empty_data(shape, dtype=np.float32):
    """
    Create an empty data object of the given shape.

    Creates an emtpy data object of the given shape, potentially of the lazy
    kind from biggus or dask, depending on the used iris version.
    """
    data = np.empty(shape, dtype=dtype)
    mask = np.empty(shape, dtype=bool)
    return np.ma.masked_array(data, mask)


def get_slice_spec(cube, ref_to_slice):
    """
    Turn a slice reference into a specification for the slice.

    Turns a slice reference into a specification comprised of the shape as well
    as the relevant dimensional and auxiliary coordinates.
    """
    slice_dims = ref_to_dims_index(cube, ref_to_slice)
    slice_shape = tuple(cube.shape[d] for d in slice_dims)
    dim_coords, aux_coords = get_associated_coords(cube, slice_dims)
    return slice_shape, dim_coords, aux_coords


def index_iterator(dims_to_slice, shape):
    """
    Return iterator for subsets of multidimensional objects.

    An iterator over a multidimensional object, giving both source and
    destination indices.
    """
    dst_slices = (slice(None, None),) * len(dims_to_slice)
    dims = [1 if n in dims_to_slice else i for n, i in enumerate(shape)]
    for index_tuple in np.ndindex(*dims):
        src_ind = tuple(
            slice(None, None) if n in dims_to_slice else i
            for n, i in enumerate(index_tuple))
        dst_ind = tuple(i for n, i in enumerate(index_tuple)
                        if n not in dims_to_slice) + dst_slices
        yield src_ind, dst_ind


def get_slice_coords(cube):
    """Return ordered set of unique coordinates."""
    slice_coords = []
    slice_set = set()
    for i in range(cube.ndim):
        coords = cube.coords(contains_dimension=i)
        for coord in coords:
            if coord not in slice_set:
                slice_coords.append(coord)
                slice_set.add(coord)
    return slice_coords


def map_slices(src, func, src_rep, dst_rep):
    """
    Map slices of a cube, replacing them with different slices.

    This method is similar to the standard cube collapsed and aggregated_by
    methods, however, where they completely remove the mapped dimensions, this
    method allows for their replacement with other dimensions.
    The new dimensions are specified with a destination representant and will
    be the last dimensions of the resulting cube, even if the removed
    dimensions are can be any of the source cubes dimensions.

    Parameters
    ----------
    src: :class:`iris.cube.Cube`
        Source cube to be mapped.
    func: callable
        Callable that takes a single cube and returns a single numpy array.
    src_rep: :class:`iris.cube.Cube`
        Source representant that specifies the dimensions to be removed from
        the source cube.
    dst_rep: :class:`iris.cube.Cube`
        Destination representant that specifies the shape of the new
        dimensions.

    Returns
    -------
    :class:`iris.cube.Cube`:
        New cube that has the shape of the source cube with the removed
        dimensions replaced with the destination dimensions.
        All coordinates that span any of the removed dimensions are removed;
        :class:`iris.coords.DimCoord` for the new dimensions are taken from
        `dst_rep`.
    """
    ref_to_slice = get_slice_coords(src_rep)
    src_slice_dims = ref_to_dims_index(src, ref_to_slice)
    src_keep_dims = list(set(range(src.ndim)) - set(src_slice_dims))
    src_keep_spec = get_slice_spec(src, src_keep_dims)
    res_shape = src_keep_spec[0] + dst_rep.shape
    dim_coords = src_keep_spec[1] + dst_rep.coords(dim_coords=True)
    dim_coords_and_dims = [(c, i) for i, c in enumerate(dim_coords)]
    dst = iris.cube.Cube(
        data=get_empty_data(res_shape, dtype=src.dtype),
        standard_name=src.standard_name,
        long_name=src.long_name,
        var_name=src.var_name,
        units=src.units,
        attributes=src.attributes,
        cell_methods=src.cell_methods,
        dim_coords_and_dims=dim_coords_and_dims,
    )
    for src_ind, dst_ind in index_iterator(src_slice_dims, src.shape):
        res = func(src[src_ind])
        dst.data[dst_ind] = res
    return dst
