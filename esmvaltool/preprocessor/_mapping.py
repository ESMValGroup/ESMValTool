# -*- coding: utf-8 -*-
"""
Provides mapping of a cube
"""

import collections
import itertools

import numpy as np
import six

import iris


def _is_single_item(testee):
    """
    Checks if testee is a single item

    Return whether this is a single item, rather than an iterable.
    We count string types as 'single', also.
    """
    return (isinstance(testee, six.string_types) or
            not isinstance(testee, collections.Iterable))


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


def ref_to_dims_index(cube, ref_to_slice):
    # Required to handle a mix between types
    if _is_single_item(ref_to_slice):
        ref_to_slice = [ref_to_slice]
    dim_to_slice = []
    for ref in ref_to_slice:
        try:
            # attempt to handle as coordinate
            coord = cube._as_list_of_coords(ref)[0]
            dims = cube.coord_dims(coord)
            if not dims:
                msg = ('Requested an iterator over a coordinate ({}) '
                       'which does not describe a dimension.')
                msg = msg.format(coord.name())
                raise ValueError(msg)
            dim_to_slice.extend(dims)
        except TypeError:
            try:
                # attempt to handle as dimension index
                dim = int(ref)
            except ValueError:
                raise ValueError('{} Incompatible type {} for '
                                 'slicing'.format(ref, type(ref)))
            if dim < 0 or dim > cube.ndim:
                msg = ('Requested an iterator over a dimension ({}) '
                       'which does not exist.'.format(dim))
                raise ValueError(msg)
            dim_to_slice.append(dim)
    dim_to_slice = np.unique(dim_to_slice)
    if len(set(dim_to_slice)) != len(dim_to_slice):
        msg = 'The requested coordinates are not orthogonal.'
        raise ValueError(msg)
    return dim_to_slice


def get_associated_coords(cube, dimensions):
    dims = []
    dim_set = set()
    for dim in dimensions:
        if dim not in dim_set:
            dims.append(dim)
            dim_set.add(dim)
    dim_coords = set(itertools.chain.from_iterable(
        [cube.coords(contains_dimension=i, dim_coords=True)
         for i in dims]
    ))
    aux_coords = set(itertools.chain.from_iterable(
        [cube.coords(contains_dimension=i, dim_coords=False)
         for i in dims]
    ))
    return list(dim_coords), list(aux_coords)


def get_empty_data(shape):
    data = np.empty(shape)
    mask = np.empty(shape, dtype=bool)
    return np.ma.masked_array(data, mask)


def get_slice_spec(cube, ref_to_slice):
    slice_dims = ref_to_dims_index(cube, ref_to_slice)
    slice_shape = tuple(cube.shape[d] for d in slice_dims)
    dim_coords, aux_coords = get_associated_coords(cube, slice_dims)
    return slice_shape, dim_coords, aux_coords


def check_slice_spec(shape, dim_coords):
    if shape is None:
        shape = tuple(c.shape[0] for c in dim_coords)
    if dim_coords is not None:
        for length, coord in zip(shape, dim_coords):
            assert length == coord.shape[0]
    return shape


def index_iterator(dims_to_slice, shape):
    dst_slices = (slice(None, None),) * len(dims_to_slice)
    dims = [1 if n in dims_to_slice else i for n, i in enumerate(shape)]
    for index_tuple in np.ndindex(*dims):
        src_ind = tuple(
            slice(None, None) if n in dims_to_slice else i
            for n, i in enumerate(index_tuple)
        )
        dst_ind = tuple(i for n, i in enumerate(index_tuple)
                        if n not in dims_to_slice) + dst_slices
        yield src_ind, dst_ind


def map_slices(src, func, src_rep, dst_rep):
    ref_to_slice = src_rep.coords(dim_coords=True)
    src_slice_dims = ref_to_dims_index(src, ref_to_slice)
    src_keep_dims = list(set(range(src.ndim)) - set(src_slice_dims))
    src_keep_spec = get_slice_spec(src, src_keep_dims)
    res_shape = src_keep_spec[0] + dst_rep.shape
    dim_coords = src_keep_spec[1] + dst_rep.coords(dim_coords=True)
    dim_coords_and_dims = [(c, i) for i, c in enumerate(dim_coords)]
    dst = iris.cube.Cube(
        data=get_empty_data(res_shape),
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
