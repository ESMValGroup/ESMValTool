
import iris
import matplotlib.collections as mcoll
import matplotlib.pyplot as plt
import numpy as np


def pcolor(c_data, y, x, c_intervals, colourbar, ColourbarPos=None):

    #---------------------------------------
    # Check the input arrays and coordinates
    #---------------------------------------

    # Check that the input data array is 2-dimensional
    if not (np.rank(c_data) == 2):
        print 'Error: input data array must be 2-dimensional'
    # Find array dimensions
    (ny, nx) = np.shape(c_data)

    # Check that the interval vector is 1-D and find its length
    c_intervals = np.atleast_1d(c_intervals)
    if not (np.rank(c_intervals) == 1):
        print 'Error: input interval vector must be 1-dimensional'
    (nc, ) = np.shape(c_intervals)

    # Check that the dimensions of the colourbar match the interval vector
    if not np.rank(colourbar) == 2:
        print 'Error: colourbar input must be 2-dimensional'
    if not np.shape(colourbar) == (nc+1, 3):
        print "Error: the input colourbar's shape doesn't match the input interval vector"

    #-----------------------------------------------------------
    # Set up the vertex coordinates and colour data for QuadMesh
    #-----------------------------------------------------------

    # Find the bounds of each coordinate
    x_corners = calc_coord_bounds(c_data, x, 1)
    y_corners = calc_coord_bounds(c_data, y, 0)

    # If the input is a cube, extract the data from it
    if isinstance(c_data, iris.cube.Cube):
        c_data = c_data.data.copy()

    # Create an array containing a list of the coordinates of the vertices
    # of every square to be rendered
    (x_mesh, y_mesh) = np.meshgrid(x_corners, y_corners)
    xy_list = np.array([np.ravel(x_mesh), np.ravel(y_mesh)]).T

    # Create a similarly unraveled version of the data array
    c_1d = np.ravel(c_data)

    # For each square in the list of points, use its value of c_data and
    # (relative to the input intervals in c) to find
    # subscripts to get corresponding colours from the colourbar...
    i_cb = np.zeros(ny*nx, dtype=int)
    i = c_1d <= c_intervals[0]
    i_cb[i] = 0
    for ic in range(1, nc):
        i = np.logical_and(c_1d > c_intervals[ic-1], c_1d <= c_intervals[ic])
        i_cb[i] = ic
    i = c_1d > c_intervals[nc-1]
    i_cb[i] = nc

    #--------------
    # Make the plot
    #--------------

    # Set the axes limits
    plt.xlim([x_corners[0], x_corners[-1]])
    plt.ylim([y_corners[0], y_corners[-1]])
    # Create a pcolor-style plot using QuadMesh, and add it to the current axes
    qm = mcoll.QuadMesh(nx, ny, xy_list, color=colourbar[i_cb, :])
    plt.gca().add_collection(qm)
    plt.draw()

    #-------------------
    # Make the colourbar
    #-------------------

    if not ColourbarPos is None:

        if not np.rank(ColourbarPos) == 1:
            print 'Error: input colourbar postion data not a vector'
        if not np.shape(ColourbarPos) == (4, ):
            print 'Error: input colourbar position data must contain 4 numbers'
        try:
            ColourbarPos = np.asarray(ColourbarPos, dtype=float)
        except ValueError:
            print 'Error: input colourbar position data must contain numeric values'

        # Create a new set of axes for the colourbar
        cbar_axes = plt.axes(ColourbarPos)
        # Define psuedo-coordinate data for a QuadMesh plot of the colourbar
        c_cbar = np.arange(nc+2)
        (c_mesh, y_mesh) = np.meshgrid(c_cbar, np.arange(2))
        cy_list = np.array([np.ravel(c_mesh), np.ravel(y_mesh)]).T
        # Create the QuadMesh of the colourbar, and add it to the axes
        qm_cbar = mcoll.QuadMesh(nc+1, 1, cy_list, color=colourbar)
        cbar_axes.add_collection(qm_cbar)
        # Add ticks and labels on the axes, marking the input intervals
        plt.setp(cbar_axes, xticks=c_cbar[1:nc+1])
        plt.setp(cbar_axes, xticklabels=np.asarray(c_intervals, dtype=str))
        plt.setp(cbar_axes, yticks=[])
        # Set colourbar axis limits
        plt.xlim([c_cbar[0], c_cbar[-1]])
        plt.ylim([0, 1])

    return qm


def pcolor_2d_colourbar(c1_data, c2_data, y, x, c1_intervals, c2_intervals,
                        colourbar2d, ColourbarPos=None):

    #---------------------------------------
    # Check the input arrays and coordinates
    #---------------------------------------

    # Check that input arrays are 2-dimensional
    if not (np.rank(c1_data) == 2 and np.rank(c2_data) == 2):
        print 'Error: input arrays must be 2-dimensional'
    # Find array dimensions and check that both inputs are the same shape
    (ny, nx) = np.shape(c1_data)
    if not np.shape(c2_data) == (ny, nx):
        print 'Error: the 2 input data arrays have different shapes'

    # Check that the interval vectors are 1-D and find their lengths
    c1_intervals = np.atleast_1d(c1_intervals)
    c2_intervals = np.atleast_1d(c2_intervals)
    if not (np.rank(c1_intervals) == 1 and np.rank(c2_intervals) == 1):
        print 'Error: input interval vectors must be 1-dimensional'
    (n_c1, ) = np.shape(c1_intervals)
    (n_c2, ) = np.shape(c2_intervals)

    # Check that the dimensions of the colourbar match the interval vectors
    if not np.rank(colourbar2d) == 3:
        print 'Error: colourbar input must be 3-dimensional'
    if not np.shape(colourbar2d) == (n_c1+1, n_c2+1, 3):
        print "Error: the input colourbar's shape doesn't match the input" +\
            " interval vectors"

    #-----------------------------------------------------------
    # Set up the vertex coordinates and colour data for QuadMesh
    #-----------------------------------------------------------

    # Find the bounds of each coordinate
    x_corners = calc_coord_bounds(c1_data, x, 1)
    y_corners = calc_coord_bounds(c1_data, y, 0)

    # If the inputs are cubes, extract the data from them
    if isinstance(c1_data, iris.cube.Cube):
        c1_data = c1_data.data.copy()
    if isinstance(c2_data, iris.cube.Cube):
        c2_data = c2_data.data.copy()

    # Create an array containing a list of the coordinates of the vertices
    # of every square to be rendered
    (x_mesh, y_mesh) = np.meshgrid(x_corners, y_corners)
    xy_list = np.array([np.ravel(x_mesh), np.ravel(y_mesh)]).T

    # Create similarly unraveled versions of the data arrays
    c1_1d = np.ravel(c1_data)
    c2_1d = np.ravel(c2_data)

    # For each square in the list of points, use its values of c1_data and
    # c2_data (relative to the input intervals in c1 and c2) to find
    # subscripts to get corresponding colours from the colourbar...
    # For c1_data:
    i_c1 = np.zeros(ny*nx, dtype=int)
    i = c1_1d < c1_intervals[0]
    i_c1[i] = 0
    for ic in range(1, n_c1):
        i = np.logical_and(c1_1d >= c1_intervals[ic-1], c1_1d < c1_intervals[ic])
        i_c1[i] = ic
    i = c1_1d >= c1_intervals[n_c1-1]
    i_c1[i] = n_c1
    # For c2_data:
    i_c2 = np.zeros(ny*nx, dtype=int)
    i = c2_1d < c2_intervals[0]
    i_c2[i] = 0
    for ic in range(1, n_c2):
        i = np.logical_and(c2_1d >= c2_intervals[ic-1], c2_1d < c2_intervals[ic])
        i_c2[i] = ic
    i = c2_1d >= c2_intervals[n_c2-1]
    i_c2[i] = n_c2

    #--------------
    # Make the plot
    #--------------

    # Set the axes limits
    plt.xlim([x_corners[0], x_corners[-1]])
    plt.ylim([y_corners[0], y_corners[-1]])
    # Create a pcolor-style plot using QuadMesh, and add it to the current axes
    qm = mcoll.QuadMesh(nx, ny, xy_list, color=colourbar2d[i_c1, i_c2, :])
    plt.gca().add_collection(qm)
    plt.draw()

    #-------------------
    # Make the colourbar
    #-------------------

    if not ColourbarPos is None:

        if not np.rank(ColourbarPos) == 1:
            print 'Error: input colourbar postion data not a vector'
        if not np.shape(ColourbarPos) == (4, ):
            print 'Error: input colourbar position data must contain 4 numbers'
        try:
            ColourbarPos = np.asarray(ColourbarPos, dtype=float)
        except ValueError:
            print 'Error:  input colourbar ' +\
                'position data must contain numeric values'

        # Create a new set of axes for the colourbar
        cbar_axes = plt.axes(ColourbarPos)
        # Define psuedo-coordinate data for a QuadMesh plot of the colourbar
        c1_cbar = np.arange(n_c1+2)
        c2_cbar = np.arange(n_c2+2)
        (c2_mesh, c1_mesh) = np.meshgrid(c2_cbar, c1_cbar)
        c2c1_list = np.array([np.ravel(c2_mesh), np.ravel(c1_mesh)]).T
        # Create the QuadMesh of the colourbar, and add it to the axes
        qm_cbar = mcoll.QuadMesh(n_c2+1, n_c1+1, c2c1_list,
                                 color=np.reshape(colourbar2d,
                                                  ((n_c2+1)*(n_c1+1), 3)))
        cbar_axes.add_collection(qm_cbar)
        # Add ticks and labels on the axes, marking the input intervals
        plt.setp(cbar_axes, xticks=c2_cbar[1:n_c2+1])
        plt.setp(cbar_axes, xticklabels=np.asarray(c2_intervals, dtype=str))
        plt.setp(cbar_axes, yticks=c1_cbar[1:n_c1+1])
        plt.setp(cbar_axes, yticklabels=np.asarray(c1_intervals, dtype=str))
        # Set colourbar axis limits
        plt.xlim([c2_cbar[0], c2_cbar[-1]])
        plt.ylim([c1_cbar[0], c1_cbar[-1]])

    return qm


def calc_coord_bounds(cdata, v, v_rank):

    nv = np.shape(cdata)[v_rank]
    v_corners = np.zeros(nv+1)

    # If the input is an iris cube, v should be a coord name string
    if isinstance(cdata, iris.cube.Cube):

        if not isinstance(v, str):
            print 'Error: for cube input, coord inputs must be strings ' +\
                  'containing cube coordinate names'

        # Check that the coordinate is a 1-D vector and find its bounds
        if np.rank(cdata.coord(v).points) == 1:
            if cdata.coord(v).bounds is None:
                cdata.coord(v).guess_bounds()
            v_corners[0:nv] = cdata.coord(v).bounds[:, 0]
            v_corners[nv] = cdata.coord(v).bounds[nv-1, 1]
        else:
            print "Error: this routine doesn't currently support " +\
                  "coordinate vectors which aren't 1-dimensional"

    # If the input isn't a cube, interpret v as coordinate vector data
    else:

        try:
            v = np.atleast_1d(v, dtype=float)
        except ValueError:
            print 'Error: for non-cube input, coord inputs must contain ' +\
                  'numeric  coordinate data'

        if np.rank(v) == 1:
            # If length of v equals nv, interpolate to find bounds
            if np.shape(v)[0] == nv:
                v_corners[1:nv] = 0.5*(v[1:nv]+v[0:nv-1])
                v_corners[0] = v_corners[1] - (v_corners[2]-v_corners[1])
                v_corners[nv] = v_corners[nv-1] + (v_corners[nv-1]-v_corners[nv-2])
            # If length of v equals nv+1, assume it is already bounds
            elif np.shape(v)[0] == nv+1:
                v_corners[:] = v[:]
            else:
                print "Error: length of coordinate doesn't match the " +\
                      "corresponding dimension of the data array"
        else:
            print "Error: this routine doesn't currently support " +\
                  "coordinate vectors which aren't 1-dimensional"

    return v_corners
