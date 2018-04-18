
import numpy as np


def diverging_one_sided(n, style='blue'):

    # Define the colours to be used
    if style == 'blue':
        c_0 = np.array([1.0, 1.0, 1.0])
        c_1 = np.array([1.0, 1.0, 0.0])
        c_2 = np.array([0.0, 1.0, 0.0])
        c_3 = np.array([0.0, 0.7, 0.7])
        c_4 = np.array([0.0, 0.0, 1.0])
        c_5 = np.array([0.0, 0.0, 0.5])
        c_offscale = np.array([1.0, 0.0, 1.0])
        skew_1_2 = -0.4
        skew_2_3 = 0.15
        skew_3_4 = -0.3
        skew_4_5 = 0.0
    elif style == 'red':
        c_0 = np.array([1.0, 1.0, 1.0])
        c_1 = np.array([0.0, 1.0, 1.0])
        c_2 = np.array([0.0, 1.0, 0.0])
        c_3 = np.array([0.8, 0.8, 0.0])
        c_4 = np.array([1.0, 0.0, 0.0])
        c_5 = np.array([0.5, 0.0, 0.0])
        c_offscale = np.array([1.0, 0.0, 1.0])
        skew_1_2 = -0.4
        skew_2_3 = 0.0
        skew_3_4 = -0.4
        skew_4_5 = 0.0
    elif style == 'radar':
        c_0 = np.array([0.0, 0.0, 0.0])
        c_1 = np.array([0.0, 0.0, 1.0])
        c_2 = np.array([0.0, 0.9, 0.0])
        c_3 = np.array([1.0, 1.0, 0.0])
        c_4 = np.array([1.0, 0.0, 0.0])
        c_5 = np.array([1.0, 0.5, 1.0])
        c_offscale = np.array([1.0, 1.0, 1.0])
        skew_1_2 = 0.3
        skew_2_3 = 0.4
        skew_3_4 = -0.4
        skew_4_5 = 0.3

    # Define the positions of the colours in the colour array
    i_1 = 1
    i_2 = int(np.rint(np.nextafter(0.25, 1.0)*(n-1))) + 1
    i_3 = int(np.rint(np.nextafter(0.5, 1.0)*(n-1))) + 1
    i_4 = int(np.rint(np.nextafter(0.75, 1.0)*(n-1))) + 1
    i_5 = n

    # Intialise the colour array
    colours = np.zeros((n+2, 3))

    colours[0, :] = c_0
    colours[n+1, :] = c_offscale

    colours[i_4, :] = c_4

    if i_4 > i_1:

        colours[i_1, :] = c_1

        if i_2 > i_1 and i_2 < i_4:

            colours[i_2, :] = c_2
            if i_2 > i_1+1:
                for i in range(i_1+1, i_2):
                    interp = skew_interp(float(i-i_1)/float(i_2-i_1), skew_1_2)
                    colours[i, :] = (1.0-interp)*c_1 + interp*c_2

            if i_3 > i_2 and i_3 < i_4:

                colours[i_3, :] = c_3
                if i_3 > i_2+1:
                    for i in range(i_2+1, i_3):
                        interp = skew_interp(float(i-i_2)/float(i_3-i_2), skew_2_3)
                        colours[i, :] = (1.0-interp)*c_2 + interp*c_3

                if i_4 > i_3+1:
                    for i in range(i_3+1, i_4):
                        interp = skew_interp(float(i-i_3)/float(i_4-i_3), skew_3_4)
                        colours[i, :] = (1.0-interp)*c_3 + interp*c_4

    if i_5 > i_4:
        colours[i_5, :] = c_5
        if i_5 > i_4+1:
            for i in range(i_4+1, i_5):
                interp = skew_interp(float(i-i_4)/float(i_5-i_4), skew_4_5)
                colours[i, :] = (1.0-interp)*c_4 + interp*c_5

    return colours


def diverging_two_sided(n):

    # Define the colours to be used
    c_0 = np.array([1.0, 1.0, 1.0])
    c_1_pos = np.array([1.0, 1.0, 0.0])
    c_2_pos = np.array([1.0, 0.0, 0.0])
    c_3_pos = np.array([0.5, 0.0, 0.0])
    c_offscale_pos = np.array([1.0, 0.0, 1.0])
    c_1_neg = np.array([0.0, 1.0, 1.0])
    c_2_neg = np.array([0.0, 0.0, 1.0])
    c_3_neg = np.array([0.0, 0.0, 0.5])
    c_offscale_neg = np.array([0.0, 1.0, 0.0])
    skew_1_2_pos = -0.3
    skew_2_3_pos = 0.0
    skew_1_2_neg = -0.3
    skew_2_3_neg = 0.0

    # Define the position of the intermediate colours c_2 in the colour array
    i_2 = int(np.rint((2.0/3.0)*n)) - 1

    # Intialise separate arrays for the positive and negative parts
    colours_pos = np.zeros((n+1, 3))
    colours_neg = np.zeros((n+1, 3))

    # Put the c_2 colours at their assigned locations
    colours_pos[i_2, :] = c_2_pos
    colours_neg[i_2, :] = c_2_neg

    # If the colourbar length is great enough for their to be space...
    if i_2 > 0:
        # Insert the c_1 colours at the lower ends
        colours_pos[0, :] = c_1_pos
        colours_neg[0, :] = c_1_neg
        # If there's space between the c_1 and c_2 colours, interpolate between
        if i_2 > 1:
            for i in range(1, i_2):
                interp0 = float(i)/float(i_2)
                interp = skew_interp(interp0, skew_1_2_pos)
                colours_pos[i, :] = (1.0-interp)*c_1_pos + interp*c_2_pos
                interp = skew_interp(interp0, skew_1_2_neg)
                colours_neg[i, :] = (1.0-interp)*c_1_neg + interp*c_2_neg

    # If the colourbar length is great enough for their to be space...
    if n-1 > i_2:
        # Put the c_3 colours at their assigned locations
        colours_pos[n-1, :] = c_3_pos
        colours_neg[n-1, :] = c_3_neg
        # If there's space between the c_2 and c_3 colours, interpolate between
        if n-1 > i_2+1:
            for i in range(i_2+1, n-1):
                interp0 = float(i-i_2)/float(n-1-i_2)
                interp = skew_interp(interp0, skew_2_3_pos)
                colours_pos[i, :] = (1.0-interp)*c_2_pos + interp*c_3_pos
                interp = skew_interp(interp0, skew_2_3_neg)
                colours_neg[i, :] = (1.0-interp)*c_2_neg + interp*c_3_neg

    # Add in the offscale colours at the far extremes
    colours_pos[n, :] = c_offscale_pos
    colours_neg[n, :] = c_offscale_neg

    colours = np.concatenate((colours_neg[::-1], c_0[None, :], colours_pos))

    # Return outputs
    return colours


def diverging_two_sided_cyclic(n):

    # Define the colours to be used
    c_0 = np.array([1.0, 1.0, 1.0])
    c_1_pos = np.array([1.0, 1.0, 0.0])
    c_2_pos = np.array([1.0, 0.0, 0.0])
    c_1_neg = np.array([0.0, 1.0, 1.0])
    c_2_neg = np.array([0.0, 0.0, 1.0])
    c_3 = np.array([0.5, 0.0, 0.5])
    skew_1_2_pos = -0.3
    skew_2_3_pos = 0.1
    skew_1_2_neg = -0.3
    skew_2_3_neg = 0.1

    # Define the position of the intermediate colours c_2 in the colour array
    i_2 = int(np.rint((2.0/3.0)*n)) - 1

    colours_pos = np.zeros((n+1, 3))
    colours_neg = np.zeros((n+1, 3))

    colours_pos[i_2, :] = c_2_pos
    colours_neg[i_2, :] = c_2_neg

    if i_2 > 0:
        colours_pos[0, :] = c_1_pos
        colours_neg[0, :] = c_1_neg
        if i_2 > 1:
            for i in range(1, i_2):
                interp0 = float(i)/float(i_2)
                interp = skew_interp(interp0, skew_1_2_pos)
                colours_pos[i, :] = (1.0-interp)*c_1_pos + interp*c_2_pos
                interp = skew_interp(interp0, skew_1_2_neg)
                colours_neg[i, :] = (1.0-interp)*c_1_neg + interp*c_2_neg

    if i_2 < n:
        for i in range(i_2+1, n):
            interp0 = float(i-i_2)/(float(n-i_2)-0.5)
            interp = skew_interp(interp0, skew_2_3_pos)
            colours_pos[i, :] = (1.0-interp)*c_2_pos + interp*c_3
            interp = skew_interp(interp0, skew_2_3_neg)
            colours_neg[i, :] = (1.0-interp)*c_2_neg + interp*c_3

    colours_pos[n] = colours_neg[n-1]
    colours_neg[n] = colours_pos[n-1]

    colours = np.concatenate((colours_neg[::-1], c_0[None, :], colours_pos))

    # Return outputs
    return colours


def cyclic(n):

    # Setup the 6 basic colours on which this is based
    n_col = 6
    cols = np.zeros([n_col+1, 3])
#    # Colours chosen as a compromise between constant intensity and ease
#    # of seeing gradients around the colourbar
#    cols[0:n_col, :] = np.asarray([[0.0, 0.0, 1.0],  # blue
#                                   [0.0, 0.5, 0.5],
#                                   [0.0, 0.8, 0.0],  # green
#                                   [0.9, 0.9, 0.0],
#                                   [1.0, 0.0, 0.0],  # red
#                                   [0.7 ,0.0, 0.7]])
#    skew = np.array([0.2, 0.2, 0.3, -0.4, 0.0, -0.2])
    # Colours chosen to maximise gradients
    cols[0:n_col, :] = np.asarray([[0.0, 0.0, 1.0],  # blue
                                   [0.0, 1.0, 1.0],
                                   [0.0, 0.5, 0.0],  # green
                                   [1.0, 1.0, 0.0],
                                   [1.0, 0.0, 0.0],  # red
                                   [1.0, 0.7, 1.0]])
    skew = np.array([0.2, 0.0, 0.2, -0.4, 0.4, -0.3])
    cols[-1, :] = cols[0, :]

    r_col = float(n_col) * np.arange(n)/float(n)
    red = np.zeros(n)
    grn = np.zeros(n)
    blu = np.zeros(n)

    for i_col in range(n_col):
        i = np.logical_and(r_col >= i_col, r_col < i_col+1)
        interp = skew_interp(r_col[i]-float(i_col), skew[i_col])
        red[i] = (1.0-interp)*cols[i_col, 0] + interp*cols[i_col+1, 0]
        grn[i] = (1.0-interp)*cols[i_col, 1] + interp*cols[i_col+1, 1]
        blu[i] = (1.0-interp)*cols[i_col, 2] + interp*cols[i_col+1, 2]

    colours = np.asarray([red, grn, blu]).T
    colours = np.concatenate((colours, colours[0, :][None, :]))

    return colours


def add_intensity_dimension(colours1d, I_vals):

    import numpy as np

    # I = 0 : black
    # I = 1 : colour unchanged
    # I = 2 : white

    # Convert I values to array with atleast 1 dimension
    I_vals = np.atleast_1d(I_vals)

    # Get array sizes
    n_c = np.shape(colours1d)[0]
    n_I, = np.shape(I_vals)

    # Initialise the output 2d colour array as a tiled mesh of colours1d
    colours2d = np.tile(colours1d[None, :, :], (n_I, 1, 1))
    # Also create a tiled mesh of the input intensity values
    I = np.tile(I_vals[:, None], (1, n_c))

    # Calculate array of intensity-scaled colours using scale_intensity
    colours2d = scale_intensity(colours2d, I)

    return colours2d


def scale_intensity(colours, I):

    # I = 0 : black
    # I = 1 : colour unchanged
    # I = 2 : white

    if not np.shape(colours) == np.shape(I)+(3, ):
        print 'Error: Input array shape mismatch.'

    # Roll the RGB axis of the input colour array to the leading subscript
    rnk = np.rank(colours)
    colours = np.rollaxis(colours, rnk-1)

    # Calculate the sum intensity of each input colour
    I_orig = np.sum(colours, axis=0)

    # Calculate a multiplier, using skew_interp, such that the sum intensity
    # of the new colours will be 0 if I=0, unchanged if I=1, and 3 if I=2
    mult = skew_interp(I/2.0, (I_orig-1.5)/1.5) * 3.0 / I_orig

    # Apply the multiplier to each RGB component of the input colours
    for i in range(3):
        colours[i] = colours[i] * mult

    # Where any colour has over-saturated, this code shifts it towards
    # the 3D diagonal by the amount needed to restore it to one, whilst
    # conserving the overall intensity R+G+B
    # For each colour...
    for i in range(3):

        # Find any points where this colour has saturated (>1.0)
        i_sat = np.nonzero(colours[i] > 1.0)
        # For those points, calculate nearest colour on the diagonal
        meancol = (colours[0][i_sat] + colours[1][i_sat] + colours[2][i_sat]) / 3.0

        # Rounding errors can cause meancol>1, and occasionally div-by-zero.
        # Prevent this by enforcing meancol<=1
        i_sat2 = np.nonzero(meancol > 1.0)
        if len(i_sat2) > 0:
            for i2 in range(3):
                colours[i2][i_sat][i_sat2] = colours[i2][i_sat][i_sat2] / meancol[i_sat2]
            meancol[i_sat2] = 1.0

        # Calculate how far along the line towards meancol we should move
        dcol = (colours[i][i_sat] - 1.0) / (colours[i][i_sat] - meancol)
        # Apply the vector increment
        for i2 in range(3):
            colours[i2][i_sat] = colours[i2][i_sat] + dcol * (meancol - colours[i2][i_sat])

    # Clip output to between 0 and 1 to remove miniscule oversaturations
    # resulting from rounding errors
    colours = np.clip(colours, 0.0, 1.0)

    # Roll the RGB axis back to the end
    colours = np.rollaxis(colours, 0, start=rnk)

    return colours


# Take a value x between 0 and 1, and skew it towards 0 or 1 by an amount skew.
# y=x for skew=0,  y=0 for skew=-1,  y=1 for skew=1,  varies smoothly between.
def skew_interp(x, skew):

    # Convert skew (between -1 and +1) to half-way value (between 0 and 1)
    y_half = 0.5*(1.0+skew)

    # Formula to apply smoothly varying skew (with clipping of the denominator
    # to avoid div-by-zero error when x=0 and y_half=1)
    y = y_half * x / np.clip(1.0-y_half+(2.0*y_half-1.0)*x, np.nextafter(0.0, 1.0), 1.0)

    # Clip the output to within the allowed range 0:1 (rounding errors in the
    # above calculation occasionally yield values slightly >1)
    y = np.clip(y, 0.0, 1.0)

    return y


def colourbar_test_plot(colours):

    import matplotlib.pyplot as plt
    import matplotlib.collections as mcoll

    if np.rank(colours) == 1:
        colours = colours[None, :]
    if np.rank(colours) == 2:
        colours = colours[None, :, :]

    ny, nx, n = np.shape(colours)

    (xc, yc) = np.meshgrid(np.asarray(np.arange(nx+1), dtype=float)/nx,
                           np.asarray(np.arange(ny+1), dtype=float)/ny)
    xcyc = np.array([np.ravel(xc), np.ravel(yc)]).T

    plt.figure()
    axeid = plt.axes()
    qm = mcoll.QuadMesh(nx, ny, xcyc, color=np.reshape(colours, (ny*nx, 3)))
    axeid.add_collection(qm)
