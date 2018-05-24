'''
Module docstring
'''

import numpy as np

def amplitude_and_phase(tseries):
    '''
    Routine docstring
    '''

    # Get time dimension (assumed to be the leading dimension)
    ntimes = np.shape(tseries)[0]

    # Calculate a typical amplitude; we will detect "in the noise"
    #  frequencies relative to this.
    #rms_value = np.std(tseries, axis=0)

    # Use FFT to do the Fourier transform
    ftrans = np.fft.rfft(tseries, axis=0)

    # Calculate the amplitude for each frequency
    amplitude = np.abs(ftrans)
    # The normalisation is by n for the mean, and n/2 for non-zero frequencies.
    amplitude[0] = amplitude[0]/ntimes
    amplitude[1:] = amplitude[1:]/(ntimes/2)

    # Calculate the phase for each frequency (change sign to get phase of the
    #  maximum from zero in the positive direction).
    phase = -np.angle(ftrans)

    # Shift negative phases to positive ones greater than pi.
    l_negative = phase < 0.0
    phase[l_negative] = phase[l_negative] + 2.0*np.pi

    return (amplitude, phase)
