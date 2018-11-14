import numpy as np

def grad_psl(msl):
    gradlatpsl = np.gradient(msl, axis=1)
    gradlonpsl = np.gradient(msl, axis=2)

    return gradlatpsl, gradlonpsl
