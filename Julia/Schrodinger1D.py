from scipy import integrate
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML
plt.rc('savefig', dpi=300)

import numpy as np

from numba import jit, njit

dx = 0.02
x = np.arange(0, 10, dx)


kx = 0.1
m = 1
sigma = 0.1
