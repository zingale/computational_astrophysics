# plot the frequencies that are represented by the FFT for N points

import numpy as np
import matplotlib.pyplot as plt

xmin = 0.0
xmax = 1.0
L = xmax - xmin

def f(k, x):
    return np.sin(2.0 * np.pi * k * x /L)


N = 128

x = np.linspace(0.0, 1.0, N, endpoint=False)  # N points, not including end point

fig, ax = plt.subplots(nrows=2, ncols=2)

ax = ax.flatten()

# low k's
ks = [1, 2, 4, 10]

for n in range(len(ks)):
    ax[n].plot(x, f(ks[n], x))
    ax[n].set_title(rf"$k = {ks[n]},\, \nu = {ks[n]}/L$")
    ax[n].set_xlabel(r"$x$")
    ax[n].set_xmargin(0.0)

fig.tight_layout()
fig.savefig("sine-k.png")

