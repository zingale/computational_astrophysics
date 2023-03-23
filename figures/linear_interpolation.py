import matplotlib.pyplot as plt
import numpy as np


def f(x):
    """the actual function we will sample"""
    return x * (0.5 * np.exp(-(x - 1.0)**2))


x1 = 0.75
x2 = 1.5

xfine = np.linspace(0.0, 2.5, 100)

fig, ax = plt.subplots()
ax.plot(xfine, f(xfine))


# interpolant

m = (f(x2) - f(x1)) / (x2 - x1)
x_interp = np.linspace(x1, x2, 10)
ax.plot(x_interp, m * (x_interp - x1) + f(x1), color="C1")
ax.scatter([x1, x2], [f(x1), f(x2)], marker="o", color="C1")

eps = 0.02
ax.text(x1, f(x1)-eps, "$f_1$", color="C1", fontsize="large")
ax.text(x2+eps, f(x2), "$f_2$", color="C1", fontsize="large")

ax.plot([x1, x1], [0, f(x1)], color="0.5", ls=":")
ax.plot([x2, x2], [0, f(x2)], color="0.5", ls=":")

ax.text(x1, -eps, "$x_1$", color="0.5", fontsize="large")
ax.text(x2, -eps, "$x_2$", color="0.5", fontsize="large")

ax.spines[['left', 'bottom']].set_position('zero')
ax.spines[['top', 'right']].set_visible(False)

fig.savefig("linear-interpolation.png")
