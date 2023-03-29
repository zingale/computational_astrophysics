import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots()


# draw a grid

ax.scatter([1, 1, 2, 2], [1, 2, 1, 2], marker="o")
ax.scatter([1.25], [1.25], marker="x")

eps = 0.075

ax.text(1+eps, 1-eps, r"$P(\rho_i, T_j)$")
ax.text(2+eps, 1-eps, r"$P(\rho_{i+1}, T_j)$")
ax.text(1+eps, 2+eps, r"$P(\rho_i, T_{j+1})$")
ax.text(2+eps, 2+eps, r"$P(\rho_{i+1}, T_{j+1})$")

ax.plot([0.5, 2.5], [1, 1], color="0.5", ls=":")
ax.plot([0.5, 2.5], [2, 2], color="0.5", ls=":")

ax.plot([1, 1], [0.5, 2.5], color="0.5", ls=":")
ax.plot([2, 2], [0.5, 2.5], color="0.5", ls=":")

ax.margins(0.0)

ax.axis("off")
ax.set_aspect("equal", "datalim")

fig.set_size_inches(4, 4)

fig.tight_layout()

fig.savefig("bilinear_eos.png", dpi=150, bbox_inches="tight")

