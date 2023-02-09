import numpy as np
import matplotlib.pyplot as plt
import grid_plot as gp

# plot two stacked fv grids of different (2x) resolution to show prolongation

#-----------------------------------------------------------------------------

nf = 2
nc = nf//2

grf = gp.FVGrid(nf)
grc = gp.FVGrid(nc, voff=2.0)


plt.clf()

grf.draw_grid(center_ticks=False, grid_top=0.5)
grc.draw_grid(center_ticks=False, grid_top=0.5)

grf.label_center(nf//2-1, r"$\tau/2$")
grf.label_center(nf//2,   r"$\tau/2$")

grc.label_center(nc//2-1, r"$\tau$")

grf.label_edge(nf//2-1, r"$y(t)$", offset=0.2, fontsize="medium")
grf.label_edge(nf//2, r"$y(t+\tau)$", right_edge=True, offset=0.2, fontsize="medium")

# connect the dots...

plt.plot([grf.xl[nf//2-1], grf.xl[nf//2-1]], [-0.25, 3.25], ":", color="0.5")
plt.plot([grf.xr[nf//2], grf.xr[nf//2]], [-0.25, 3.25], ":", color="0.5")


plt.axis([grf.xmin-0.5*grf.dx,grf.xmax+0.5*grf.dx, -0.5, 3.5])
plt.axis("off")

plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

f = plt.gcf()
f.set_size_inches(6.0, 2.5)

plt.savefig("adaptive_timestep.png", dpi=200)

