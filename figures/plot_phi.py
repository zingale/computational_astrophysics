# plot the analytic solution at various times and t0

import matplotlib.pyplot as plt
import numpy as np

def phi_a(t, x, xc, k, t0, phi1, phi2):
    """ analytic solution """

    return (phi2 - phi1)*np.sqrt(t0/(t + t0)) * \
        np.exp(-0.25*(x-xc)**2/(k*(t + t0))) + phi1


phi1 = 1
phi2 = 2

k = 1

npts = 16384
xmin = 0.0
xmax = 1.0

xc = 0.5*(xmin + xmax)
dx = (xmax - xmin)/npts

x = (np.arange(npts) + 0.5)*dx + xmin


#-----------------

t0 = 1.e-8
p8 = phi_a(0.0, x, xc, k, t0, phi1, phi2)

plt.plot(x, p8, label="$t_0 = 10^{-8}$")


t0 = 1.e-6
p6 = phi_a(0.0, x, xc, k, t0, phi1, phi2)

plt.plot(x, p6, label="$t_0 = 10^{-6}$")


t0 = 1.e-4
p4 = phi_a(0.0, x, xc, k, t0, phi1, phi2)

plt.plot(x, p4, label="$t_0 = 10^{-4}$")

plt.xlim(0.4,0.6)

plt.xlabel("$x$")
plt.ylabel("$\phi_0$")

plt.legend(frameon=False)

plt.savefig("phi_analytic.png", dpi=150)


#-----------------

plt.clf()

t0 = 1.e-4

pa = phi_a(0.0, x, xc, k, t0, phi1, phi2)
plt.plot(x, pa, label="$t = 0$")

pa = phi_a(0.0005, x, xc, k, t0, phi1, phi2)
plt.plot(x, pa, label="$t = 0.0005$")

pa = phi_a(0.005, x, xc, k, t0, phi1, phi2)
plt.plot(x, pa, label="$t = 0.005$")



plt.xlim(0.3,0.7)

plt.xlabel("$x$")
plt.ylabel("$\phi_0$")

plt.legend(frameon=False)

plt.savefig("phi_analytic_t.png", dpi=150)
