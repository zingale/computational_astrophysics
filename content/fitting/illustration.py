# simple example of linear regression.  We make up some "experimental"
# data by calling a random number generator to perturb a line and then
# fit to that line
#
# M. Zingale

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

def y_experiment(a1, a2, sigma, x):
    """ return the experimental data in a linear + random fashion a1
        is the intercept, a2 is the slope, and sigma is the error """

    N = len(x)

    # randn gives samples from the "standard normal" distribution
    np.random.seed(1)
    r = np.random.randn(N)

    y = a1 + a2*x + sigma*r

    return y


def linear_regression(x, y, sigma, sigmaa=0):

    N = len(x)

    S = np.sum(1.0/sigma**2)

    xi_1 = np.sum(x/sigma**2)
    xi_2 = np.sum(x*x/sigma**2)

    eta = np.sum(y/sigma**2)
    mu = np.sum(x*y/sigma**2)

    a2 = (S*mu - xi_1*eta)/(xi_2*S - xi_1**2)
    a1 = (eta*xi_2 - mu*xi_1)/(xi_2*S - xi_1**2)

    chisq = np.sum( (a1 + a2*x - y)**2/sigma**2)
    chisq /= N-2

    if not sigmaa:
        return a1, a2, chisq

    else:
        sigmaFit = np.array([xi_2/(S*xi_2 - xi_1**2),
                                S/(S*xi_2 - xi_1**2)])
        sigmaFit = np.sqrt(sigmaFit)

        return a1, a2, chisq, sigmaFit


N = 10
x = np.linspace(0.0, 100.0, N)


# make up the experimental data with errors
sigma = 25.0*np.ones(N)

y = y_experiment(10.0, 3.0, sigma, x)

#plt.scatter(x, y, color="C0")
plt.errorbar(x, y, yerr=sigma, fmt="o", color="C0", capsize=5)


# label a few points
eps = 0.025*np.max(x)
plt.text(x[0]+eps, y[0], "$(x_1, y_1)$", horizontalalignment="left", verticalalignment="bottom")
plt.text(x[N//2]+eps, y[N//2], "$(x_i, y_i)$", horizontalalignment="left", verticalalignment="bottom")
plt.text(x[N-1]+eps, y[N-1], "$(x_N, y_N)$", horizontalalignment="left", verticalalignment="bottom")


# do the linear regression
a1, a2, chisq, sigmaFit = linear_regression(x, y, sigma, sigmaa=1)

plt.plot(x, a1 + a2*x, color="C1")

plt.text(x[N//2+1], a1 + a2*x[N//2+1], "$Y(x;\{a_j\})$", color="C1", fontsize="medium", horizontalalignment="right", verticalalignment="bottom")


# plot a distance
plt.plot([x[N//2], x[N//2]], [y[N//2], a1 + a2*x[N//2]], color="C6", lw=2, zorder=100)
plt.text(x[N//2]-eps, 0.5*(y[N//2] + a1 + a2*x[N//2]), "$\Delta_i$", color="C6", horizontalalignment="right", verticalalignment="center")
dx = x[1] - x[0]
plt.xlim(x[0]-0.5*dx, x[-1]+1.5*dx)

plt.tight_layout()
plt.savefig("fitting-illustration.png", dpi=150, bbox_inches="tight")
