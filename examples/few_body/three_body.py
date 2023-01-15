# this comes from Newman 8.16

import numpy as np
import matplotlib.pyplot as plt

class Star(object):
    """ a container to hold the star positions """

    def __init__(self):
        self.t = []
        self.x = []
        self.y = []
        self.u = []
        self.v = []


class ThreeBody(object):
    """ model the gravitational interaction of three stars """

    def __init__(self, M1, pos1, M2, pos2, M3, pos3, SMALL=1.e-30):
        """ pos is a tuple of x-y coords of the stars """

        self.G = 1  # dimensionless units

        self.pos0 = []
        self.pos0 += list(pos1)
        self.pos0 += list(pos2)
        self.pos0 += list(pos3)
        self.pos0 = np.asarray(self.pos0, dtype=np.float64)

        self.M = []
        self.M.append(M1)
        self.M.append(M2)
        self.M.append(M3)
        self.M = np.array(self.M, dtype=np.float64)

        self.stars = []

        self.N = len(self.M)

        self.SMALL = SMALL

    def time(self):
        """return the time that the steps are defined at --- it is the same
        for all stars"""
        return self.stars[0].t

    def npts(self):
        """ return the number of integration points """
        return len(self.stars[0].t)

    def center_of_mass(self):
        """return the center of mass of the system at all time points"""
        Mtot = np.sum(self.M)

        x_cm = np.zeros(self.npts(), dtype=np.float64)
        y_cm = np.zeros(self.npts(), dtype=np.float64)

        for n in range(self.N):
            M = self.M[n]
            x_cm[:] = x_cm[:] + M*np.asarray(self.stars[n].x[:])
            y_cm[:] = y_cm[:] + M*np.asarray(self.stars[n].y[:])

        x_cm /= Mtot
        y_cm /= Mtot

        return x_cm, y_cm

    def energy(self):
        """return the energy of the system vs time"""

        E = np.zeros(self.npts(), dtype=np.float64)

        for n in range(self.npts()):

            # kinetic energy

            KE = [0.5 * self.M[i] * (self.stars[i].u[n]**2 +
                                     self.stars[i].v[n]**2)
                  for i in range(self.N)]

            # potential energy -- we need all pairs, but order doesn't matter

            PE = [-self.G * self.M[i] * self.M[j] /
                  np.sqrt((self.stars[i].x[n] - self.stars[j].x[n])**2 +
                          (self.stars[i].y[n] - self.stars[j].y[n])**2)
                  for i in range(self.N) for j in range(self.N) if i > j]

            E[n] = sum(KE) + sum(PE)

        return E

    def rhs(self, xvec, vvec):
        """ x vec is (x1, y1, x2, y2, x3, y3) """

        xdot = []
        vdot = []
        for n in range(len(self.M)):

            xs = xvec[2*n]
            ys = xvec[2*n+1]

            us = vvec[2*n]
            vs = vvec[2*n+1]

            # compute the space derivatives
            xdot.append(us)
            xdot.append(vs)

            # compute the velocity derivatives
            vxdot = 0.0
            vydot = 0.0

            for m in range(len(self.M)):
                if n == m:
                    continue
                else:
                    xo = xvec[2*m]
                    yo = xvec[2*m+1]

                    r = np.sqrt((xo - xs)**2 + (yo - ys)**2) + self.SMALL
                    vxdot += self.G*self.M[m]*(xo - xs)/r**3
                    vydot += self.G*self.M[m]*(yo - ys)/r**3

            vdot.append(vxdot)
            vdot.append(vydot)

        return np.asarray(xdot), np.asarray(vdot)


    def single_step(self, xvec, vvec, t, dt):
        """ take a single RK-4 timestep from t to t+dt for the system
        ydot = rhs """

        # get the RHS at several points
        xdot1, vdot1 = self.rhs(xvec, vvec)
        xdot2, vdot2 = self.rhs(xvec + 0.5*dt*xdot1, vvec + 0.5*dt*vdot1)
        xdot3, vdot3 = self.rhs(xvec + 0.5*dt*xdot2, vvec + 0.5*dt*vdot2)
        xdot4, vdot4 = self.rhs(xvec + dt*xdot3, vvec + dt*vdot3)

        # advance
        xvec_new = xvec + (dt/6.0)*(xdot1 + 2.0*xdot2 + 2.0*xdot3 + xdot4)
        vvec_new = vvec + (dt/6.0)*(vdot1 + 2.0*vdot2 + 2.0*vdot3 + vdot4)

        return xvec_new, vvec_new

    def integrate(self, dt, err, tmax):
        """ integrate the equations of motion using 4th order R-K
            method with an adaptive stepsize, to try to achieve the
            relative error err.  dt here is the initial timestep

            if err < 0, then we don't do adaptive stepping, but rather
            we always walk at the input dt
            """

        # safety parameters
        S1 = 0.9
        S2 = 4.0

        # initial conditions
        t = 0.0
        xvec = np.asarray(self.pos0)
        vvec = np.zeros_like(xvec)

        for n in range(len(self.M)):
            star = Star()
            star.t.append(t)
            star.x.append(xvec[2*n])
            star.y.append(xvec[2*n+1])
            star.u.append(vvec[2*n])
            star.v.append(vvec[2*n+1])
            self.stars.append(star)

        # start with the old timestep
        dt_new = dt

        n_reset = 0

        while t < tmax:

            # adaptive stepping iteration loop -- keep trying to take
            # a step until we achieve our desired error
            rel_error = 1.e10

            n_try = 0
            while rel_error > err:
                dt = dt_new
                if t+dt > tmax:
                    dt = tmax-t

                # take 2 half steps
                xvec_tmp, vvec_tmp = self.single_step(xvec, vvec, t, 0.5*dt)
                xvec_new, vvec_new = self.single_step(xvec_tmp, vvec_tmp, t+0.5*dt, 0.5*dt)

                # now take just a single step to cover dt
                xvec_single, vvec_single = self.single_step(xvec, vvec, t, dt)

                # xvec_new, vvec_new should be more accurate than
                # xvec_single, vvec_single, since it used smaller
                # steps.

                # estimate the relative error now
                ex = np.abs((xvec_new - xvec_single)/(xvec_single + self.SMALL))
                ev = np.abs((vvec_new - vvec_single)/(vvec_single + self.SMALL))
                rel_error = max(ex.max(), ev.max())

                # adaptive timestep algorithm from Garcia (Eqs. 3.30
                # and 3.31)
                dt_est = dt*abs(err/rel_error)**0.2
                dt_new = min(max(S1*dt_est, dt/S2), S2*dt)

                n_try += 1

            if n_try > 1:
                # n_try = 1 if we took only a single try at the step
                n_reset += (n_try-1)

            # successful step
            t += dt

            # set for the next step
            xvec[:] = xvec_new[:]
            vvec[:] = vvec_new[:]

            # store
            for n in range(len(self.M)):
                self.stars[n].t.append(t)
                self.stars[n].x.append(xvec[2*n])
                self.stars[n].y.append(xvec[2*n+1])
                self.stars[n].u.append(vvec[2*n])
                self.stars[n].v.append(vvec[2*n+1])


if __name__ == "__main__":

    tb = ThreeBody(150., (3, 1),
                   200., (-1, -2),
                   250., (-1, 1))

    tb.integrate(0.05, 1.e-8, 2.0)

    print(f"number of steps: {tb.npts()}")

    xmin = min([min(q.x) for q in tb.stars])
    xmax = max([max(q.x) for q in tb.stars])
    ymin = min([min(q.y) for q in tb.stars])
    ymax = max([max(q.y) for q in tb.stars])

    fig = plt.figure()
    ax = fig.add_subplot()

    for i in range(tb.N):
        ax.plot(tb.stars[i].x, tb.stars[i].y)
        ax.scatter([tb.stars[i].x[0]], [tb.stars[i].y[0]], marker="x")

    fig.savefig("three_body.png")

    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot()

    t = tb.time()
    x_cm, y_cm = tb.center_of_mass()

    ax.plot(t, x_cm, label="x center of mass")
    ax.plot(t, y_cm, label="y center of mass")

    fig.savefig("center_of_mass.png")

    plt.close(fig)

    # timestep plot

    t = np.array(tb.stars[0].t)
    dt = t[1:] - t[:-1]
    ts = 0.5 * (t[1:] + t[:-1])

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.plot(ts, dt)
    ax.set_yscale("log")

    ax.set_xlabel("t")
    ax.set_ylabel("dt")

    fig.savefig("three_body_dts.png")

    plt.close(fig)

    # energy

    fig = plt.figure()
    ax = fig.add_subplot()

    E = tb.energy()

    ax.plot(t, np.abs(E))

    ax.set_yscale("log")

    ax.set_xlabel("t")
    ax.set_ylabel("|E|")

    fig.tight_layout()
    fig.savefig("three_body_energy.png")
