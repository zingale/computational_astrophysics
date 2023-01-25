import numpy as np
import matplotlib.pyplot as plt

class State:
    """ a container to hold a single star's position """

    def __init__(self, x, y, u, v):

        self.x = x
        self.y = y
        self.u = u
        self.v = v

    def __add__(self, other):
        return State(self.x + other.x, self.y + other.y,
                     self.u + other.u, self.v + other.v)

    def __sub__(self, other):
        return State(self.x - other.x, self.y - other.y,
                     self.u - other.u, self.v - other.v)

    def __mul__(self, other):
        return State(other * self.x, other * self.y,
                     other * self.u, other * self.v)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __str__(self):
        return f"{self.x:10.6f} {self.y:10.6f} {self.u:10.6f} {self.v:10.6f}"


class ThreeBody(object):
    """ model the gravitational interaction of three stars """

    def __init__(self, M0, pos0, M1, pos1, M2, pos2, SMALL=1.e-30):
        """ pos is a tuple of x-y coords of the stars """

        self.G = 1  # dimensionless units

        star0 = State(pos0[0], pos0[1], 0.0, 0.0)
        star1 = State(pos1[0], pos1[1], 0.0, 0.0)
        star2 = State(pos2[0], pos2[1], 0.0, 0.0)

        # this history will be stored in a list of tuples, such that
        # stars[n][istar] is star istar and timestep n

        self.stars = [(star0, star1, star2)]

        self.M = np.array([M0, M1, M2])

        self.SMALL = SMALL

        self.time = [0.0]

    def npts(self):
        """ return the number of integration points """
        return len(self.time)

    def center_of_mass(self, n):
        """return the center of mass of the system at time n"""

        Mtot = np.sum(self.M)

        x_cm = 0.0
        y_cm = 0.0

        for mass, star in zip(self.M, self.stars[n]):
            x_cm += mass * star.x
            y_cm += mass * star.y

        x_cm /= Mtot
        y_cm /= Mtot

        return x_cm, y_cm

    def energy(self, n):
        """return the energy of the system at time n """

        # kinetic energy
        KE = [0.5 * self.M[i] * (self.stars[n][i].u**2 +
                                 self.stars[n][i].v**2)
              for i in range(len(self.M))]

        # potential energy -- we need all pairs, but order doesn't matter

        PE = [-self.G * self.M[i] * self.M[j] /
              np.sqrt((self.stars[n][i].x - self.stars[n][j].x)**2 +
                      (self.stars[n][i].y - self.stars[n][j].y)**2)
              for i in range(len(self.M)) for j in range(len(self.M)) if i > j]

        return sum(KE) + sum(PE)

    def rhs(self, star_states):
        """ here, star states is (State, State, ...) """

        ydots = []

        for istar in range(len(self.M)):

            # compute the space derivatives

            dxdt = star_states[istar].u
            dydt = star_states[istar].v

            # compute the velocity derivatives
            dudt = 0.0
            dvdt = 0.0

            for jstar in range(len(self.M)):
                if istar == jstar:
                    continue
                else:
                    dx = star_states[jstar].x - star_states[istar].x
                    dy = star_states[jstar].y - star_states[istar].y

                    r = np.sqrt(dx**2 + dy**2) + self.SMALL

                    dudt += self.G * self.M[jstar] * dx / r**3
                    dvdt += self.G * self.M[jstar] * dy / r**3

            ydots.append(State(dxdt, dydt, dudt, dvdt))

        return ydots


    def single_step(self, state_old, dt):
        """ take a single RK-4 timestep from t to t+dt for the system
        ydot = rhs """

        ydot1 = self.rhs(state_old)

        state_temp = [state_old[istar] + 0.5 * dt * ydot1[istar]
                      for istar in range(len(self.M))]
        ydot2 = self.rhs(state_temp)

        state_temp = [state_old[istar] + 0.5 * dt * ydot2[istar]
                      for istar in range(len(self.M))]
        ydot3 = self.rhs(state_temp)

        state_temp = [state_old[istar] + dt * ydot3[istar]
                      for istar in range(len(self.M))]
        ydot4 = self.rhs(state_temp)

        state_new = [state_old[istar] + (dt/6.0) * (ydot1[istar] + 2.0*ydot2[istar] +
                                                    2.0*ydot3[istar] + ydot4[istar])
                     for istar in range(len(self.M))]

        return state_new

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

        # start with the old timestep
        dt_new = dt

        n_reset = 0

        while t < tmax:

            state_old = self.stars[-1]

            # adaptive stepping iteration loop -- keep trying to take
            # a step until we achieve our desired error
            rel_error = 1.e10

            n_try = 0
            while rel_error > err:
                dt = dt_new
                if t + dt > tmax:
                    dt = tmax-t

                # take 2 half steps
                state_tmp = self.single_step(state_old, 0.5*dt)
                state_new = self.single_step(state_tmp, 0.5*dt)

                # now take just a single step to cover dt
                state_single = self.single_step(state_old, dt)

                # state_new should be more accurate than state_single,
                # since it used smaller steps.

                # estimate the relative error now
                pos_err = max([max(abs((state_new[istar].x - state_single[istar].x) / state_single[istar].x),
                                   abs((state_new[istar].y - state_single[istar].y) / state_single[istar].y))
                               for istar in range(len(self.M))])

                vel_err = max([max((abs(state_new[istar].u - state_single[istar].u) / state_single[istar].u),
                                   (abs(state_new[istar].v - state_single[istar].v) / state_single[istar].v))
                               for istar in range(len(self.M))])

                rel_error = max(pos_err, vel_err)

                # adaptive timestep algorithm from Garcia (Eqs. 3.30
                # and 3.31)
                dt_est = dt * abs(err / rel_error)**0.2
                dt_new = min(max(S1*dt_est, dt/S2), S2*dt)

                n_try += 1

            if n_try > 1:
                # n_try = 1 if we took only a single try at the step
                n_reset += (n_try-1)

            # successful step
            t += dt

            # set for the next step
            self.time.append(t)
            self.stars.append(state_new)

            print(t, end="")
            for s in state_new:
                print(s, end="")
            print("\n")


if __name__ == "__main__":

    tb = ThreeBody(150., (3, 1),
                   200., (-1, -2),
                   250., (-1, 1))

    tb.integrate(0.05, 1.e-8, 2.0)

    print(f"number of steps: {tb.npts()}")

    fig = plt.figure()
    ax = fig.add_subplot()

    for istar in range(len(tb.M)):
        x = [star[istar].x for star in tb.stars]
        y = [star[istar].y for star in tb.stars]

        ax.plot(x, y)
        ax.scatter([x[0]], [y[0]], marker="x")

    fig.savefig("three_body.png")

    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot()

    t = tb.time
    x_cm = []
    y_cm = []
    for n in range(tb.npts()):
        _xcm, _ycm = tb.center_of_mass(n)
        x_cm.append(_xcm)
        y_cm.append(_ycm)

    ax.plot(t, x_cm, label="x center of mass")
    ax.plot(t, y_cm, label="y center of mass")

    fig.savefig("center_of_mass.png")

    plt.close(fig)

    # timestep plot

    t = np.asarray(t)
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

    E = []
    for n in range(tb.npts()):
        E.append(tb.energy(n))

    ax.plot(t, np.abs(E))

    ax.set_yscale("log")

    ax.set_xlabel("t")
    ax.set_ylabel("|E|")

    fig.tight_layout()
    fig.savefig("three_body_energy.png")
