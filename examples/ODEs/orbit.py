"""
orbital motion.  We consider low mass objects orbiting the Sun.  We
work in units of AU, yr, and solar masses.  From Kepler's third law:

4 pi**2 a**3 = G M P**2

if a is in AU, P is in yr, and M is in solar masses, then

a**3 = P**2

and therefore

4 pi**2 = G

we work in coordinates with the Sun at the origin

"""

import numpy as np
import matplotlib.pyplot as plt

# global parameters
GM = 4.0 * np.pi**2  #(assuming M = 1 solar mass)


class OrbitState:
    # a container to hold the star positions
    def __init__(self, x, y, u, v):
        self.x = x
        self.y = y
        self.u = u
        self.v = v

    def __add__(self, other):
        return OrbitState(self.x + other.x, self.y + other.y,
                          self.u + other.u, self.v + other.v)

    def __sub__(self, other):
        return OrbitState(self.x - other.x, self.y - other.y,
                          self.u - other.u, self.v - other.v)

    def __mul__(self, other):
        return OrbitState(other * self.x, other * self.y,
                          other * self.u, other * self.v)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __str__(self):
        return f"{self.x:10.6f} {self.y:10.6f} {self.u:10.6f} {self.v:10.6f}"


class Orbit:
    """ hold the initial conditions of a planet/comet/etc. orbiting
        the Sun and integrate """

    def __init__(self, a, e):
        """ a = semi-major axis (in AU),
            e = eccentricity """

        self.a = a
        self.e = e

        # perihelion velocity (see C&O Eq. 2.33 for ex)
        x0 = 0.0          # start at x = 0 by definition
        y0 = a * (1.0 - e)  # start at perihelion

        u0 = -np.sqrt((GM/a)* (1.0 + e) / (1.0 - e))
        v0 = 0.0

        self.state0 = OrbitState(x0, y0, u0, v0)

    def kepler_period(self):
        """ return the period of the orbit in yr """
        return np.sqrt(self.a**3)

    def rhs(self, state):
        """ RHS of the equations of motion."""

        # current radius
        r = np.sqrt(state.x**2 + state.y**2)

        # position
        xdot = state.u
        ydot = state.v

        # velocity
        udot = -GM * state.x / r**3
        vdot = -GM * state.y / r**3

        return OrbitState(xdot, ydot, udot, vdot)

    def int_euler(self, dt, tmax):
        """ integrate the equations of motion using Euler's method.
            We integrate until t = tmax """

        # initial conditions
        times = [0.0]
        history = [self.state0]

        t = 0.0

        while t < tmax:

            state_old = history[-1]

            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t + dt > tmax:
                dt = tmax-t

            # get the RHS
            ydot = self.rhs(state_old)

            # advance
            state_new = state_old + dt * ydot
            t += dt

            # store
            times.append(t)
            history.append(state_new)

        return times, history

    def int_euler_cromer(self, dt, tmax):
        """ integrate the equations of motion using the Euler-Cromer method.
            We integrate until t = tmax """

        # initial conditions
        times = [0.0]
        history = [self.state0]

        t = 0.0

        while t < tmax:

            state_old = history[-1]

            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t + dt > tmax:
                dt = tmax-t

            # get the RHS
            ydot = self.rhs(state_old)

            # do the Euler update
            unew = state_old.u + dt * ydot.u
            vnew = state_old.v + dt * ydot.v

            xnew = state_old.x + dt * unew
            ynew = state_old.y + dt * vnew

            state_new = OrbitState(xnew, ynew, unew, vnew)
            t += dt

            # store the state
            times.append(t)
            history.append(state_new)

        return times, history

    def int_rk2(self, dt, tmax):
        """ integrate the equations of motion using 2nd order R-K
            method (also know as the midpoint method) We integrate
            until we reach tmax """

        # initial conditions
        times = [0.0]
        history = [self.state0]

        t = 0.0

        while t < tmax:

            state_old = history[-1]

            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t + dt > tmax:
                dt = tmax-t

            # get the RHS
            ydot = self.rhs(state_old)

            # advance to dt/2 to get intermediate position
            state_tmp = state_old + 0.5* dt * ydot

            # get the RHS with the intermediate state
            ydot = self.rhs(state_tmp)

            state_new = state_old + dt * ydot

            t += dt

            # store
            times.append(t)
            history.append(state_new)

        return times, history

    def int_velocity_verlet(self, dt, tmax):
        """ integrate the equations of motion using the Velocity Verlet
            method (which is a kick-drift-kick Leapfrog method) """

        # initial conditions
        times = [0.0]
        history = [self.state0]

        t = 0.0

        while t < tmax:

            state_old = history[-1]

            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t+dt > tmax:
                dt = tmax-t

            ydot_old = self.rhs(state_old)

            uhalf = state_old.u + 0.5 * dt * ydot_old.u
            vhalf = state_old.v + 0.5 * dt * ydot_old.v

            xnew = state_old.x + dt * uhalf
            ynew = state_old.y + dt * vhalf

            # this works only because the acceleration does not depend
            # on velocity
            ydot_new = self.rhs(OrbitState(xnew, ynew, uhalf, vhalf))

            unew = uhalf + 0.5 * dt * ydot_new.u
            vnew = vhalf + 0.5 * dt * ydot_new.v

            state_new = OrbitState(xnew, ynew, unew, vnew)

            t += dt

            # store
            times.append(t)
            history.append(state_new)

        return times, history

    def int_rk4(self, dt, tmax):
        """ integrate the equations of motion using 4th order R-K
            method.  """

        # initial conditions
        times = [0.0]
        history = [self.state0]

        t = 0.0

        while t < tmax:

            state_old = history[-1]

            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t+dt > tmax:
                dt = tmax-t

            # get the RHS at several points
            ydot1 = self.rhs(state_old)

            state_tmp = state_old + 0.5 * dt * ydot1
            ydot2 = self.rhs(state_tmp)

            state_tmp = state_old + 0.5 * dt * ydot2
            ydot3 = self.rhs(state_tmp)

            state_tmp = state_old + dt * ydot3
            ydot4 = self.rhs(state_tmp)

            # advance
            state_new = state_old + (dt / 6.0) * (ydot1 + 2.0 * ydot2 +
                                                  2.0 * ydot3 + ydot4)
            t += dt

            # store
            times.append(t)
            history.append(state_new)

        return times, history


def plot(history, ax=None, label=None):
    """make a plot of the solution.  If ax is None we setup a figure
    and make the entire plot returning the figure object, otherwise, we
    just append the plot to a current axis"""

    fig = None

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # draw the Sun
        ax.scatter([0], [0], marker=(20,1), color="y", s=250)

    # draw the orbit
    xs = [q.x for q in history]
    ys = [q.y for q in history]

    ax.plot(xs, ys, label=label)

    if fig is not None:
        ax.set_aspect("equal")
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")

    return fig

def error_radius(history):

    # define the error to be distance from (0, 0) at end compared to start
    R_orig = np.sqrt(history[0].x**2 + history[0].y**2)
    R_new = np.sqrt(history[-1].x**2 + history[-1].y**2)
    e = np.abs(R_new - R_orig)

    return e

def error_radius(history):

    # define the error to be distance from (0, 0) at end compared to start
    R_orig = np.sqrt(history[0].x**2 + history[0].y**2)
    R_new = np.sqrt(history[-1].x**2 + history[-1].y**2)
    e = np.abs(R_new - R_orig)

    return e

def error_position(history):
    """return the difference in the distance from the Sun"""
    dx = history[0].x - history[-1].x
    dy = history[0].y - history[-1].y
    return np.sqrt(dx**2 + dy**2)

def energy(history):
    """return a vector of energy for all times"""
    E = [0.5 * (state.u**2 + state.v**2) - GM / np.sqrt(state.x**2 + state.y**2)
         for state in history]
    return E

