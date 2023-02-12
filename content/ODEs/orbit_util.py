import numpy as np
import matplotlib.pyplot as plt

GM = 4*np.pi**2

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

def rhs(state):
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

def initial_conditions():
    x0 = 0
    y0 = 1
    u0 = -np.sqrt(GM / y0)
    v0 = 0

    return OrbitState(x0, y0, u0, v0)

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
