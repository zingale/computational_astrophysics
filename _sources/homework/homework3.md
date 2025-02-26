# Homework 3

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

Consider a simple pendulum.  The equations of motion
are:

   $$\dot{\theta} = \omega$$
   $$\dot{\omega} = -\frac{g}{L} \sin \theta$$

   where $\theta$ is the angular displacement from vertical and
   $\omega$ is the angular velocity.  The angular acceleration in this
   case is $\alpha = -(g/L) \sin\theta$.

```{note}
In your intro mechanics
class, you probably solved this using the small-angle approximation
an found the solution has the form:

$$\theta(t) = A \cos(\tfrac{2\pi}{T} t + \phi)$$

where the period (in the small-angle approximation) is:

$$T = 2\pi \sqrt{\frac{L}{g}}$$
```

We will solve this system without making the small-angle
approximation.  For large-amplitude displacements, the period takes
the form:

$$T = 2\pi \sqrt{\frac{L}{g}} \left ( 1 + \frac{1}{16}\theta_0^2 + \ldots \right )$$

(see, e.g, [Wikipedia](https://en.wikipedia.org/wiki/Pendulum_(mechanics)#Power_series_solution_for_the_elliptic_integral)), where $\theta_0$ is the initial angular displacement of the pendulum.

The total energy of the system is

   $$E = \frac{1}{2} m L^2 \omega^2  - mg L \cos \theta$$

```{admonition} your task
Write a program to integrate this system *using fixed timesteps*
with the following methods:

* Euler-Cromer.
* Velocity-verlet (_kick-drift-kick_).  As we
  saw in class, this applies when the acceleration term ($\alpha$ in
  our case) does not depend on the velocity ($\omega$ for us).
* 4th order Runge-Kutta.

Pick $L = 10$ m and $g = 10$ m/s$^2$ for all runs.

You may use the code from class as your starting point.

1. Run twice with each solver with different initial values for
   $\theta$: $\theta_0 = 10^\circ$ and $\theta_0 = 100^\circ$.
   Choose your timestep, $\tau$, so you have at least 20 steps
   in a period (the small-angle approximate period would be $2\pi~\mathrm{s}$).
   Run for at least 10 periods.

   Make a plot of $\theta(t)$ vs. $t$ for each run (2 initial conditions $\times$ 3 solvers)
   and estimate the periods from your numerical solution.

2. Plot the total energy vs. time and comment on how the methods
   compare.

3. Now run the $\theta_0 = 100^\circ$ case with each solver, using timesteps of 
   $\tau = 2\pi/20$, $2\pi/40$, and $2\pi/80$, and plot the energy vs. time
   for each solver and comment on the convergence.
