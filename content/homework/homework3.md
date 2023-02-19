# Homework 3

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

1. Consider a simple pendulum.  The equations of motion
   are:

   $$\dot{\theta} = \omega$$
   $$\dot{\omega} = -\frac{g}{L} \sin \theta$$

   where $\theta$ is the angular displacement from vertical and
   $\omega$ is the angular velocity.  The angular acceleration in this
   case is $\alpha = -(g/L) \sin\theta$.  In your intro mechanics
   class, you probably solved this using the small-angle approximation
   an found the solution has the form:
  
   $$\theta(t) = A \cos(\tfrac{2\pi}{T} t + \phi)$$
  
   where the period (in the small-angle approximation) is:
  
   $$T = 2\pi \sqrt{\frac{L}{g}}$$

   We will solve this system without making the small-angle
   approximation.  For larger amplitude displacements, the period takes
   the form:
  
   $$T = 2\pi \sqrt{\frac{L}{g}} \left ( 1 + \frac{1}{16}\theta_m^2 + \ldots \right )$$
   
   (see, e.g, [Wikipedia](https://en.wikipedia.org/wiki/Pendulum_(mechanics)#Power_series_solution_for_the_elliptic_integral)), where $\theta_0$ is the initial angular displacement of the pendulum.
  
   The total energy of the system is
  
   $$E = \frac{1}{2} m L^2 \omega^2  - mg L \cos \theta$$

   * Solve the pendulum system using both the Euler and
     Euler-Cromer methods.  Pick $L = 10$ m and $g = 10$ m/s$^2$.
     Compare the solutions for $\theta = 10^\circ$ and $100^\circ$.
     Make a plot of both $\theta(t)$ vs. $t$ and $E$ vs. $t$.  Notice
     the stark difference in the behavior between the two methods.

     In each case, estimate the period from your numerical solution.

   * Now consider the velocity Verlet method 
     (_kick-drift-kick_) we discussed in class.  As we
     saw in class, this applies when the acceleration term ($\alpha$ in
     our case) does not depend on the velocity ($\omega$ for us).

     Implement this method for the simple pendulum and estimate its
     convergence by comparing the total energy after several periods to
     the initial energy for a variety of choices of $\tau$.

  Use a fixed timestep for all of these integrations.


