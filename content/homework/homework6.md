# Homework 6

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

1. Let's consider two different first-order finite-difference methods for Burgers'
   equation:

   $$u_t + u u_x = 0$$

   As we saw in class, this is a nonlinear equation which can be written in conservative
   form as:

   $$u_t + \left [ \frac{1}{2} u^2 \right ]_x = 0$$

   Although these forms are analytically equivalent, the numerical solutions
   will differ.


   * Using first-order upwinding, we can difference the first form of Burgers' equation
     as:

     $$u_{i}^{n+1} = u_i^n - \frac{\Delta t}{\Delta x} u_i^n (u_i^n - u_{i-1}^n)$$

     This is valid as long as $u_i^n > 0$.  Write a program to solve this
     equation with initial conditions:

     $$u(x,t=0) = \left \{ \begin{array}{cc} 2 & \mathrm{if~} x < 0.5 \\
                                        1 & \mathrm{if~} x \ge 0.5 \\
                      \end{array} \right .$$

     These are the initial conditions for a shock that we explored in class.

     Important: This is a non-linear equation---make sure you use the
     proper CFL condition for this method.  This means: pick a value for $C$,
     then evaluate $\Delta t = C \Delta x / \max_i\{u\}$

   * Now difference the conservative form as:

     $$u_{i}^{n+1} = u_i^n - \frac{\Delta t}{\Delta x} \left ( \frac{1}{2} (u_i^n)^2
         - \frac{1}{2} (u_{i-1}^n)^2 \right )$$

     Write a program to solve this, and run it on the same initial data as
     in before.  Note that both of these discretizations assume that
     $u > 0$.

   * How do the solutions differ for the two discretizations?  Measure the speed of the shock for each
     method.

     How does the shock speed differ with resolution?

     Which method gets the shock speed correct?  Why?

     Present your results at 3 different grid resolutions.
