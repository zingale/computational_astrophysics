# Homework 2

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```


1. In class we showed how Simpsons rule is derived from fitting a
   parabola to 3 points and integrating under the parabola.  Our
   composite method worked for an even number of intervals.  Now we'll
   consider the case where there is an odd number of intervals.

   If $N$ is odd, then we can use the standard composite Simpsons
   integration method for the first $N-1$ intervals (treating them in
   pairs).  For the very last interval, we need a new form of Simpsons
   rule.  Denote the last 3 points in the domain as $x_{N-2}$,
   $x_{N-1}$, $x_N$&mdash;we want to do the final integration over $x
   \in [x_{N-1}, x_N]$.

   a. Construct the form of the integral in this case by using the
      same parabolic interpolant, and evalute it by integrating only
      over $[x_{N-1}, x_N]$.

   b. Write a general composite Simpsons rule integrator that works for
      even or odd $N$ (you can start with the code from class).

   c. Integrate
        
      $$I = \int_0^5 x \sin(2\pi x) dx$$

      from $[0, 5]$ using $N = \{3, 5, 9, 17, 33\}$ and compute the error compared
      to the analytic solution,

      $$I_\mathrm{analytic} = -\frac{5}{2\pi}$$

      How does the error converge?

2. The particles in an ideal gas a well-described by the Maxwell-Boltzmann distribution:

   $$n(p) d^3 p \rightarrow 4\pi n(p) p^2 dp = \frac{n_I}{(2\pi m_I k_B T)^{3/2}} e^{-p^2/(2m_I k_BT)} 4\pi p^2 dp$$

   The average velocity of this distribution, defined as:

   $$\langle v \rangle = \frac{1}{n_I} \int_0^\infty 4\pi n(p) p^2 \left (\frac{p}{m_I}\right ) dp$$

   Evaluate this integral numerically for a gas of protons ($m_I = m_p)$ and the conditions
   in the center of the Sun ($T = 1.5\times 10^7~\mathrm{K}$).

   Your integral should agree with the analytic result:

   $$\langle v \rangle = 2 \sqrt{\frac{2k_B T}{\pi m_I}}$$

   Be sure to vary the number of intervals in your numerical
   integration to show that you get an accurate answer (i.e., you
   converge).
