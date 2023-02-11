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
   $x_{N-1}$, $x_N$, and we want to do the final integration over $x
   \in [x_{N-1}, x_N]$.  Construct the form of the integral in this
   case by using the same parabolic interpolant, but now only
   integrating over the right half of it.
   
   Write a general composite Simpsons rule integrator that works for
   even or odd $N$ (you can start with the code from class).
   
   
2. Integrate Maxwellian?

3. Integrate reaction rate <sigma v>?

