# Course Project

```{note}
The final project counts as 30% of your total course grade.  You are expected
to work independently on the project.

The project is due on the last day of class, Thurs. May 8.  There will be _no_ extensions,
since I need to turn in grades shortly afterwards.
```

For your final course project, you are to pick an interesting
astrophysical problem or algorithm and explore it in some depth, to
demonstrate you understand how the numerical method works and how it
is used in research.  You should pick a topic that interests
you&mdash;some suggestions are given below but feel free to pick
something else.

If you are implementing the algorithm yourself, then you will turn in
your own code.  If you are using an existing code or library, you will
give detailed instructions on what you did, together with any code
additions that were necessary.

```{important}
On the last week of class, each student will git a 3 minute / 1-2 slide summary
of their project to the class
```

## You will turn in:

1. You code

   * This can be a Jupyter notebook, C++ code, etc. if your wrote your
     own implementation of a method.

     Make sure I can build your code easily.

   * This may also be detailed instructions on how to obtain, modify,
     and reproduce your results using an external code if you chose to
     use that.

2. A PDF write-up (abstract, introduction, description of what you did
   / results, figures, references).  This is something that you can use
   for the upper-level writing requirement, if desired.

   This should include:

   * An Explanation of how the algorithm works, referring to the
     topics we learned in class.

   * A demonstration of convergence (if applicable) or controlling
     the error.

   * Some example tests where we know the right behavior (if applicable)

   You should aim for 10 pages (including figures and references).
   
3. 1-2 slides for the in-class presentation (details will be shared
   toward the end of the semester).

## Some inspiration for projects:

* Integrate the pressure and energy for an electron gas (building on
  the example we did for the degeneracy parameter), and build a table
  that you can interpolate from.

* Implement the Barnes-Hut tree algorithm for approximating N-body gravity.
  https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation

  Demonstrate that it works by comparing to the direct sum, $\mathcal{O}(N^2)$
  algorithm.

* Any of the ODE example applications listed here:
  https://zingale.github.io/computational_astrophysics/ODEs/more-applications.html

* Extending our linear advection solver to two-dimensions.

* Reproduce the analysis of any astro paper that has publicly available data.

* Download an existing simulation code and run some simulations
  (existing examples or your own sets of initial conditions).  Some codes that you
  might read about include:

  * REBOUND: an N-body integrator with a python interface:
    https://rebound.readthedocs.io/en/latest/

    Some example applications:

    * Reproduce some of the results from [_Dynamics of Systems of Two
      Close Planets_](https://ui.adsabs.harvard.edu/abs/1993Icar..106..247G/abstract) by B. Gladman

  * pyro: a python hydrodynamics code written by my research group:
    https://python-hydro.github.io/pyro2/

  * MESA: a stellar evolution code:
    https://docs.mesastar.org/en/release-r22.11.1/

* Read about the _smoothed particle hydrodynamics_ method (SPH) for
  solving PDEs and implement a simple advection solver using this.

* Read about some astrophysical applications of machine learning and
  try to reproduce their results or do some simple classification of
  astrophysical data on your own.

## Rubric for grading

Project scores will be assigned with the following weighting (out of 30)

* clearly commented, working code the produces the results in your
  write-up _or_ clear instructions + additional custom code needed to
  use an external code / library to reproduce your results.  (12 pts.)

* a thorough write-up that demonstrates that you understand:

  * the basic ideas of the algorithm you are working with
  * the astrophysical applications of your algorithm
  * the connections to the topics we covered in class
  * the limitations of the method

  and also indicates that you spent some time exploring the method in-depth.

  overall (15 pts.)

* an in-class summary of your project (3 pts.)
