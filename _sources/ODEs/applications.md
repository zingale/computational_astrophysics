Application: Few Body Problem
=============================

Consider a system of masses (e.g., stars) interacting gravitationally.  Each
experiences a gravitational force from the others.  For $N$ masses,
there are $N(N-1)/2$ interactions, which means that the work involved
in evaluating all of the forces scales like $\mathcal{O}(N^2)$ for
large $N$.

We'll solve the "few body" problem under the simplifying assumption
that all of the motion is in the $x-y$ plane.

The evolution equations for the system are:

$$\frac{d^2 {\bf x}_i}{dt^2} = \sum_{j = 1; j \ne i}^{N}
                G m_j \frac{{\bf x}_j - {\bf x}_i}{|{\bf x}_j - {\bf x}_i|^3}$$


where $i = 1, 2, \ldots, N$, and ${\bf x}_i = (x_i, y_i)$ are the coordinates
of star $i$.  We can write this as a system of first-order ODEs by
introducing the velocity,

$${\bf v}_i = \frac{d{\bf x}_i}{dt}$$

There is one additional bit&mdash;if two stars get too close, then the
denominator can become very small, and this will greatly impact our
timestep.  The standard way of dealing with this is to introduce a
[softening
length](https://en.wikipedia.org/wiki/N-body_simulation#Softening),
changing the equations to:

\begin{align*}
\frac{d{\bf x}_i}{dt} &= {\bf v}_i \\
\frac{d{\bf v}_i}{dt} &= \sum_{j = 1; j \ne i}^{N}
                G m_j \frac{{\bf x}_j - {\bf x}_i}{(|{\bf x}_j - {\bf x}_i|^2 + \epsilon^2)^{3/2}}
\end{align*}

where $\epsilon$ is chosen to be small.


Implementation
--------------

Here's an implementation for $N = 3$.  Extending it to more stars is
straightforward.  We work in units with $G = 1$.

```{note}
You should try extending this to an arbitrary $N$.  You'll want to change the
interface to allow you to pass in a list of initial particles, or create a
method that randomly initializes them
```

We'll solve this with an adaptive RK4 method.  The code implementing
this is shown below:

````{tab-set}
```{tab-item} python
```{literalinclude} ../../examples/few_body/three_body.py
:language: python
```

```{tab-item} C++
```{literalinclude} ../../examples/few_body/three_body.cpp
:language: c++
```
````


```{note}
If we start with zero velocity, then the center of mass is at rest, and we can
compute the center of mass of the system at all times to see how well it is
preserved (it should not change, since there is no outside force acting on the system).

Likewise, we can look at the total energy of the system as a diagnostic.
How well is the energy conserved in this calculation?
```


Going Further
-------------

For a very large number of particles, the $N^2$ work becomes very
computationally expensive.  The [Barnes-Hut tree
algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)
reduces the work to $\mathcal{O}(N\log N)$.
