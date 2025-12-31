# Advection

In astrophysics, many systems (stars, the ISM, gas in galaxies) is modeled as
a fluid, following the equations of hydrodynamics:

\begin{align*}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho {\bf U}) &= 0 \\
\frac{\partial (\rho {\bf U})}{\partial t} + \nabla \cdot (\rho {\bf U}{\bf U}) + \nabla p 
&= 0 \\
\frac{\partial (\rho E)}{\partial t} + \nabla \cdot (\rho {\bf U} E + {\bf U}p) &= 0
\end{align*}

These express conservation of mass, momentum, and energy.  If we expand out the diveragence
in the first equation, we get:

$$\frac{\partial \rho}{\partial t} + {\bf U} \cdot \nabla \rho + \rho \nabla \cdot {\bf U} = 0 $$

This has the form of an advection equation.  

```{note}
Our model equation for this, in one-dimension
will be the _linear advection equation_:

$$\frac{\partial a}{\partial t} + u \frac{\partial a}{\partial x} = 0$$
```

Looking at the second equation (momentum), and expanding out and simplifying, we have:

$$\frac{\partial {\bf U}}{\partial t} + {\bf U} \cdot \nabla {\bf U} + \frac{1}{\rho} \nabla p = 0$$

This is a nonlinear equation (because of the ${\bf U}\cdot \nabla {\bf
U}$ term).  This nonlinearity admits a rich diversity of solutions:
rarefactions, shocks, and turbulence.

```{note} 
Our model equation for this will be the [inviscid Burgers'
equation](https://en.wikipedia.org/wiki/Burgers%27_equation):

$$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = 0$$
```
