# Homework 4

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

Equations of state (EOSs) are often formulated in terms of density, $\rho$,
and temperature, $T$.  But it is common to need to use the equation
of state with different inputs (like pressure, $p$) and recover the density 
or temperature that yields thermodynamic consistency.  This is a root-finding problem (sometimes we say that we need to *invert* the EOS).

Here we will consider an ideal gas and radiation EOS:

\begin{align*}
p &= \frac{\rho k T}{\mu m_u} + \frac{1}{3} a T^4 \\
e &= \frac{3}{2} \frac{k T}{\mu m_u} + \frac{a T^4}{\rho}
\end{align*}

where $e$ is the specific internal energy (energy / mass)
and the constants are:

* k: Boltzmann's constant ($k = 1.38\times 10^{-16}~\mathrm{erg/K}$)
* $\mu$: the mean molecular weight of the matter (we'll just take $\mu = 4$)
* $m_u$: the atomic mass unit ($m_u = 1.66\times 10^{-24}~\mathrm{g}$)
* $a$: the radiation constant ($a = 7.56\times 10^{-15}~\mathrm{erg/cm^3/K^4}$).

Imagine that we know the value of $p$, and $e$---we'll call them $p_\star$ and $e_\star$.  We can find the density and temperature that give these by solving the system:

\begin{align*}
p(\rho, T) &= p_\star \\
e(\rho, T) &= e_\star
\end{align*}

Let's define $\Psi$ as:

$$\Psi = \left ( \begin{array}{c} p(\rho, T) - p_\star \\
                                  e(\rho, T) - e_\star \end{array} \right )$$

and imagine that we have an initial guess for density, $\rho_0$, and temperature, $T_0$.  We then want to find the corrections, $\delta \rho$
and $\delta T$, such that:

$$\Psi(\rho_0 + \delta\rho, T_0 + \delta T) = 0$$

Tayloring expanding, we get:

$$\Psi(\rho_0 + \delta_\rho, T_0 + \delta T) = 0 = \Psi(\rho_0, T_0) + {\bf J} \left ( \begin{array}{c} \delta \rho \\ \delta T \end{array} \right ) + \ldots$$

where ${\bf J}$ is the Jacobian,

$${\bf J} = \left ( \begin{array}{cc} \partial p/\partial \rho |_T & \partial p/\partial T |_\rho \\ 
                                      \partial e/\partial \rho |_T & \partial e/\partial T |_\rho \end{array} \right )$$

these partial derivatives can be computed analytically from the EOS.

We can then find the correction by solving the system:

$${\bf J} \left (\begin{array}{c} \delta \rho \\ \delta T \end{array} \right ) = - \Psi(\rho_0, T_0)$$

We would need to iterate until the correction is small.

```{admonition} Your task
Write a function that finds $\rho$ and $T$ given an input $p_\star$ and
$e_\star$.  The function should do the following:

* Pass in $p_\star$ and $e_\star$, as well as initial guesses for
  density and temperature, $(\rho_0, T_0)$

* Loop until the size of the corrections $\delta\rho$ and $\delta T$
  are small:

  * Call the equation of state and compute $p_0 = p(\rho_0, T_0)$ and $e_0 = e(\rho_0, T_0)$ as well as the 4 derivatives needed in the Jacobian

  * Construct the Jacobian, ${\bf J}$
  * Construct $\Psi_0 = (p_0 - p_\star, e_0 - e_\star)$
  * Solve the linear system:

    $${\bf J} \left (\begin{array}{c} \delta \rho \\ \delta T \end{array} \right ) = - \Psi_0$$

    (you can use any method you want---one of the functions
    from class, a library, etc.)
    
  * Correct your initial guesses:

    \begin{align*}
    \rho_0 &\rightarrow \rho_0 + \delta \rho \\
    T_0 &\rightarrow T_0 + \delta T
    \end{align*}

Test this out by finding the $(\rho, T)$ that gives $p_\star = 2.3\times 10^{10}~\mathrm{erg/cm^3}$
and $e_\star = 3.87\times 10^{13}~\mathrm{erg/g}$.

To check your answer, use the $(\rho, T)$ you found to
evaluate $p$ and $e$ via the EOS, and see that they
should be close to $p_\star$ and $e_\star$.
```
                                      

                                  
