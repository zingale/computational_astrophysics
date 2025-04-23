# Homework 9

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

```{important}
Choose **one** of the problems below.  You do not need to do both!!!
```

## Option 1: Diffusion

Let's redo the implicit diffusion solve using a second-order
accurate discretization in time: the [Crank-Nicolson
method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method).
The discretized diffusion equation with this method appears as:

   $$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta t} =
           \frac{k}{2} \left ( \frac{\phi_{i+1}^n - 2 \phi_i^n + \phi_{i-1}^n}{\Delta x^2} +
                               \frac{\phi_{i+1}^{n+1} - 2 \phi_i^{n+1} + \phi_{i-1}^{n+1}}{\Delta x^2} \right )$$

   The main difference here with what we worked out in class is that the righthand side
   is now centered in time.  If we define:

   $$\alpha \equiv \frac{k\Delta t}{\Delta x^2}$$

   and put all of the "$n+1$" terms on one side, we have:

   $$-\frac{\alpha}{2} \phi_{i+1}^{n+1} + (1 + \alpha) \phi_i^{n+1} - \frac{\alpha}{2} \phi_{i-1}^{n+1} =
          \phi_i^n + \frac{\alpha}{2} (\phi_{i+1}^n - 2 \phi_i^n + \phi_{i-1}^n)$$

   We'll write the righthand side as $f_i$:

   $$f_i = \phi_i^n + \frac{\alpha}{2} (\phi_{i+1}^n - 2 \phi_i^n + \phi_{i-1}^n)$$

   Now, just like in class, we need to see what happens at the boundary.  For the left
   most point at the boundary, $\mathrm{lo}$, with Neumann boundary conditions,
   we have $\phi_\mathrm{lo-1} = \phi_\mathrm{lo}$ and the expression becomes:

   $$\left ( 1 + \frac{\alpha}{2} \right ) \phi_\mathrm{lo}^{n+1} - \frac{\alpha}{2} \phi_\mathrm{lo+1}^{n+1} = f_\mathrm{lo}$$

   and similarly at the right boundary.  With this, our linear system takes the form:

   $$\left (
   \begin{array}{ccccccc}
   1+\tfrac{\alpha}{2} &   -\tfrac{\alpha}{2} &           &        &         &           &          \\
   -\tfrac{\alpha}{2}  & 1+\alpha & -\tfrac{\alpha}{2}   &        &         &           &          \\
            & -\tfrac{\alpha}{2}   & 1+\alpha & -\tfrac{\alpha}{2}&         &           &          \\
            &           & \ddots    & \ddots & \ddots  &           &          \\
            &           &           & \ddots & \ddots  & \ddots    &          \\
            &           &           &        & -\tfrac{\alpha}{2} & 1+\alpha &-\tfrac{\alpha}{2}   \\
            &           &           &        &         & -\tfrac{\alpha}{2}   &1+\tfrac{\alpha}{2} \\
   \end{array}\right )
   \left ( \begin{array}{c}
   \phi_\mathrm{lo}^{n+1} \\
   \phi_\mathrm{lo+1}^{n+1} \\
   \phi_\mathrm{lo+2}^{n+1} \\
   \vdots \\
   \vdots \\
   \phi_\mathrm{hi-1}^{n+1} \\
   \phi_\mathrm{hi}^{n+1} \\
   \end{array} \right )
   =
   \left (\begin{array}{c}
   f_\mathrm{lo}^{n} \\
   f_\mathrm{lo+1}^{n} \\
   f_\mathrm{lo+2}^{n} \\
   \vdots \\
   \vdots \\
   f_\mathrm{hi-1}^{n} \\
   f_\mathrm{hi}^{n} \\
   \end{array} \right )
   $$

   ```{admonition} Your task
   Solve this system by modifying the direct solve example from class
   (or implementing your own version).  Alternately, you can choose to do it via
   relaxation, again following what we did in class as a guide.

   Vary the resolution, and demonstrate that this method converges better than the
   first-order accurate (in time) method that we implemented in class.
   ```

## Option 2: FFTs for Poisson

Consider the 1d Poisson equation:

$$\phi^{\prime\prime} = 2 (1 - 3 x)$$

on $[0, 1]$ with periodic boundary conditions.  This
has the solution:

$$\phi(x) = x^2 (1 - x)$$

We'll call the righthand side of the equation $f$,

$$f = 2 (1 - 3 x)$$

Now, we can express $\phi$ and $f$ in terms of their Fourier transforms:

$$\phi(x) = \int \Phi(k) e^{-2\pi i kx} dk$$

$$f(x) = \int F(k) e^{-2\pi i k x} dk$$


Note that:

$$\frac{d^2 \phi(x)}{dx^2} = - 4 \pi^2 k^2 \int \Phi(k) e^{2\pi i kx}dk$$

and since the Poisson equation is linear, the different
modes will not mix.  

Putting these into our Poisson equation we have:

$$-4\pi^2 k^2 \Phi(k) = F(k)$$

This allows us to easily solve for the Fourier transform
of $\phi(x)$, $\Phi(k)$, algebraically:

$$\Phi(k) = - \frac{F(k)}{4\pi^2 k^2}$$

We can then transform back to get the solution to the
original Poisson equation in real space.

```{admonition} Your task
Write a code to solve our Poisson equation using Fourier
transforms.  You will do the following:

* Compute the FFT of $f(x)$
* Compute the Fourier components of $\Phi(k)$ using
  our algebraic expreession: $\Phi(k) = - {F(k)}/(4\pi^2 k^2)$

* Compute the *inverse* FFT of $\Phi(k)$ to get $\phi(x)$
* Compute the error with respect to the exact solution.

Do this for several grid sizes, $N = 32, 64, 128$ and
measure the convergence.

```



