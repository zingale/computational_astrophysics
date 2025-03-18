# Homework 5

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

## 1. Pendulum

In [homework 3](https://zingale.github.io/computational_astrophysics/homework/homework3.html)
you integrated a pendulum without using the small-angle approximation.

Now let's take the FFT!

```{admonition} Your task
* Recreate solution for the $\theta_0 = 100^\circ$ evolution,
  $\theta(t)$, by integrating the system for *1 period* (you can
  use the expression for the finite-amplitude period or just the
  small-angle value, $T = 2\pi$).

  This will result in a series of points
  in the time-domain that we'll call $\theta_n$.

* Compute the FFT of $\theta(t)$,

  $$\Theta_k = \mathcal{F}(\theta_n)$$

  You can use a library (like `numpy.fft.rfft`) or my DFT code from class (I
  provided a python and C++ implementation).
  
* Plot $|\Theta_k|$ vs. $\nu_k$, where $\nu_k$ is the frequency---is the power
  where you expect it should be?
  
* Repeat this for the same problem, but now integrating for 10 periods---how does
  the FFT change?
```

You are free to use my homework 3 solutions for the integration if you had difficulty
with the assignment originally.

## 2. Git practice
