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

Let's get some practice with git.  Here is a git repository
on github that we looked at in class.

https://github.com/ast390-sbu/test-repo

You should make sure that you have a github account.


```{admonition} Your task
* *Fork* this repository: https://github.com/ast390-sbu/test-repo

  This will create a version of the git repo under your account.

* *Clone* your fork.

* Create a change to the repo, for example, add your
  name to the `README.md` file or add a file of your own
  to the repo.  Whatever you wish.

* *add* and *commit* your change and *push* your change back
  to your fork.

* From the github page for your fork, create a *pull request*
  to the original repository.
```

```{tip}
We are essentially doing this sequence from class:
https://zingale.github.io/computational_astrophysics/git/pull-requests.html

If you did it in class while following along then you are done.
```
