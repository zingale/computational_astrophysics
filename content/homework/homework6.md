# Homework 6

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

## 1. Cepheid $P$-$L$ relation

The paper [Calibrating the Cepheid period-luminosity relation from the infrared surface brightness technique. I. The p-factor, the Milky Way relations, and a universal K-band relation](https://ui.adsabs.harvard.edu/abs/2011A%26A...534A..94S/abstract) (Storm et al. 2011) presents a table of [Cepheid variable stars](https://en.wikipedia.org/wiki/Cepheid_variable).  Cepheid's are useful in astronomy because there is a period-luminosity ($P$-$L$) relationship that allows them to be used as distance indicators.

The file [cepheids.txt](cepheids.txt) contains the data from Table 7 of that paper,
with a few stars commented out, because they were disregarded in the fit.  We will
concern ourselves with 3 columns:

* $\log P$ : the base-10 log of the period (in days)
* $M_k$ : the K-band absolute magnitude.
* $\sigma(m-M)$ : the uncertainty in the absolute magnitude due to the distance
  uncertainty.  We'll take this to be the uncertainty in $M_k$.

```{warning}
There are a few NaNs in the data file that represent missing data.  Depending on how you
read in the file you will need to deal with these.
```

```{admonition} Your task
* Use the generalized linear least squares method we did [in class](https://zingale.github.io/computational_astrophysics/fitting/generalized-linear-least-squares-tests.htm) to fit the data to a line of the form:

  $$M = a (\log_{10}(P) - 1.0) + b$$

  (this is what is used for Table 8 in Storm et al. 2011).

  Find the values of $a$ and $b$.

* Plot the data together with the fit you found.

* Using the matrix ${\bf A}$ from the generalized linear least squares, we can
  estimate the errors in the fit parameters as:

  $${\bf C} = ({\bf A}^\intercal {\bf A})^{-1}$$

  and then the uncertainty on fit parameter $j$ is:

  $$\sigma_j = \sqrt{C_{jj}}$$

  This will require that you compute the inverse of ${\bf A}^\intercal {\bf A}$:

  * If you are using python, you can use `np.linalg.inv()`.

  * If you are using C++, I've added a function in `inverse.H` to the regression code
    from class that contains a function to compute the inverse.  See:
    https://github.com/zingale/computational_astrophysics/tree/main/examples/fitting

```

## 2. Tabulating Reaction Rates

Previously we looked at the [temperature sensitivity of the 3-$\alpha$ rate](https://zingale.github.io/computational_astrophysics/basics/diff-int/application-rate-temperature-sensitivity.html).  Now we want to tabulate this rate at a small number of temperature points and use cubic interpolation to approximate the rate.

The temperature-dependent portion of the rate takes the form:

$$\lambda(T) = 5.09\times 10^{11} T_8^{-3} e^{-44.027/T8}$$

where $T_8 = T/(10^8~\mathrm{K})$.

Imagine that we tabulate this at $N$ temperature points (equally spaced in log).
When we did Lagrange interpolation, we tried to fit a polynomial to all $N$ points.

Now instead, we will find the 4 points that are closest to our desired (log) temperature, $\log(T_0)$, 
and fit a cubic polynomial to those and evaluate the rate from this interpolant.

```{tip}
Usually this will be 2 on each side of $\log T_0$, but near the ends of
our data sampling, we might need to use 1 point on one side of $\log T_0$
and 3 on the other.
```

```{admonition} Your task
* Sample $\log\lambda$ at 10 temperatures equally spaced in log between
  $10^8~\mathrm{K}$ and $5\times 10^9~\mathrm{K}$.  This will give you 10
  pairs of the form $((\log T)_i, (\log\lambda)_i)$.  Note: $\log()$ here
  means base-10 logarithm.

* Now randomly pick 25 temperature values in the same range
  $[10^8~\mathrm{K}, 5\times 10^9~\mathrm{K}]$.  For each of these points,
  do the following:

  * Call the random temperature we are working with $T_0$.  This is where we
    want to approximate the value of $\lambda$ via interpolation.
    
    For this $T_0$, find the 4 points that are closest to $\log T_0$.  Call
    these $j$, $j+1$, $j+2$, $j+3$.  We will use the data at those 4 points
    to construct our cubic interpolation.

  * Construct the Lagrange polynomial passing through those 4 data points,
    $((\log T)_j, (\log\lambda)_j), \ldots, ((\log T)_{j+3}, (\log\lambda)_{j+3})$.
    and evaluate $\log\lambda_0$ corresponding to $\log T_0$.
  
  * Now compute the value of $\lambda_0$ from the original expression and compute
    the error 

* Report the largest error you observe for all 25 points you tried.