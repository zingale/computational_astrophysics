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

  This will require that you compute the inverse of ${\bf A}^\intercal {\bf A}$.
  If you are using python, you can use `np.linalg.inv()`.  If you are using a
  different language, I can provide you with a function if you let me know.

  
```

## 2. Tabulating Reaction Rates

