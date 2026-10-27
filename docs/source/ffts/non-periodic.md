# Non-uniform Data

The FFT assumes that the input data is uniformly spaced.  But if you are working with
a collection of observational data, it may be sampled at irregular time intervals.
In this case, the [Lomb-Scargle periodogram](https://en.wikipedia.org/wiki/Least-squares_spectral_analysis#The_generalized_Lomb%E2%80%93Scargle_periodogram) is often used.

```{tip}
A nice introduction to the Lomb-Scargle periodogram is provided in the
paper [Understanding the Lomb-Scargle Periodogram](https://arxiv.org/abs/1703.09824) by
J. VanderPlas.
```

SciPy provides the [`scipy.signal.lombscargle` function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lombscargle.html)
which implements this algorithm.
