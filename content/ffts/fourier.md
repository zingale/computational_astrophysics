# Fourier Transforms

The [Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform) converts a function from a physical-space (or time) representation
into a frequency-space representation.  The two representations are equivalent, but give
alternate views into how the function behaves.

The general form of a Fourier transform is:

$$F(k) = \int_{-\infty}^{\infty} f(x) e^{-2\pi i k x} dx$$

where $f(x)$ is the real-space form of the function.
If we think of this in terms of space, $x$ has dimensions of length and $k$ has dimensions of 1/length, and can be interpreted as a wavenumber.

You can think of $F(k)$ as being the amount of the function $f$ represented by a frequency $k$.

If we think of it in terms of time, then $f(t)$ is the time-domain function and the transform is:

$$F(\nu) = \int_{-\infty}^{\infty} f(t) e^{-2\pi i \nu t} dx$$

and $\nu$ is now a frequency (units of 1/time).

The inverse transform is then:

$$f(x) = \int_{-\infty}^\infty F(k) e^{2\pi ikx} dk$$

```{note}
Alternate definitions of the transform exist, including swapping the minus sign in the exponental between the forward and inverse transform.
```

## Discrete Fourier transform

The discrete Fourier transform (DFT) works on discrete data, e.g.:

* a time-series from an experiment
* simulation data for a velocity field

We'll work on the assumption that our data is evenly-spaced.

The discrete analogue of the Fourier transform is:

$$F_k = \sum_{n = 0}^{N-1} f_n e^{-2\pi i n k / N}$$

where we have $N$ samples with values $f_n$.

The inverse transform is:

$$f_n = \frac{1}{N} \sum_{k=0}^{N-1} F_k e^{2\pi i n k / N}$$

```{note}
The convention we choose here is the same as used by the NumPy library in [numpy.fft](https://numpy.org/doc/stable/reference/routines.fft.html)
```

```{note}
The discrete Fourier transform is closely related to the Fourier series.
```

Sometimes we write the DFT operation as:

$$F_k = \mathcal{F}(f_n)$$

and the inverse as:

$$f_n = \mathcal{F}^{-1}(F_k)$$

### Significance of real and imaginary parts

Recall that we are integrating with a factor $e^{-2\pi i k x}$.  Euler's formula
tells us that

$$e^{ix} = \cos(x) + i \sin(x)$$

This means that the real part of the transform represents the cosine terms (symmetric functions) and the imaginary part represents the sine terms (antisymmetric).

If $f_n$ is real, then:

$$\mathrm{Re}(F_k) = \sum_{n=0}^{N-1} f_n \cos \left ( \frac{2 \pi n k}{N} \right )$$
$$\mathrm{Im}(F_k) = \sum_{n=0}^{N-1} f_n \sin \left ( \frac{2 \pi n k}{N} \right )$$


```{note}
$k = 0$ is special:

$$\mathrm{Re}(F_0) = \sum_{n=0}^{N-1} f_n \cos \left ( \frac{2 \pi n (0)}{N} \right ) = \sum_{n=0}^{N-1} f_n$$
$$\mathrm{Im}(F_0) = \sum_{n=0}^{N-1} f_n \sin \left ( \frac{2 \pi n (0)}{N} \right ) = 0$$

This is sometimes called the _DC offset_.

```

### Normalization

We haven't yet motivated the $1/N$ in front of the inverse DFT.  This can be shown to arise
from the discrete form of the [Plancherel theorem](https://en.wikipedia.org/wiki/Plancherel_theorem)&mdash;namely that both the real-space
and frequency-space form of the function contain the same power.

But we can motivate it by a simple example.

Consider $f(x) = 1$:

$$ f(x) = 1 \Longleftrightarrow F(k) = \delta(k)$$

If we look at the discrete transform then:

$$F_k = \sum_{n=0}^{N-1} 1 \cdot e^{-2\pi i k n /N}
= \sum_{n=0}^{N-1} \left [ \cos \left (\frac{2\pi kn}{N}\right) + i \sin \left (\frac{2\pi kn}{N}\right ) \right ]$$

* For $k = 0$: $\cos(0) = 1$ and $\sin(0) = 0$, so $F_0 = N$
* For $k > 0$: we are essentially doing a sum over the cosine and sine using equally spaced points, and the sum is always over a multiple of a full wavelength $\rightarrow$ $F_k = 0$

With these frequency-space $F_k$'s, we can compute the inverse transform:

$$f_n = \frac{1}{N} \sum_{k=0}^{N-1} F_k e^{2\pi i nk/N}
 = \frac{1}{N} N e^{2\pi i n 0/N} = 1$$

So we see that the $1/N$ is needed to get back the original value of $f(x)$.