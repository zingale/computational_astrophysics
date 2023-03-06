# Homework 4

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```


1. A convolution is defined as:

   $$(f \star g)(t) \equiv \int_{-\infty}^{\infty} f(\tau) g(t - \tau) d\tau$$

   It is easy to compute this with FFTs, via the {\em convolution theorem},

   $$\mathcal{F}\{f \star g\} = \mathcal{F}\{f\} \, \mathcal{F}\{g\}$$

   That is the Fourier transform of the convolution of $f$ and $g$ is simply
   the product of the individual transforms of $f$ and $g$.  This allows us
   to compute the convolution via multiplication in Fourier space and then take
   the inverse transform, $\mathcal{F}^{-1}\{\}$, to recover the convolution
   in real space:

   $$f \star g = \mathcal{F}^{-1}\{ \mathcal{F}\{f\} \, \mathcal{F}\{g\}\}$$

   The file [`signal.txt`](signal.txt) contains data of a
   function polluted with noise.  We want to remove the noise to
   recover the original function.  The three columns in the file are:
   $x$, $f^\mathrm{(orig)}(x)$, and $f^\mathrm{(noisy)}(x)$.

   Consider the following kernel:

   $$q^\mathrm{gauss}(x) = \frac{1}{\sigma \sqrt{2 \pi}}\, e^{-\frac{1}{2} (x/\sigma)^2}$$

   * Make the kernel periodic on the domain defined by the $x$ in
     the `signal.txt` file.  You can do this simply by left-right
     flipping the definitions above and applying them at the far end of
     the domain.

     Make sure that your kernel function by ensuring that it sums to $1$
     on the domain.  You might need to sum up the values and divide by the sum.

   * Plot the noisy function, $f^\mathrm{(noisy)}(x)$, and the kernel together.

   * Take the FFT (or DFT) of $f^\mathrm{(noisy)}(x)$ and the FFT of the kernel and
     plot them.

   * Compute the convolution of $f^\mathrm{(noisy)}(x)$ and $q(x)$
     in Fourier space and transform back to real space, and plot the
     {\em de-noised} function together with the original signal from
     the `signal.txt` (i.e., $f^\mathrm{(orig)}(x)$).

   * Experiment with the tunable parameters in the kernels ($L$ or $\sigma$)
     to see how clean you can get the noisy data and comment on what you see.

  This process is used a lot in image processing both to remove noise
  and to compensate for the behavior of cameras to sharpen images.


