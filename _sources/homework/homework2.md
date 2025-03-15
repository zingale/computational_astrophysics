# Homework 2

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

1. Stellar colors are expressed as the difference in [astronomical
   magnitudes](https://en.wikipedia.org/wiki/Magnitude_(astronomy))
   with the flux measured through different filters.

   The standard [photometric system](https://en.wikipedia.org/wiki/Photometric_system)
   defines U, B, V, and more filters, where U is ultraviolet, B is blue, and V
   from which we can define a [color index](https://en.wikipedia.org/wiki/Color_index)
   like $B-V$ (where here, $B$ is now the magnitude of the star through the B filter).
   It is $B-V$ that traditionally appears as the $x$-axis in the
   [Hertzsprung-Russel diagram](
   https://en.wikipedia.org/wiki/Hertzsprung%E2%80%93Russell_diagram).

   $B-V$ is a difference in magnitudes, which can be expressed as:

   $$B - V = -2.5 \log_{10}(f_B/f_V)$$

   where $f_B$ is the flux received from a star through the B filter and $f_V$
   is the flux from a star through the V filter.

   We can treat a star as a perfect blackbody (narrator: they are not), and then
   take

   $$f_X \approx \int_0^\infty B_\lambda(T) s_X(\lambda) d\lambda$$

   where $B_\lambda(T)$ is the Planck function for temperature $T$, and $s_X$ is the filter X's response
   (I know, too many things are getting called B...).
   The "$\approx$" here is because the distance to the star, its radius, and the fact that we need to
   integrate the Planck function over the outward angles all introduce a proportionality constant.  However, since
   we are taking the ratio of fluxes, this proportionality constant will cancel.
   
   We will take the filter to be a box-filter, with:

   $$s_X(\lambda) = \left \{ \begin{array}{cc} 1 & |\lambda - \lambda_{0,X}| \le  \frac{d\lambda_X}{2} \\
                                               0 & \mbox{otherwise} \end{array} \right . $$

   (i.e., the filter only lets light through in a narrow band of width $d\lambda_X$
   centered on a wavelength $\lambda_{0,X}$.)  Reasonable values for the filter properties
   are:

   | filter   | $\lambda_{0,X}$  | $d\lambda_{X}$  |
   |----------|------------------|-----------------|
   |  B       |    445 nm        |     90 nm       |
   |  V       |    550 nm        |     90 nm       |
   

   There is one further complication.  Magnitudes have a zero-point (calibration).  It used to be
   that Vega was taken to have a magnitude of 0 in all filters (but Vega is variable, so a different
   procedure is used now).  We will assume something Vega-like, and define $B - V = 0$ for $T = 10000~\mathrm{K}$.
   This means that you will need to shift your results by a constant offset.
   
   ```{admonition} Your task
   Using either the trapezoid rule or Simpson's rule, write a code that computes
   $B-V$ for 20 temperatures evenly spaced in $\log_{10}(T)$ between 4000 K and 25000 K, and
   plot $B-V$ vs. $T$.

   Make sure you demonstrate that you used enough bins in your integral (perhaps by doubling the
   number of bins and showing that the change in the answer is small).
   
   The add the following data points of $B-V$ from real stars (data taken from Carroll & Ostlie):

   | spectral type |  $T$    | $B - V$  |
   |:-------------:|--------:|-----------
   |   B0          | 30000   | -0.30    |
   |   A0          |  9520   | -0.02    |
   |   F0          |  7200   | +0.30    |
   |   G0          |  6030   | +0.58    |
   |   K0          |  5250   | +0.81    |
   |   M0          |  3850   | +1.40    |
   
   