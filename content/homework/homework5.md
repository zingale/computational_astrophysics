# Homework 5

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```

1. We are going to create a simple two-dimensional table of pressure
   (in terms of density and temperature) and interpolate from it to get the pressure
   at any thermodynamic point.

   We'll use a simple expression for pressure, but this will mimic what is done in
   real stellar evolution codes:

   $$P = \frac{k_B \rho T}{\mu_I m_u} + \frac{1}{3} a T^4 + K \left (\frac{\rho}{1~\mathrm{g~cm^{-3}}} \right )^{5/3} \mu_e^{-5/3}$$

   this is ideal gas for ions (first term), radiation (middle term), and non-relativistic electrons (last term).  The composition is specified by $\mu_I$
   (the mean molecular weight per ion) and $\mu_e$ (the mean molecular weight per electron).


   The values of the constants are:

   * $k_b = 1.38\times 10^{-16}~\mathrm{K}$
   * $a = 7.56\times 10^{-15}~\mathrm{erg~cm^{-3}}$
   * $\mu_I = 1.26$
   * $\mu_e = 1.15$
   * $K = 10^{13}~\mathrm{erg~cm^{-3}}$

   Our goal is to create a two-dimensional table of pressure at discrete
   values of density and temperature.  We'll construct this in the following way:

   * Use $N_\rho$ points for density between $1~\mathrm{g~cm^{-3}}$ and $10^4~\mathrm{g~cm^{-3}}$, spaced equally logarithmically.

     ```{note}
     In most programming languages, the `log()` function is the natural log.  We want
     The base-10 log, so you'll want to use `log10()`
     ```
     
   * Use $N_T$ points for temperature between $10^6~\mathrm{K}$ and $10^8~\mathrm{K}$,
     spaced equally logarithmically.

   * Compute the pressure at each of those points.  This will create a grid that
     looks like below:

   ![simple 2-d grid showing pressure](bilinear_eos.png)

   Now we need to intepolate in two-dimensions.  Imagine that the orange "$\times$"
   is wher we want to know the pressure.  We will use the 4 tabulated points surrounding it: $(\rho_i, T_j)$, $(\rho_{i+1}, T_j)$, $(\rho_i, T_{j+1})$, $(\rho_{i+1}, T_{j+1})$

   * Write our interpolating function as:

     $$P(\rho, T) = a (\rho - \rho_i)(T - T_j) + b (\rho - \rho_i) + c (T - T_j) + d$$

     This is called _bilinear interpolation_.
     Compute the 4 unknowns, $a$, $b$, $c$, and $d$ using the 4 data points surrounding
     our $\times$.

   * Write a function to return the interpolated value of pressure given $(\rho, T)$.
