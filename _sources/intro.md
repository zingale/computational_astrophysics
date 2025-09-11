AST 390: Computational Astrophysics
===================================

Michael Zingale

(Spring 2025)

This is a collection of notebooks on computational (astro)physics.
Starting at the beginning, these notebooks introduce numerical methods
for derivatives, integration, rooting finding, ODEs, and linear algebra
and then move onto applications in astrophysics.


Course Outline
--------------

We'll start with an [An Overview of Python](https://zingale.github.io/computational_astrophysics/python/jupyter.html)
and [Version Control with Git](https://zingale.github.io/computational_astrophysics/git/version-control.html) and then
move on to core numerical methods.

You should follow the outline in the navigation panel to the left.

```{note}
This course assumes that you are already familiar with a programming language.
```


Astrophysical Applications
--------------------------

Throughout the course, we'll see some applications to interesting problems in astrophysics.  Here's a listing
to some of them:

::::{grid}
:gutter: 3

:::{grid-item-card} Reaction Rate $T$ Sensitivity
{bdg-info-line}`diff`
An example of using numerical differentiation to
[estimate the 3-$\alpha$ reaction $T$ sensitivity](https://zingale.github.io/computational_astrophysics/basics/diff-int/application-rate-temperature-sensitivity.html).
:::

:::{grid-item-card} Blackbody Radiation
{bdg-info-line}`int`
An example of integrating to infinity by
[integrating the Planck function over wavelength](https://zingale.github.io/computational_astrophysics/basics/diff-int/application-blackbody.html).
:::

:::{grid-item-card} Wien's Law
{bdg-info-line}`roots`
Demonstrating root finding by
[numerically deriving Wien's law](https://zingale.github.io/computational_astrophysics/basics/roots/application-wiens.html).
:::


::::

::::{grid}
:gutter: 3

:::{grid-item-card} Degeneracy
{bdg-info-line}`int` {bdg-info-line}`roots`
Combining integration over the Fermi-Dirac
distribution and root-finding to [find
the electron degeneracy parameters](https://zingale.github.io/computational_astrophysics/basics/roots/application-degeneracy.html).
:::

:::{grid-item-card} Few-Body Problem
{bdg-info-line}`ODEs`
Using adaptive stepping in ODE integration
to [solve the few-body problem](https://zingale.github.io/computational_astrophysics/ODEs/application-few-body.html).
:::

:::{grid-item-card} Polytropes
{bdg-info-line}`roots` {bdg-info-line}`ODEs`
Shooting methods for two-point boundary
value problelms
applied to the [Lane-Emden equation for polytropes](https://zingale.github.io/computational_astrophysics/ODEs/application-lane-emden.html).
:::

::::


::::{grid}
:gutter: 3


:::{grid-item-card} Stability of Planetary Systems
{bdg-info-line}`ODEs`
Investigating the [longterm stability of planetary systems](https://zingale.github.io/computational_astrophysics/ODEs/application-planetary-stability.html)
using symplectic integrators.
:::

:::{grid-item-card} One-Zone X-ray Burst Model
{bdg-info-line}`roots` {bdg-info-line}`ODEs`
ODE integration + root finding
to explore limit-cycles in a 
[one-zone model of an X-ray burst](https://zingale.github.io/computational_astrophysics/ODEs/application-xrb-one-zone.html)
:::

:::{grid-item-card} Integrating the CNO Cycle
{bdg-info-line}`ODEs`
Using stiff-ODE solvers to
[integrate an CNO reaction network](https://zingale.github.io/computational_astrophysics/ODEs/application-CNO-network.html)
:::

::::


::::{grid}
:gutter: 3


:::{grid-item-card} Lorenz System Stationary States
{bdg-info-line}`linalg`
A demonstration of using Newton's method to find
the [stationary states of the Lorenz system](https://zingale.github.io/computational_astrophysics/basics/linear-algebra/application-lorenz.html).
:::

:::{grid-item-card} Interpolating Reaction Rates
{bdg-info-line}`interp`
Using bilinear interpolation to [evaluate tabulated electron-capture rates](https://zingale.github.io/computational_astrophysics/interpolation/application-interpolating-reaction-rates.html).
:::

:::{grid-item-card} Radial Velocity in Frequency Space
{bdg-info-line}`FFTs`
Using FFTs to [find planets from stellar radial velocity curves](https://zingale.github.io/computational_astrophysics/ffts/application-planets-fft.html).
:::


::::

::::{grid}
:gutter: 3

:::{grid-item-card} X-ray timing
{bdg-info-line}`FFTs`
Using FFTs on time-series data to
[study low mass X-ray binaries](https://zingale.github.io/computational_astrophysics/ffts/application-lightcurves.html).
:::

:::{grid-item-card} Hubble's Constant
{bdg-info-line}`fitting`
Using linear regression to [estimate $H_0$ from
Type Ia supernova](https://zingale.github.io/computational_astrophysics/fitting/application-snia-h0.html).
:::

:::{grid-item-card} Galaxy Classification
{bdg-info-line}`ML/NN`
Using machine learning / neural networks to [classify galaxy images](https://zingale.github.io/computational_astrophysics/machine-learning/application-galaxies-intro.html).
:::

::::
