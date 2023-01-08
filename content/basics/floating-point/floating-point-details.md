Floating Point Format
=====================

```{note}
   Reading: [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://dl.acm.org/doi/10.1145/103162.103163)
```

Storage overview
----------------

We can think of a floating point number as having the form:

$$
   \mbox{significand} \times 10^\mbox{exponent}
$$

Most computers follow the IEEE 754 standard for floating point, and we commonly work
with 32-bit and 64-bit floating point numbers (single and double precision).  These bits
are split between the signifcand and exponent as well as a single bit for the sign:

```{figure} 1024px-IEEE_754_Double_Floating_Point_Format.svg.png
---
width: 80%
---
(source: wikipedia)
```

Since the number is stored in binary, we can think about expanding a number in powers of 2:

$$
   0.1 \sim (1 +
             1 \cdot 2^{-1} +
             0 \cdot 2^{-2} +
             0 \cdot 2^{-3} +
             1 \cdot 2^{-4} +
             1 \cdot 2^{-5} + \ldots) \times 2^{-4}
$$

In fact, ``0.1`` cannot be exactly represented in floating point:

C++:

```{literalinclude} ../../../examples/floating_point/simple_roundoff.cpp
:language: c++
```

python:

```python
print(f"{0.1:30.20f}")
```

Precision
---------

With 52 bits for the significand, the smallest number compared to ``1`` we can represent is

$$
   2^{-52} \approx 2.22\times 10^{-16}
$$

but the IEEE 754 format always expresses the significant such that the
first bit is ``1`` (by adjusting the exponent) and then doesn't need
to store that ``1``, giving us an extra bit of precision, so the
machine epsilon is

$$
   2^{-53} \approx 1.11\times 10^{-16}
$$


```{note}
This is a relative error, so for a number like ``1000`` we could only add
``1.1e-13`` to it before it became indistinguishable from ``1000``.

$$
   \mbox{relative roundoff error} = \frac{|\mbox{true number} - \mbox{computer representation} |}
      {|\mbox{true number}|} \le \epsilon
$$
```

There are [varying definitions of machine epsilon](https://en.wikipedia.org/wiki/Machine_epsilon#Variant_definitions) which differ by a factor of 2.

Range
-----

Now consider the exponent, we use 11 bits to store it in double
precision.  Two are reserved for special numbers, so out of the 2048
possible exponent values, one is 0, and 2 are reserved, leaving 2045
to split between positive and negative exponents.  These are set as:

$$ 2^{-1022} \mbox{ to } 2^{1023} $$

converting to base 10, this is

$$ \sim 10^{-308} \mbox{ to } \sim 10^{308} $$


```{note}
In C++, you can get these ranges via ``std::numeric_limits``:
```{literalinclude} ../../../examples/floating_point/limits.cpp
:language: c++
```






Roundoff vs. truncation error
==============================

Consider the Taylor expansion of :math:`f(x)` about some point :math:`x_0`:

.. math::

   f(x) = f(x_0 + \Delta x) = f(x_0) + \left . \frac{df}{dx} \right |_{x_0} \Delta x + \mathcal{O}(\Delta x^2)

where :math:`\Delta x = x - x_0`

We can solve for the derivative to find an approximation for the first derivative:

.. math::

   \left . \frac{df}{dx} \right |_{x_0} = \frac{f(x_0 + \Delta x) - f(x_0)}{\Delta x} + \mathcal{O}(\Delta x)

This shows that this approximation for the derivative is first-order accurate in :math:`\Delta x` -- that is the truncation error of the approximation.

We can see the relative size of roundoff and truncation error by using this approximation
to compute a derivative for different values of :math:`\Delta x`:

.. literalinclude:: ../../examples/floating_point/truncation_vs_roundoff.cpp
   :language: c++
   :caption: ``truncation_vs_roundoff.cpp``

It is easier to see the behavior if we make a plot of the output:

.. figure:: error_plot.png
   :align: center
   :width: 80%

Let's discuss the trends:

* Starting with the largest value of :math:`\Delta x`, as 
  we make :math:`\Delta x` smaller, we see that the error decreases.
  This is following the expected behavior of the truncation error
  derived above.

* Once our :math:`\Delta x` becomes really small, roundoff error starts
  to dominate.  In effect, we are seeing that:

  .. math::

     (x_0 + \Delta x) - x_0 \ne 0

  because of roundoff error.

* The minimum error here is around :math:`\sqrt{\epsilon}`, where :math:`\epsilon` is
  machine epsilon.


Testing for equality
====================

Because of roundoff error, we should never exactly compare two floating point numbers,
but instead ask they they agree within some tolerance, e.g., test equality as:

.. math::

   | x - y | < \epsilon

For example:

.. literalinclude:: ../../examples/floating_point/comparing.cpp
   :language: c++
   :caption: ``comparing.cpp``



Minimizing roundoff
===================

Consider subtracting the square of two numbers -- taking the difference of two very close-in-value numbers is a prime place where roundoff can come into play.

Instead of doing:

.. math::

   x^2 - y^2

we can instead do:

.. math::

   (x - y)(x + y)

by factoring this, we are subtracting more reasonably sized numbers, reducing the roundoff.

We can see this directly by doing this with single precision (`float`) and comparing to an answer computed via double precious (`double`)

Here's an example:

.. literalinclude:: ../../examples/floating_point/subtraction.cpp
   :language: c++
   :caption: ``subtraction.cpp``


As another example, consider computing [#f1]_:

.. math::

   \sqrt{x + 1} - \sqrt{x}

We can alternately rewrite this to avoid the subtraction of two close numbers:

.. math::

   \sqrt{x + 1} - \sqrt{x} = (\sqrt{x + 1} - \sqrt{x})
        \left ( \frac{\sqrt{x+1} + \sqrt{x}}{\sqrt{x+1} + \sqrt{x}} \right )
        = \frac{1}{\sqrt{x+1} + \sqrt{x}}

Again we'll compare a single-precision calculation using each of these methods
to a double precision "correct answer".  To ensure that we use the
single-precision version of the ``std::sqrt()`` function, we will use single
precision literal suffix, e.g., ``1.0f`` tells the compiler that this is a
single-precision constant.

.. literalinclude:: ../../examples/floating_point/squareroots.cpp
   :language: c++
   :caption: ``squareroots.cpp``

Notice that we get several more significant digits correct when we compute it with
the second formulation compared to the original form.

Summation algorithms
--------------------

Summing a sequence of numbers is a common place where roundoff error
comes into play, especially if the numbers all vary in magnitude and
you do not attempt to add them in a sorted fashion.  There are a
number of different summation algorithms that keep track of the loss
due to roundoff and compensate the sum, for example `the Kahan summation algorithm <https://en.wikipedia.org/wiki/Kahan_summation_algorithm#See_also>`_.

Special numbers
===============

IEEE 754 defines a few special quantities:

* ``NaN`` (not a number) is the result of ``0.0/0.0`` or ``std::sqrt(-1.0)``

* ``Inf`` (infinity) is the result of ``1.0/0.0``

* ``-0`` is a valid number and the standard says that ``-0`` is equivalent to ``0``


Trapping floating point exceptions
==================================

What happens when we do something bad?  Consider this example:

.. literalinclude:: ../../examples/floating_point/undefined.cpp
   :language: c++
   :caption: ``undefined.cpp``

Here, we pass ``-1`` to ``trouble()`` which then takes the square root
of it -- this results in a NaN.  But if we run the code, it goes
merrily about its way, using that result in the later computations.

Unix uses `signals <https://en.wikipedia.org/wiki/Signal_(IPC)>`_ to
indicate that a problem has happened during the code execution.  If a
program created a *signal handler* then that signal can be trapped and
any desired action can be taken.

.. note::

   This example was only tested on a Linux machine with GCC.  Other OSes
   or compilers might have slightly different headers or functionality.

There are a few parts to trapping a floating point exception (FPE).  First we need
to enable exception trapping via:

.. code:: c++

   feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

That catches 3 different types of floating point exceptions -- invalid, divide-by-zero, and overflows.

Next we need to add a handler to deal with the exception:

.. code:: c++

   signal(SIGFPE, fpe_handler);

Here, ``SIGFPE`` is the standard name for a floating point exception,
and ``fpe_handler`` is the name of a function that will be called when
we detect a ``SIGFPE``.

In our handler, we use the Linux ``backtrace()`` function to access the stack
of our program execution.  This is really a C-function, so we need to use C-style
arrays here.

Here's the new version of our code:

.. literalinclude:: ../../examples/floating_point/undefined_trap.cpp
   :language: c++
   :caption: ``undefined_trap.cpp``

When we compile the code, we want to add the ``-g`` option to store the
symbols in the code -- this allows us to understand where problems arise:

.. prompt:: bash

   g++ -g -o undefined_trap undefined_trap.cpp



Now when we run this, the program aborts and we see:

::

    floating point exception, signal 8
    0: ./undefined_trap() [0x401261]
    1: /lib64/libc.so.6(+0x42750) [0x7f3dc35dc750]
    2: /lib64/libm.so.6(+0x1435c) [0x7f3dc37d335c]
    3: ./undefined_trap() [0x4012ff]
    4: ./undefined_trap() [0x401347]
    5: /lib64/libc.so.6(+0x2d560) [0x7f3dc35c7560]
    6: /lib64/libc.so.6(__libc_start_main+0x7c) [0x7f3dc35c760c]
    7: ./undefined_trap() [0x401145]
    Aborted (core dumped)

This is the call stack for our program.  In the brackets are the address in
the program where the execution was when the FPE occurred.  These are ordered
such that the calling function is below the function where the execution is.
So it usually is best to look at the addresses near the top.

We can turn those into line numbers using ``addr2line``:

.. prompt:: bash

  addr2line  -e undefined_trap 0x4012ff

gives:

::

   /home/zingale/classes/phy504/examples/floating_point/undefined_trap.cpp:23

and that line is precisely where the ``sqrt()`` is!

We can get slightly nicer output (including the function name) by doing:

.. prompt:: bash

   addr2line  -C -f -i -p -e undefined_trap 0x4012ff

which gives:

::

    trouble(double) at /home/zingale/classes/phy504/examples/floating_point/undefined_trap.cpp:23


.. note::

   On the MathLab machines, the stack trace seems to include an offset, like:

   ::

      floating point exception, signal 8
      0: ./undefined_trap(+0xc03) [0x561d8799dc03]
      1: /lib/x86_64-linux-gnu/libc.so.6(+0x3ef10) [0x7f5461e3df10]
      2: /lib/x86_64-linux-gnu/libm.so.6(+0x11397) [0x7f5462201397]
      3: ./undefined_trap(+0xccf) [0x561d8799dccf]
      4: ./undefined_trap(+0xd21) [0x561d8799dd21]
      5: /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xe7) [0x7f5461e20c87]
      6: ./undefined_trap(+0xaaa) [0x561d8799daaa]
      Aborted (core dumped)

   and we need to use that offset instead with ``addr2line``, like:

   .. prompt:: bash

      addr2line -a -f -e ./undefined_trap +0xcd5


.. [#f1] this example is based on Yakowitz & Szidarovszky
