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


