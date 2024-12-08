C++ Implementation
==================

Here's an implementation of Gaussian elimination in C++, with partial pivoting.
This uses a two-dimensional multidimensional array class, ``Array`` defined
in the header ``array.H``.  This allows us to do things like:

```c++
Array A{{{1, 1, 1},
         {-1, 2, 0},
         {2, 0, 1}}};
```

To initialize a 2-d array.  We can then get the number of rows and
columns via ``A.nrows()`` and ``A.ncols()``.  We can also index the
array using ``A(i,j)`` to get the ``i``-th row and ``j``-th column.

Some details of this class can be found here: https://zingale.github.io/phy504/contiguous_array.html
(Although the implementation here is slightly different).

Here's the class header ``array.H``:

```{literalinclude} ../../../examples/linear-algebra/array.H
:language: c++
```

A function to do a matrix-vector multiply is contained in ``matmul.H``:

```{literalinclude} ../../../examples/linear-algebra/matmul.H
:language: c++
```

The Gaussian elimination function is in the header ``gauss.H``:

```{literalinclude} ../../../examples/linear-algebra/gauss.H
:language: c++
```

Finally, a test driver ``gauss_test.cpp`` is here:

```{literalinclude} ../../../examples/linear-algebra/gauss_test.cpp
:language: c++
```

To run this, we would do:

```bash
g++ -I . --std=c++17 -o gauss_test gauss_test.cpp
./gauss_test
```

Try it
------

Let's test this out on the system we used when we introduced pivoting:

\begin{alignat*}{4}
 \epsilon x &+ y    &+ z  &= 6 \\
-x &+ 2y   &     &= 3 \\
2x &+      &+ z  &= 5
\end{alignat*}
