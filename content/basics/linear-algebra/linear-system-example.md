Solving A System of Linear Equations
====================================

Consider the system:

\begin{alignat*}{4}
 x &+ y    &+ z  &= 6 \\
-x &+ 2y   &     &= 3 \\
2x &+      &+ z  &= 5
\end{alignat*}

We will solve this in steps.

1. Use the first equation to eliminate $x$ in the second two equations, resulting in:

   \begin{alignat*}{4}
    x &+ y    &+ z  &= 6 \\
      &+ 3y   &+ z  &= 9 \\
      &- 2y   &- z  &= -7
   \end{alignat*}

2. Use the second equation to eliminate $y$ in the final equation, resulting in:

   \begin{alignat*}{4}
    x &+ y    &+ z  &= 6 \\
      &+ 3y   &+ z  &= 9 \\
      &       &- \tfrac{1}{3}z  &= -1
   \end{alignat*}

   At this point, we completed _forward elimination_ &mdash; the last
   equation has a single unknown, $z$, and once we solve for it, we can
   substitute the value of $z$ in the equation above, leaving it with a
   single unknown, $y$, and so forth.

3. Now we do _back substitution_.

   From the last equation, we see $z = 3$. Putting this into the
   equation above, we get $y = 2$.  Putting both of these into the
   first equation, we finally find $x = 1$.

The process we just did is the basis of _Gaussian elimination_.

```{note}
The only source of error in the procedure we did would be round-off error.  There is no trunction error.
```

```{note}
You may have learned [Cramer's rule](https://en.wikipedia.org/wiki/Cramer's_rule) for solving a linear system
in school.  But for practical purposes, this method is usually much more computationally expensive than
Gaussian elimination.
```


Matrix form
-----------

We can redo this same procedure in matrix form, starting with writing the system as:

$$
{\bf A}{\bf x} = {\bf b}
$$

with

$$
{\bf A} = \left ( \begin{array}{ccc}
                     1  &  1  &  1 \\
                    -1  &  2  &  0 \\
                     2  &  0  &  1 \end{array} \right )
$$

and

$$
{\bf b} = \left ( \begin{array}{c} 6 \\ 3 \\ 5 \end{array} \right )
$$

We will write an [augmented matrix](https://en.wikipedia.org/wiki/Augmented_matrix) $({\bf A}|{\bf b})$ as

$$
({\bf A}|{\bf b}) = \left ( \begin{array}{ccc|c}
                     1  &  1  &  1 & 6 \\
                    -1  &  2  &  0 & 3 \\
                     2  &  0  &  1 & 5 \end{array} \right )
$$

This allows us to easily do the same steps to both ${\bf A}$ and ${\bf b}$ together.

Doing the same sequence of steps as above, the augmented matrix looks as:


$$
({\bf A}|{\bf b}) = \left ( \begin{array}{ccc|c}
                     1  &  1  &  1 & 6 \\
                    -1  &  2  &  0 & 3 \\
                     2  &  0  &  1 & 5 \end{array} \right )
$$

$$
({\bf A}|{\bf b}) = \left ( \begin{array}{ccc|c}
                     1  &  1  &  1 & 6 \\
                     0  &  3  &  1 & 9 \\
                     0  & -2  & -1 & -7 \end{array} \right )
$$

$$
({\bf A}|{\bf b}) = \left ( \begin{array}{ccc|c}
                     1  &  1  &  1 & 6 \\
                     0  &  3  &  1 & 9 \\
                     0  &  0  &  \tfrac{1}{3} & -1 \end{array} \right )
$$

At this point, the matrix is in [row echelon](https://en.wikipedia.org/wiki/Row_echelon_form), and it
clearly shows the structure needed for back-substitution.


Pivoting
--------

The above case worked well because all of the matrix elements were
roughly the same magnitude.  Consider now the case where one of the
matrix elements is small.


