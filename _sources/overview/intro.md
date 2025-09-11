# Course Overview

This is a special topics class on computational astrophysics.  The
main prerequisite is that you have taken PHY 277: *Computation for
Physics and Astronomy*.

## Goals

### Understanding algorithms

The main goal of this class is to understand how the core numerical
methods that we use to solve problems in astrophysics work.  Today,
you can find libraries that implement much of these algorithms for us,
with nice, easy-to-use interfaces.  For example, the [SciPy
library](https://scipy.org/) for python.

However, before we use a "black box" implementation of a method, we
should understand the basics of the algorithm, how it works, and
importantly, what are the underlying assumptions.

For this reason, we will implement all of the methods here from
scratch, first, before moving onto any libraries.

### Topics

There are a large number of interesting numerical methods used in astrophysics,
and we will not have time to cover everything.  Some references will be given
throughout on how you can learn about other methods.

Topics in the table of contents marked with an "✱" will likely not be covered, but
are self-contained and can be read through at your own pace.

### Programming experience


The second goal of this class is to get everyone more comfortable
with programming, so you should ask questions about programming as
they arise during the semester.

Usually, PHY 277 teaches Fortran and C++.  In this course, I will also
use python, since it makes interactive coding during the lectures a
lot easier.

```{note}
You are free to use whichever programming language you are
most comfortable with for assignments.  I will do interactive
demonstrations primarily in python and I will provide C++
implementations of all of the main algorithms we discuss in class.
```

```{warning}
I will use C++17.  For the GNU compilers, you can specify this
via `--std=c++17`.
```


Some resources:

* My python class notes are available at: https://sbu-python-class.github.io/python-science/Introduction.html

* My C++ class notes are available at: https://zingale.github.io/phy504/

```{note}
We will start the semester with a "crash course" on python, and we
will learn more about the language as the semester goes on and we
implement the core numerical methods and solve interesting problems.

You are *strongly encouraged* to install python on your own computer.
If you don't already have python, I suggest installing the 
(freely-available) Anaconda python distribution:

https://www.anaconda.com/download

This will install everything you will need for this course.
```
