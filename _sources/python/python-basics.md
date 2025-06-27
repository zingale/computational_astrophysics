Python Basics
=============

Python differs from C++ or Fortran in several ways:

* Python is interpreted: There is no compilation step that turns the source code
  into machine instructions.  Instead, you run python programs through an interpreter,
  usually called `python`.

  In our case, Jupyter will be our interface to the python interpreter.

* Python is dynamically typed: You don't need to declare what type of data a variable
  will hold ahead of time.

* Python manages memory automatically: you don't have to explicitly free memory
  when it is no longer needed.  Python uses *reference counting* to determine if
  memory is still needed.

  Note: with modern C++, if you use the containers defined in the
  standard library (like `vector`), then memory management is largely
  taken care of for you as well.  But it is still easy to get yourself
  in trouble in C++ if you are not careful.

* Python code is in general slower than compiled C++ or Fortran code.
  This is because it is interpreted.  However, many libraries that we
  use through python are built in C and if we use them properly, we
  can recover a lot of the speed difference compared to a compiled
  language.

