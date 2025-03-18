# C++ Implemntation of MNIST

Here we describe a C++ implementation of the neural network (with
hidden layer) for the MNIST digits.

  * {download}`matrix.H <../../examples/machine_learning/mnist/matrix.H>`

    This builds on the previous `Matrix` class by implementing
    additional operators, including an element-by-element
    multiplication (using the `%` operator) and multiplication by
    a scalar.

    ```{note}
    This is the same operator that the [armadillo library](https://arma.sourceforge.net/)
    uses for element-by-element multiplication.
    ```

    The `Matrix` class also implements addition and subtraction of matrices
    and element-by-element addition / subtraction of a scalar, and
    the compound operators `+=` and `-=`.

    Finally, an `apply()` method can be used to apply a function to
    a matrix on an element-by-element basis.

    For example:

    ```
    Matrix B{{1, 2}, {3, 4}};

    B.apply([](double x) -> double {return 2*x+1;});
    ```

* {download}`mnist.H <../../examples/machine_learning/mnist/mnist.H>`

  This provides functions for reacting the test and training data sets
  (CSV files) as well as a `MNISTDigit` class that manages the data
  and representation for a single digit.

  For example, we the training set and look at the first digit as:

  ```
  auto training_set = mnist::read_training_set();
  training_set[0].display();
  ```

  This will display:

  ```






                  ▒▒▓ ▓██▒
            ░▓▓██████▓██▓░
         ░██████████░░░░
          ██████▓▓██
          ░▓▒███   ▓
             ▓█░
             ▒█▓
              ▓█░
               ██▓▒
               ░███▒
                 ▓██▒
                  ░██▓
                   ███░
                 ▒▓███
               ▒█████▓
             ▒█████▓░
           ░█████▓░
         ▓█████▓░
      ░▓██████▒
      ▒████▒▒


  ```


* {download}`neural_network.H <../../examples/machine_learning/mnist/neural_network.H>`

  This is the main network, which provides a `NeuralNetwork` class that manages the
  weights (including a single hidden layer).

* {download}`test_nn.cpp <../../examples/machine_learning/mnist/test_nn.cpp>`

  This is a test driver that first trains the network using the training set
  and then evaluates how we do by reporting the fraction of the test set we get
  correct.

  It will also display the first 10 digits we get wrong, along with our prediction
  and the right answer.

  ````{note}
  This requires C++23, and should be compiled with optimization like:

  ```
  g++ -O3 -DNDEBUG -std=c++23 -o test_nn test_nn.cpp
  ```
  ````

  ```{note}
  To run, the unzipped CSV files need to be in the same directory as the
  executable.
  ```
