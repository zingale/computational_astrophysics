Multivariate Root Finding
=========================

Imagine a vector function,

$${\bf f}({\bf x}) = \left ( \begin{array}{c} f_1({\bf x}) \\ f_2({\bf x}) \\ \vdots \\ f_N({\bf x}) \end{array} \right )$$

where

$${\bf x} = \left ( \begin{array}{c} x_1 \\ x_2 \\ \vdots \\ x_N \end{array} \right )$$

We want to find the ${\bf x}$ that zeros ${\bf f}({\bf x})$, i.e.,

$${\bf f}({\bf x}) = 0$$

We'll use a generalization of Newton's method to systems of equations.

Solution procedure
------------------

```{note}
This is the simplest case of a multivariable root finding algorithm.  More
sophisticated methods exist, but this method will give you a feel for what
is involve
```

Start with an initial guess, ${\bf x}^{(0)}$

Taylor expand our function seeking a correction $\delta {\bf x}$.  Each component has the form:
  
$$f_i ({\bf x}^{(0)} + \delta {\bf x}) \approx 0 = f_i({\bf x}^{(0)}) + \sum_{j=1}^N \frac{\partial f_i}{\partial x_j} \delta x_j + \ldots$$
  
and we can write the vector form as:
  
$${\bf f}({\bf x}^{(0)} + \delta {\bf x}) = {\bf f}({\bf x}^{(0)}) + {\bf J} \delta {\bf x} + \ldots \approx 0$$
  
where ${\bf J}$ is the [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant):
  
$${\bf J} = \left ( \begin{array}{cccc} \frac{\partial f_1}{\partial x_1} &
                                        \frac{\partial f_1}{\partial x_2} &
                                        \cdots &
                                        \frac{\partial f_1}{\partial x_N} \\
%
                                        \frac{\partial f_2}{\partial x_1} &
                                        \frac{\partial f_2}{\partial x_2} &
                                        \cdots &
                                        \frac{\partial f_2}{\partial x_N} \\
%
                                        \vdots & \vdots & \ddots & \vdots \\
%
                                        \frac{\partial f_N}{\partial x_1} &
                                        \frac{\partial f_N}{\partial x_2} &
                                        \cdots &
                                        \frac{\partial f_N}{\partial x_N} \end{array} \right )$$


This can be expresses as the linear system:

$${\bf J}\delta {\bf x} = - {\bf f}({\bf x}^{(0)})$$

Our solution procedure is iterative:

* Iterate over $k$, seeking a new guess at the solution ${\bf x}^{k+1}$:

  * Solve the linear system:

    $${\bf J}\delta {\bf x} = - {\bf f}({\bf x}^{(k)})$$
  
  * Correct our initial guess:

    $${\bf x}^{(k+1)} = {\bf x}^{(k)} + \delta {\bf x}$$
  
  * Stop iteration if

    $$\| \delta {\bf x} \| < \epsilon \| {\bf x}^{(k)} \|$$
  
    where $\| \cdot \|$ is a suitable vector norm.
