# Deriving the Learning Correction

For gradient descent, we need to derive the update to the matrix
${\bf A}$ based on training on a set of our data, $({\bf x}^k, {\bf y}^k)$.

Let's start with our cost function:

$$f(A_{ij}) = \sum_{i=1}^{N_\mathrm{out}} \left [ g\left (\sum_{j=1}^{N_\mathrm{in}} A_{ij} x^k_j \right ) - y^k_i \right ]^2$$

and we'll compute the derivative with respect to a single matrix
element, $A_{pq}$:

$$\frac{\partial f}{\partial A_{pq}} =
  2 \sum_{i=1}^{N_\mathrm{out}} (z_i - y^k_i) \frac{\partial z_i}{\partial A_{pq}}$$
  
Now we apply the chain rule:

$$\frac{\partial z_i}{\partial A_{pq}} = \frac{\partial g(b_i)}{\partial A_{pq}} = \left . \frac{\partial g}{\partial \xi} \right |_{\xi=b_i} \frac{\partial b_i}{\partial A_{pq}}$$

with

$$\frac{\partial b_i}{\partial A_{pq}} = \sum_{j=1}^{N_\mathrm{in}} \frac{\partial A_{ij}}{\partial A_{pq}} x^k_j = \sum_{j=1}^{N_\mathrm{in}} \delta_{ip} \delta_{jq} x^k_j = \delta_{ip} x^k_q$$

and for $g(\xi)$, we will assume the sigmoid function,

$$g(\xi) = \frac{1}{1 + e^{-\xi}}$$

so

$$\frac{\partial g}{\partial \xi} = - (1 + e^{-\xi})^2 (- e^{-\xi})
  = g(\xi) \frac{e^{-\xi}}{1+ e^{-\xi}} = g(\xi) (1 - g(\xi))$$

which gives us:

$$\frac{\partial z_i}{\partial A_{pq}} = g(b_i)(1 - g(b_i)) \delta_{ip} x^k_q = z_i (1- z_i) \delta_{ip} x^k_q$$

Then finally, we have:

\begin{align}
\frac{\partial f}{\partial A_{pq}} &= 2 \sum_{i=1}^{N_\mathrm{out}}
   (z_i - y^k_i) z_i (1 - z_i) \delta_{ip} x^k_q \\
   &= 2 (z_p - y^k_p) z_p (1- z_p) x^k_q
\end{align}
   
where we used the fact that the $\delta_{ip}$ means that only a single term contributes to the sum.

Note that:

* $e_p^k \equiv (z_p - y_p^k)$ is the error on the output layer,
  and the correction is proportional to the error (as we would
  expect).

* The $k$ superscripts here remind us that this is the result of
  only a single pair of data from the training set.
  
Now ${\bf z}$ and ${\bf y}^k$ are all vectors of size $N_\mathrm{out} \times 1$ and ${\bf x}^k$ is a vector of size $N_\mathrm{in} \times 1$, so we can write this expression for the matrix as a whole as:

$$\frac{\partial f}{\partial {\bf A}} = 2 ({\bf z} - {\bf y}^k) \circ {\bf z} \circ (1 - {\bf z}) \cdot ({\bf x}^k)^\intercal$$

where the operator $\circ$ represents _element-by-element_ multiplication (the [Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)).

