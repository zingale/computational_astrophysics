# Hidden Layers


 We can get better performance from a neural network by adding a hidden layer:

![hidden layers](nn_fig_hidden.png)

The size of the hidden layer is independent of the size of the input and output layers.  In this case, we have a hidden layer that is larger
than either the input or output layers.

Now we have an additional matrix ${\bf B}$ to train.  This can all be done together using the same algorithm described above.  Where we now minimze:

$$f(A_{lm}, B_{ij}) = \sum_{l=1}^{N_\mathrm{out}} (z_l - y_l)^2$$

$$\tilde{z}_i = g \left ( \sum_{j=1}^{N_\mathrm{in}} B_{ij} x_j \right )$$

$$z_l = g \left ( \sum_{m=1}^{N_\mathrm{hidden}} A_{lm} \tilde{z}_m \right )$$
