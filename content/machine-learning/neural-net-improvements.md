# Improvements

There are many ways we could improve the performance of our network, but these will add
a lot of complexity to the simple class that we wrote.  Fortunately there are a lot of
machine learning libraries that provide these features, and work efficiently, so for
real applications we would want to use one of those libraries (we'll explore these next).

## Batching

Right now, we did our training as:

* Loop over the $T$ pairs $({\bf x}^k, {\bf y}^k)$ for $k = 1, \ldots, T$

  * Propagate $({\bf x}^k, {\bf y_k})$ through the network
  * Compute the corrections $\partial f/\partial {\bf A}$, $\partial f/\partial {\bf B}$
  * Update the matrices:
  
    $${\bf A} \leftarrow {\bf A} + \eta \frac{\partial f}{\partial {\bf A}}$$

    $${\bf B} \leftarrow {\bf B} + \eta \frac{\partial f}{\partial {\bf B}}$$

In this manner, each training pair sees slightly different
matrices ${\bf A}$ and ${\bf B}$, as each previous pair
updates it immediately.

We could instead divide our training set into $N$ batches,
each with $\tau = T/N$ training pairs and do our update as:

* Loop over $N$ batches

  * Loop over the $\tau$ pairs $({\bf x}^k, {\bf y}^k)$ for $k = 1, \ldots, \tau$ in the current batch 

    * Propagate $({\bf x}^k, {\bf y_k})$ through the network
    * Compute the corrections $\partial f/\partial {\bf A}^k$, $\partial f/\partial {\bf B}^k$ from the current pair
    
    * Accumulate the corrections:
  
      $$\frac{\partial f}{\partial {\bf A}} = \frac{\partial f}{\partial {\bf A}} + \frac{\partial f}{\partial {\bf A}^k}$$
      
  * Apply a single update to the matrices for this batch:

    $${\bf A} \leftarrow {\bf A} + \eta \frac{\partial f}{\partial {\bf A}}$$

    $${\bf B} \leftarrow {\bf B} + \eta \frac{\partial f}{\partial {\bf B}}$$

The advantage of this is that the $\tau$ trainings in a batch
can all be done in parallel now, spread across many CPU cores
or GPU cores.  This greatly accelerates the training time.
  

## Different activation or cost functions

We used a simple cost function: the sum of the square of the errors.  This is analogous to the $L_2$ norm we discussed previously.  But there are a lot of other cost functions
we could explore.  Changing the cost function will require
us to recompute our derivatives.

Likewise, there are a wide number of activation functions,
some of which are not differentiable.  The choice of activation
function can depend on what type of data you are using.  You
might also want to use a different activation function
on each layer.  Again, this would require us to redo
our derivatives.


## Use a different minimization technique


## Different types of layers / connections


## More hidden layers


## Auto-differentiation libraries

