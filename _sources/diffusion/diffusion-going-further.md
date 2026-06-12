# Going Further

## Crank-Nicolson

The [Crank-Nicolson](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method) method uses a time-centered
discretization of the diffusion term:

$$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta x} = \frac{k}{2}
         \left (\frac{\phi_{i+1}^n - 2\phi_i^n + \phi_{i-1}^n}{\Delta x^2} +
                \frac{\phi_{i+1}^{n+1} - 2\phi_i^{n+1} + \phi_{i-1}^{n+1}}{\Delta x^2}\right )$$
                
This is still an implicit update, but the coefficients and righthand side
in the linear system will be different

## State-dependent transport coefficients

If we have a non-constant conductivity that depends on $\phi$ itself,
then our discretization is more complicated:

$$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta t} =
        \frac{ \{ k \nabla \phi \}_{i+1/2} -
               \{ k \nabla \phi \}_{i-1/2}}{\Delta x}$$
               
Now we need $k$ at the interfaces.  There are several ways to do this, which
should be motivated by the physics of what the coefficient $k$ represents:

$$k_{i+1/2} = \frac{1}{2} (k_i + k_{i+1})$$

or

$$\frac{1}{k_{i+1/2}} = \frac{1}{2} \left ( \frac{1}{k_i} + \frac{1}{k_{i+1}} \right )$$

Then a C-N update appears as:

$$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta t} =
        \frac{1}{2} \left \{
               \nabla \cdot [ k(\phi^n) \nabla \phi^n ]_i +
               \nabla \cdot [ k(\phi^{n+1}) \nabla \phi^{n+1} ]_i
               \right \}$$

We can update this using a _predictor-corrector_ scheme:

$$\frac{\phi_i^\star - \phi_i^n}{\Delta t} =
        \frac{1}{2} \left \{
               \nabla \cdot [ k(\phi^n) \nabla \phi^n ]_i +
               \nabla \cdot [ k(\phi^n) \nabla \phi^\star ]_i
               \right \}$$
               
$$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta t} =
        \frac{1}{2} \left \{
               \nabla \cdot [ k(\phi^n) \nabla \phi^n ]_i +
               \nabla \cdot [ k(\phi^\star) \nabla \phi^{n+1} ]_i
               \right \}$$
               
This will be second-order accurate in time.


