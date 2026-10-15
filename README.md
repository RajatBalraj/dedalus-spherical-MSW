# dedalus-spherical-MSW
Moist shallow water simulations on a sphere using the [Dedalus]([url](https://dedalus-project.readthedocs.io/en/latest/index.html)) spectral solver. 
Can be run in parallel.

Precipitation and Evaporation parameterised as simple Heaviside, Betts-Miller type functions.

Some example simulation plots are added.


The following SW system is simulated. A $\nabla^4$ hyperviscosity is added to each prognostic equation. Easy to add source terms:


$$ \frac{\partial \vec{u}}{\partial t} + \vec{u} \cdot \vec{\nabla} \vec{u} + f \vec{k} \times \vec{u} = -g \vec{\nabla}h \\[1em]

\frac{\partial h}{\partial t} + \vec{\nabla} \cdot \left(\vec{u} h\right) = -\chi (P-E) \\[1em]

\frac{\partial q}{\partial t} + \vec{\nabla} \cdot \left(\vec{u} q\right) = -P + E\\[1em]

h = H+h'\\[1em]

P = \frac{(q - q_s)}{\tau_c} \Theta \left[q-q_s\right] \\[1em]
E = \frac{(q_s - q)}{\tau_e} \Theta \left[q_s - q\right] \\[1em]

f = 2 \Omega \sin \phi

$$
