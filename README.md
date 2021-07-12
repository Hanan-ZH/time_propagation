#  About:
This is a simple code that illustrates how time propagation works in a simplified non-linear time-dependent Hamiltonian $`H(t)`$ of the following form:

```math
{H}(t) = \sigma . [  \mu_{\beta}  [  \vec{B_0} +\vec{B_1}(t) ] + \lambda_{||} \langle \vec{\sigma} \rangle +  \lambda_{\perp} \langle \vec{\sigma} \rangle \times \mu_{\beta}  [ \vec{B_0} +\vec{B_1}(t) ]   ] 
```

Where $`\vec{B}_0`$ and $`\vec{B}_1(t)`$ are the static and oscillatory magnetic fields respectively, $`\lambda_{||}`$ and $`\lambda_{||}`$ are constant coefficients, and $`\sigma`$ is the spin pauli vector.

Here, an initial wavefunction (wave-vector) is propagated using the fourth-order embedded Runge-Kutta scheme (**imRK4**). 

# Method references:
[1] Butcher, J.C.The numerical analysis of ordinary differential equations: Runge-Kutta andgeneral linear methods(Wiley-Interscience, 1987).

[[2] ](https://core.ac.uk/download/pdf/196658294.pdf)Rang, J.Apdative Timestep Control for Fully Implicit Runge-Kutta Methods of Higher Order(UnivBibl., 2014).



# Compile:
`cmake ../ -DPLATFORM=system [-DCOMPILER=compiler -DDEBUG=ON]`

`make`

# Run:
To run on Linux:

`mpirun -np 2 ~/time_propagation/bin/tp.exe`

# Test:
An example calculation is in model1 file.

