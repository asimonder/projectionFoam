# projectionFoam 
A non-iterative OpenFOAM solver to simulate laminar and turbulent single-phase flows. The solver uses the projection algorithm for velocity-pressure decoupling.

## 2023 
Updated the solver for OpenFOAM v2006 and above. 

Benchmarks are in the works.

## 2013
The pisoFoam solver to solve incompressible transient flows is iterative. Here, it is replaced with an incremental projection scheme which is non-iterative therefore costs less than pisoFoam. In this classical scheme, we first solve the momentum equation using the old time step pressure. Then, unsoleinoidal part of the velocity is projected out by solving the Poisson equation. The scheme is mass-conservative and has an additional second order segregation error in time for velocity. Therefore, the accuracy order is consistent with second order time integrators in OpenFOAM such as Crank-Nicolson and backward schemes.

Additional improvement was made on the linearization of the convective flux (phi). A second order extrapolation is implemented in time to replace the first order lagged phi. Hence, the overall accuracy of time stepping is increased to second order.

The Rhie-Chow interpolation is applied in a very similar way with pisoFoam. The correction term related to the time discretization ,i.e. fvc::ddtPhiCorr(...), has been found unnecessary and therefore dropped out. No checkerboard pressure oscillations are observed in any case.

It is tested for turbulent and laminar flows. We observed that the projection solver is more stable than pisoFoam and therefore higher Courant numbers are feasible.

## Publications:
[1] Önder, A., & Meyers, J. (2014). HPC realization of a controlled turbulent round jet using OpenFOAM. arXiv preprint arXiv:1406.7231.

[2] Önder, A., & Meyers, J. (2014). Modification of vortex dynamics and transport properties of transitional axisymmetric jets using zero-net-mass-flux actuation. *Physics of Fluids*, 26(7), 075103.

[3] Önder, A. (2014) Active control of turbulent axisymmetric jets using zero-net-mass-flux actuation. *PhD Thesis*. KU Leuven.

[4] Önder, A., & Meyers, J. (2016). Optimal control of a transitional jet using a continuous adjoint method. *Computers & Fluids*, 126, 12-24.


## Author
Asim Önder (asim.onder@gmail.com)
