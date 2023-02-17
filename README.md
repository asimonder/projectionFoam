# 2023 (in progress)
Updated the solver for OpenFOAM v2006. 

Benchmarks are in the works.

# 2013
The pisoFoam solver to solve incompressible transient flows is iterative. Here it is replaced with an incremental projection scheme which is non-iterative therefore costs less than pisoFoam. In this classical scheme we first solve the momentum equation using the old time step pressure. Then, unsoleinoidal part of the velocity is projected out by solving the Poisson equation. The scheme is mass-conservative and has an additional second order segregation error in time for velocity. Therefore, the accuracy order is consistent with second order time integrators in OpenFOAM such as Crank-Nicolson and backwards schemes.

Additional improvement was made on the linearization of the convective flux (phi). A second order extrapolation is implemented in time to replace the first order lagged phi. Hence, the overall accuracy of time stepping is increased to second order.

The Rhie-Chow interpolation is applied in a very similar way with pisoFoam. The correction term related to the time discretization ,i.e. fvc::ddtPhiCorr(...), has been found unnecessary and therefore dropped out. No checkerboard pressure oscillations are observed in any case.

It is tested for turbulent and laminar flows. We observed that the projection solver is more stable than pisoFoam and therefore higher Courant numbers are feasible.

References:
Onder, A., Meyers, J. (2013). HPC realization of a controlled turbulent jet using OpenFOAM. Open Source CFD International Conference 2013. Hamburg, 24-25 October 2013. (will be available online soon)
