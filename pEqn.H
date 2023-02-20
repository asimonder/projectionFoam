dimensionedScalar dt=runTime.deltaT();

U += dt*fvc::grad(p);
U.correctBoundaryConditions();

phi= (fvc::interpolate(U) & mesh.Sf());
      
MRF.makeRelative(phi);

adjustPhi(phi, U, p);

while (piso.correctNonOrthogonal())
  {
    // Pressure projection
    
    fvScalarMatrix pEqn
      (
       fvm::laplacian(dt, p) == fvc::div(phi)
       );
    
    pEqn.setReference(pRefCell, pRefValue);
    
    pEqn.solve(mesh.solver(p.select(piso.finalNonOrthogonalIter())));
    
    if (piso.finalNonOrthogonalIter())
      {
	phi -= pEqn.flux();
      }
  }
	  
#include "continuityErrs.H"

U -= dt*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
