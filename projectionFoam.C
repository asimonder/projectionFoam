/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM 
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Asim Onder
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    projectionFoam

Description
    Non-iterative transient solver for incompressible flow based on an incremental
    projection scheme. The accuracy of the flux term is increased using a second
    order extrapolation in time.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    References
    \verbatim
    "HPC realization of a controlled turbulent jet using OpenFOAM"
    Asim Onder,
    Johan Meyers
    Open Source CFD International Conference 2013
    \verbatim

    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
        #include "readPISOControls.H"
        #include "CourantNo.H"

        // Non-iterative projection scheme
        {
	    //linearization of convective flux using second order extrapolation  	  
	    surfaceScalarField phi_o(0.5*(3.0*phi.oldTime()-phi.oldTime().oldTime()));
	  
	    // Momentum step

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi_o, U)
              + turbulence->divDevReff(U)
            );

	    solve(UEqn == -fvc::grad(p));
 
            // Projection step

	    dimensionedScalar dt=runTime.deltaT();
	    scalar rDeltaT=1.0/dt.value();

	    U += dt*fvc::grad(p);
	    U.correctBoundaryConditions();
      
	    phi= (fvc::interpolate(U) & mesh.Sf());
      
	    adjustPhi(phi, U, p);
	  
	    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
	      { 
		fvScalarMatrix pEqn
		  (
		   fvm::laplacian(dt, p) == fvc::div(phi)
		   );
	  
		pEqn.setReference(pRefCell, pRefValue);
	  
		if(nonOrth == nNonOrthCorr)
		  {
		    pEqn.solve(mesh.solver("pFinal"));
		  }
		else
		  {
		    pEqn.solve();
		  }
	  
		if (nonOrth == nNonOrthCorr)
		  {
		    phi -= pEqn.flux();
		  }  
	      }

            #include "continuityErrs.H"

	    U -= dt*fvc::grad(p);
	    U.correctBoundaryConditions();
	}       
	
	turbulence->correct();

	runTime.write();
	    
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
    }

    Info<< "End\n" << endl;
	
    return 0;
}


// ************************************************************************* //
