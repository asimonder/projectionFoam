/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM 
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2023 Asim Onder
     \\/     M  anipulation |
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
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow,"
        " using the PROJECTION algorithm."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
        #include "CourantNo.H"

        // Non-iterative projection scheme
        {
	  // Momentum step
          #include "UEqn.H"
	  
	  // Projection step
	  #include "pEqn.H"
	}       
	
        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
	
    return 0;
}


// ************************************************************************* //
