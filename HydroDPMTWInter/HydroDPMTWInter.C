/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#define MPPIC
#include "basicKinematicCollidingCloud.H"
#define basicKinematicTypeCloud basicKinematicCollidingCloud




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
     argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    const scalar Z = 0.0;
    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
       
            
            volScalarField alphapOld = 0 + kinematicCloud.theta();
            alphapOld = max(alphapOld,0.63);
            mu = rho*mixture.nu();
            
            Info<< "Evolving " << kinematicCloud.name() << endl;
            
            kinematicCloud.evolve();
            
            volScalarField alphap = 0 + kinematicCloud.theta();
            alphap = max(alphap,0.63);
            
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

           mixture.correct();
            
            
            // Update continuous phase volume fraction field
            alphac = max((alpha1+alpha2), alphacMin);
            alphac.correctBoundaryConditions();
            alphacf = fvc::interpolate(alphac);
            alphaPhic = alphacf*phi;
    
            fvVectorMatrix cloudSU(kinematicCloud.SU(U));
            volVectorField cloudVolSUSu
            (
                IOobject
                (
                    "cloudVolSUSu",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "0",
                    cloudSU.dimensions()/dimVolume,
                    Zero
                ),
                zeroGradientFvPatchVectorField::typeName
            );
    
            cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
            cloudVolSUSu.correctBoundaryConditions();
            cloudSU.source() = Zero;
            Info<< "Entering U Solver\n" << endl;
            #include "UEqn.H"
            Info<< "Entering P Solver\n" << endl;
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            
        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
