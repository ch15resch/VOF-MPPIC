{
   volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
         + (phicForces)
        )
    );

    // Update the pressure BCs to ensure flux consistency
     MRF.makeRelative(phiHbyA);     //For 
    adjustPhi(phiHbyA, U, p_rgh);
  
    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
   
    phiHbyA += phig;
    constrainPressure(p_rgh, U, phiHbyA, rAUf);
    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(alphacf*rAUf, p_rgh)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );
        Info<< " Done constraint pressure"  << endl;
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux()/alphacf;
            p_rgh.relax();
            U = HbyA + rAU*fvc::reconstruct((phig + (phicForces) - p_rghEqn.flux()/alphacf)/rAUf);
            U.correctBoundaryConditions();
        }
    }
}

#include "continuityErrs.H"
