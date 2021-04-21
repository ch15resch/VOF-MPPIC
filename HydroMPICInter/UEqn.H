    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U) + cloudSU
    );

    UEqn.relax();

    
    fvOptions.constrain(UEqn);
    
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    surfaceScalarField phicForces
    (
       fvc::flux((rAU*cloudVolSUSu)) + rAUf*fvc::interpolate(kinematicCloud.rhoEff()-rho)*(g & mesh.Sf())
    );

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
           
             fvc::reconstruct
                (
                  phicForces/rAUf + 
                  (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                  ) * mesh.magSf()
            
                )
        );

        fvOptions.correct(U);
    }
