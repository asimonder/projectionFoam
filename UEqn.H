MRF.correctBoundaryVelocity(U);

//linearization of convective flux using second order extrapolation  	  
phi=0.5*(3.0*phi-phi.oldTime());

fvVectorMatrix UEqn
(
 //fvm::ddt(U)+fvm::div(phi_o, U)
 fvm::ddt(U)+fvm::div(phi, U)
 + MRF.DDt(U)
 + turbulence->divDevReff(U)
 ==
 fvOptions(U)
 );

fvOptions.constrain(UEqn);

solve(UEqn == -fvc::grad(p));
fvOptions.correct(U);
