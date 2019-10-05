/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    ajointShapeOptimizationFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas the
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/
#include "volFields.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 // Cost function value
      scalar changeInControl = 10.0;
    scalar J = 0;
    scalar Jold = 0;
    scalar Jk = 0;
    //scalar volField = 0;
   // scalar alpha = 0;
    // Compute cost function value
    #include "costFunctionValue.H"

    std::ofstream file("results.csv");
    file << 0 << "," << J << nl;
    file.close();

    Info << "\nStarting time loop\n"
         << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop()&&  fabs(J - Jold) > tol)
    {
       alpha.storePrevIter();
     
        Info<< "Time = " << runTime.timeName() << nl << endl;
J=Jold;
        zeroCells(alpha, inletCells);
 //alpha +=
           // mesh.fieldRelaxationFactor("alpha")
          // *(min(max(alpha - lambda*(Ua & U), zeroAlpha), alphaMax) - alpha);


   
       
          
            #include "stateEquation.H"
            #include "adjointEquation.H"
            

   
     #include "costFunctionValue.H"
        Jk = J;


       alphak = alpha;
        bool gammaFound = false;

        // calculate derivative^2 integrate((lambda*u + beta*p)^2 dv). Why??
 //scalar phip0 = gSum(volField * (Ua.internalField()&U.internalField()));

scalar phip0 = gSum( volField *  Foam::pow(Ua.internalField() & U.internalField() ,2));

dimensionedScalar gd = dimensionedScalar("gd", dimless * dimTime / sqr(dimLength), 1.0);

  while ((!gammaFound) && (gamma > tol))
        {
            alpha = alpha - gamma * gd * (Ua & U);

            // truncate u for constrained control set
            forAll(alpha, i)
            {
                alpha[i] = min(alpha[i], alphaMax[i]);
                alpha[i] = max(alpha[i], alphaMin[i]);
            }
            alpha.correctBoundaryConditions();

            // get new y
            //solve(fvm::laplacian(k, y) -fvm::Sp(1.0, y) + beta * u + f);
        //  for (int kk = 0; kk < 100; kk++)
          //{
            //U.storePrevIter();
           // #include "stateEquation.H"
         
           

          
            // get new cost
			#include "costFunctionValue.H"

            if (J <= Jk - c1 * gamma * phip0)
            {
                Info << "alpha found, alpha = " << gamma << ", J = " << J << ", phip0" << phip0 << endl;

                gammaFound = true;
            }
            else
            {
                Info << "alpha NOT found, alpha = " << gamma << endl;
                gamma = c2 * gamma;
            }
            Info<<J<<endl;
        }

Info << "Iteration no. " << runTime.timeName() << " - "
             << "Cost value " << J
             << " - "
             << "Cost variation" << fabs(J - Jold) << endl;

        file.open("cost.csv",std::ios::app);
        file << runTime.value() << "," << J << nl;
        file.close();

           
 file.open("results.csv", std::ios::app);
        file << runTime.value() << "," << J << nl;
        file.close();

       // laminarTransport.correct();
       // turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);

     
       // changeInControl = gSum(Foam::pow(mag(alpha - alphak), 1) * volField) / (gSum(Foam::pow(mag(alphak), 1) * volField) + SMALL);
        //Info << "change in control = " << changeInControl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
