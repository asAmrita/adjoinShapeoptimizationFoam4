/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    laplaceAdjointFoam

Description

\*---------------------------------------------------------------------------*/
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
// Main program:

int main(int argc, char *argv[])
{
/*#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"*/

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"
    #include "CostFunctionValue.H"
    turbulence->validate();

    // Disable solvers performance output
    lduMatrix::debug = 0;
    solverPerformance::debug = 0;

    // Cost function value
    scalar J = 0;
    scalar Jold = 0;
    scalar Jk = 0;
    scalar erroru = 0;
 scalar errory = 0;
 scalar errorp = 0;
// Compute cost function value
//#include "CostFunctionValue.H"

    std::ofstream file("cost.csv");
    file << 0 << "," << J << "," << 0 << nl;

    std::ofstream errorFile("error.csv");

        runTime.write();

    while (runTime.loop() && fabs(J - Jold) > tol)
    {
        // save old cost value
        Jold = J;
     //#include "adjointShapeOptimizationFoam.H"
        alphak = alpha;

        // calculate current cost
		#include "CostFunctionValue.H"
        Jk = J;


        bool gammaFound = false;

        // calculate derivative^2 integrate((lambda*u + beta*p)^2 dv). Why??
        scalar phip0 = gSum(volField * Foam::pow(lambda * alphak.internalField() + beta * Ua.internalField(), 2));
        scalar phip1 = gSum(volField * Foam::pow(lambda * alphak.internalField() + beta * pa.internalField(), 2));

        while ((!gammaFound) && (gamma > tol))
        {
            alpha = alphak - gamma * (lambda * uk + beta * p);

            // truncate u for constrained control set
            forAll(alpha, i)
            {
                alpha[i] = min(alpha[i], alphaMax[i]);
                alpha[i] = max(alpha[i], alphaMin[i]);
            }

            alpha.correctBoundaryConditions();

            // get new y
            //solve(fvm::laplacian(k, y) -fvm::Sp(1.0, y) + beta * u + f);

            // get new cost
			#include "CostFunctionValue.H"

            if (J <= Jk - c1 * gamma * phip0)
            {
                Info << "alpha found, alpha = " << gamma << ", J = " << J << ", phip0" << phip0 << endl;

                gammaFound = true;
            }
            else
            {
                Info << "gamma NOT found, gamma = " << gamma << endl;
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

       /* // calculate the L2 norm of the error in control variables
        erroru = Foam::sqrt(gSum(volField * (Foam::pow(u.internalField() - ud.internalField(), 2))));
         errory = Foam::sqrt(gSum(volField * (Foam::pow(y.internalField() - yd.internalField(), 2))));
         errorp = Foam::sqrt(gSum(volField * (Foam::pow(p.internalField() - pd.internalField(), 2))));
        
        errorFile.open("error.csv",std::ios::app);
        errorFile << runTime.value() << "," << erroru << nl;

        errorFile << runTime.value() << "," << errory << nl;
        errorFile << runTime.value() << "," << errorp << nl;
        errorFile.close();

        Info << "Iteration no. " << runTime.timeName() << " - " << "error u " << erroru << nl;
         Info << "Iteration no. " << runTime.timeName() << " - " << "error y " << errory << nl;
         Info << "Iteration no. " << runTime.timeName() << " - " << "error p " << errorp << nl;
       /* uDiff = u - ud;
        forAll(uDiff,i)
        {
            uDiff[i] = fabs(u[i] - ud[i]);
        }
      */
  /*  uDiff=mag(u-ud);
    pDiff=mag(p-pd);
    yDiff=mag(y-yd);*/
        runTime.write();
    }

   /* file.close();
    errorFile.close();

    runTime++;
    u.write();
    p.write();
    u.write();
    uDiff.write();
    pDiff.write();
    yDiff.write();
    ud.write();
    yd.write();
    //    uc.write();
    //    udiff.write();*/

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << endl;
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "\nEnd\n"
         << endl;
    return 0;
}

// ************************************************************************* //
