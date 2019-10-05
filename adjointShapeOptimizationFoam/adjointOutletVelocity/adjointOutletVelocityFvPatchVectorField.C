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

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "turbulentTransportModel.H"
#include "RASModel.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocityFvPatchVectorField::
adjointOutletVelocityFvPatchVectorField
(
    const adjointOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Uap =
        patch().lookupPatchField<surfaceScalarField, scalar>("Ua");
 const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");
 const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

 //   scalarField Un(mag(patch().nf() & Up));
   // vectorField UtHat((Up - patch().nf()*Un)/(Un + SMALL));

   // vectorField Uan(patch().nf()*(patch().nf() & patchInternalField()));

    const incompressible::turbulenceModel & turb = db().lookupObject<incompressible::turbulenceModel>("turbulenceProperties");
    
    scalarField nueff = turb.nuEff()().boundaryField()[patch().index()];
    const scalarField & deltainv = patch().deltaCoeffs();
    vectorField Up_n = (Up & patch().nf())*patch().nf();

    scalarField Up_ns = phip/patch().magSf();
  
    vectorField Up_t = Up - (phip * patch().Sf())/(patch().magSf()*patch().magSf());
 
    vectorField Uaneigh = Uap.patchInternalField();

    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf();

    vectorField Uaneigh_t = Uaneigh - Uaneigh_n;

    vectorField Uap_t = (nueff*deltainv*Uaneigh_t + Up_ns*Up_t   )/(Up_ns + nueff*deltainv);

    vectorField Uap_n = (phiap * patch().Sf())/(patch().magSf()*patch().magSf());

    vectorField::operator= (Uap_t + Uap_n);
// vectorField::operator=(phiap*patch().Sf()/sqr(patch().magSf()) + UtHat);
    //vectorField::operator=(Uan + UtHat);
 
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
