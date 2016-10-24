/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "muSgsFireWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "muSgsBuoyantWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

/*
scalar muSgsFireWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar muSgsFireWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label muSgsFireWallFunctionFvPatchScalarField::maxIters_ = 10;
*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
void muSgsFireWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "muSgsFireWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muSgsFireWallFunctionFvPatchScalarField::
muSgsFireWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    //checkType();
    //read();
}


muSgsFireWallFunctionFvPatchScalarField::
muSgsFireWallFunctionFvPatchScalarField
(
    const muSgsFireWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


muSgsFireWallFunctionFvPatchScalarField::
muSgsFireWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    //checkType();
    //read();
}


muSgsFireWallFunctionFvPatchScalarField::
muSgsFireWallFunctionFvPatchScalarField
(
    const muSgsFireWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{
    //checkType();
}


muSgsFireWallFunctionFvPatchScalarField::
muSgsFireWallFunctionFvPatchScalarField
(
    const muSgsFireWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{
    //checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void muSgsFireWallFunctionFvPatchScalarField::read()
{

    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");
    beta_ = readScalar(sgs.thermo().lookup("beta"));
    const IOdictionary& environmentalProperties =
        db().lookupObject<IOdictionary>
        (
            "environmentalProperties"
        );
    dimensionedVector g(environmentalProperties.lookup("g"));
    magG_ = mag(g).value();
}
*/

void muSgsFireWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{

    // Get info from the SGS model
    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    // Field data
    const label patchI = patch().index();

    const scalarField& muw = sgs.mu().boundaryField()[patchI];
    scalarField& muSgsw = *this;

    const scalarField& alphaw = sgs.alpha().boundaryField()[patchI];
    const scalarField& alphaSgsw = sgs.alphaSgs()().boundaryField()[patchI];

    muSgsw = alphaSgsw*(muw/alphaw);
    Info<<"muSgsFireWallFunction!"<<endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    muSgsFireWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
