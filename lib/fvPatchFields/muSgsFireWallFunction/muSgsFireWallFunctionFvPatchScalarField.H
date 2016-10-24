/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Class
    muSgsFireWallFunctionFvPatchScalarField

Description
    Buoyant Thermal wall function for turbulent thermal diffusivity based on
    "Asymptotic analysis of the near-wall region of turbulence natural
    convection flows"
    Holling and Herwig. J. Fluid Mech (2005), vol 541, pp 383-397

SourceFiles
    muSgsJayatillekeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef muSgsFireWallFunctionFvPatchScalarField_H
#define muSgsFireWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
     Class muSgsFireWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class muSgsFireWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

    // Private member functions

        //- Check the type of the patch
        //void checkType();

        //- Read
        //void read();
public:

    //- Runtime type information
    TypeName("muSgsFireWallFunction");


    // Constructors

        //- Construct from patch and internal field
        muSgsFireWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        muSgsFireWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //  muSgsFireWallFunctionFvPatchScalarField
        //  onto a new patch
        muSgsFireWallFunctionFvPatchScalarField
        (
            const muSgsFireWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        muSgsFireWallFunctionFvPatchScalarField
        (
            const muSgsFireWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new muSgsFireWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        muSgsFireWallFunctionFvPatchScalarField
        (
            const muSgsFireWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new muSgsFireWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::blocking
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
