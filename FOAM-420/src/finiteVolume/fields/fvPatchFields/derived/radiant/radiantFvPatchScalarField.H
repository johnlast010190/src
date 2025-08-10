/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2010-2012 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::radiantFvPatchScalarField

Description
    Boundary for directional radiant fields like solar radiation

    Example of the boundary condition specification:
    \verbatim
        inlet
        {
            type            radiant;
            transmissivity  uniform 0.9;          // clear glass
            phi             Jsolar;
            inletValue uniform 500;      //W/m2
            value      uniform 500;      // placeholder
        }
    \endverbatim

SourceFiles
    radiantFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef radiantFvPatchScalarField_H
#define radiantFvPatchScalarField_H

#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchFields.H"
#include "mappedPatches/mappedPolyPatch/mappedPatchBase.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class radiantFvPatch Declaration
\*---------------------------------------------------------------------------*/

class radiantFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- transmissivity of this patch
        scalarField transmissivity_;

        //- radiation flux W/m2 for the radiation source
        scalar radiationFlux_;


public:

    //- Runtime type information
    TypeName("radiant");


    // Constructors

        //- Construct from patch and internal field
        radiantFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        radiantFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  radiantFvPatchScalarField
        //  onto a new patch
        radiantFvPatchScalarField
        (
            const radiantFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        radiantFvPatchScalarField
        (
            const radiantFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new radiantFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        radiantFvPatchScalarField
        (
            const radiantFvPatchScalarField&,
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
                new radiantFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Access

            //- radiation flux
            inline scalar& radiationFlux()
            {
                return radiationFlux_;
            }

            //- radiation flux
            inline scalarField& transmissivity()
            {
                return transmissivity_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
