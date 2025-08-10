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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::totalPressurePressureFvPatchScalarField

Description

    This boundary condition provides a total pressure condition.  Four
    variants are possible:

    1. p/rho formulation [m2/s2]
        \f[
            p_p = p_0 - 0.5 |U|^2
        \f]
        where
        \vartable
            p_p     | incompressible pressure at patch [m2/s2]
            p_0     | incompressible external static pressure [m2/s2]
            U       | velocity
        \endvartable

    2. p formulation [pa]
        \f[
            p_p = p_0 - 0.5 \rho |U|^2
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | external static pressure [Pa]
            \rho    | density [kg/m3]
            U       | velocity
        \endvartable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            totalPressurePressure;
        p0              uniform 1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    totalPressurePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressurePressureFvPatchScalarField_H
#define totalPressurePressureFvPatchScalarField_H

#include "derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class totalPressurePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalPressurePressureFvPatchScalarField
:
    public basePressureFvPatchScalarField
{

    // Private data

        //- Static pressure head
        scalarField pStaticHead_;

        //- Static head calculator class
        staticHead staticHead_;

        //- free stream velocity to offset dynamic head
        autoPtr<vector> freeStreamVelocity_;


    // Private Functions

        virtual void computeOtherSources();

public:

    //- Runtime type information
    TypeName("totalPressurePressure");


    // Constructors

        //- Construct from patch and internal field
        totalPressurePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalPressurePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalPressurePressureFvPatchScalarField
        //  onto a new patch
        totalPressurePressureFvPatchScalarField
        (
            const totalPressurePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalPressurePressureFvPatchScalarField
        (
            const totalPressurePressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new totalPressurePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalPressurePressureFvPatchScalarField
        (
            const totalPressurePressureFvPatchScalarField&,
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
                new totalPressurePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

        // Coefficients - these are overidden by the derived BCs

            //- C0 Field
            virtual tmp<scalarField> C0Field() const;

            //- C1 Field
            virtual tmp<scalarField> C1Field() const;

            //- C2 Field
            virtual tmp<scalarField> C2Field() const;

            //- C3 Field
            virtual tmp<scalarField> C3Field() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
