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
    (c) 2019 Esi Ltd.

Class
    Foam::timeVaryingMappedTotalPressureFvPatchScalarField

Group

Description
    This boundary condition is a mix of totalPressure and
    timeVaryingMappedFixedValue boundary condition.
    The user specifies p0 in points and these are mapped on the patch.
    The resulted field is used in the following four cases to compute p_p.

    Four
    variants are possible:

    1. incompressible subsonic:
        \f[
            p_p = p_0 - 0.5 |U|^2
        \f]
        where
        \vartable
            p_p     | incompressible pressure at patch [m2/s2]
            p_0     | incompressible total pressure [m2/s2]
            U       | velocity
        \endvartable

    2. compressible subsonic:
        \f[
            p_p = p_0 - 0.5 \rho |U|^2
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | total pressure [Pa]
            \rho    | density [kg/m3]
            U       | velocity
        \endvartable

    3. compressible transonic (\f$\gamma = 1\f$):
        \f[
            p_p = \frac{p_0}{1 + 0.5 \psi |U|^2}
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | total pressure [Pa]
            G       | coefficient given by \f$\frac{\gamma}{1-\gamma}\f$
        \endvartable

    4. compressible supersonic (\f$\gamma > 1\f$):
        \f[
            p_p = \frac{p_0}{(1 + 0.5 \psi G |U|^2)^{\frac{1}{G}}}
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | total pressure [Pa]
            \gamma  | ratio of specific heats (Cp/Cv)
            \psi    | compressibility [m2/s2]
            G       | coefficient given by \f$\frac{\gamma}{1-\gamma}\f$
        \endvartable

    The modes of operation are set by the dimensions of the pressure field
    to which this boundary condition is applied, the \c psi entry and the value
    of \c gamma:
    \table
        Mode                    | dimensions | psi   | gamma
        incompressible subsonic | p/rho      |       |
        compressible subsonic   | p          | none  |
        compressible transonic  | p          | psi   | 1
        compressible supersonic | p          | psi   | > 1
    \endtable


Usage
    \table
        Property     | Description                | Required | Default value
        U            | Velocity field name        | no       | U
        phi          | Flux field name            | no       | phi
        rho          | Density field name         | no       | rho
        psi          | Compressibility field name | no       | none
        gamma        | (Cp/Cv)                    | no       | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            timeVaryingMappedTotalPressure;
        perturb         0;
        fieldTable      p0;
        mapMethod       nearest;
        value           uniform 100000;
        rho             rho;
        psi             none;
        gamma           1;
    }
    \endverbatim

See also
    Foam::timeVaryingMappedTotalPressureFvPatchScalarField
    Foam::totalPressureFvPatchScalarField

SourceFiles
    timeVaryingMappedTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedTotalPressureFvPatchScalarField_H
#define timeVaryingMappedTotalPressureFvPatchScalarField_H

#include "fields/fvPatchFields/derived/timeVaryingMappedFixedValue/timeVaryingMappedFixedValueFvPatchFields.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class timeVaryingMappedTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMappedTotalPressureFvPatchScalarField
:
    public timeVaryingMappedFixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Static head calculator class
        staticHead staticHead_;

        //- free stream velocity to offset dynamic head
        autoPtr<vector> freeStreamVelocity_;

public:

    //- Runtime type information
    TypeName("timeVaryingMappedTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedTotalPressureFvPatchScalarField
        //  onto a new patch
        timeVaryingMappedTotalPressureFvPatchScalarField
        (
            const timeVaryingMappedTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedTotalPressureFvPatchScalarField
        (
            const timeVaryingMappedTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingMappedTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedTotalPressureFvPatchScalarField
        (
            const timeVaryingMappedTotalPressureFvPatchScalarField&,
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
                new timeVaryingMappedTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch total pressure and velocity fields
            virtual void updateCoeffs
            (
                const scalarField& p0p,
                const vectorField& Up
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
