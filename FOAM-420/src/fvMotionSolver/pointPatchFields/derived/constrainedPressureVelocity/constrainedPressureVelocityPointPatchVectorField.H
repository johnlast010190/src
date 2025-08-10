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
    (c) 2013 Esi Ltd.

Class
    Foam::constrainedPressureVelocityPointPatchVectorField

Description
    Foam::constrainedPressureVelocityPointPatchVectorField

    Moving surface based on acoustic pressure input
    constrained to a spherical subsection bounded
    by an annular buffer zone with cubic smoothing interface

SourceFiles
    constrainedPressureVelocityPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef constrainedPressureVelocityPointPatchVectorField_H
#define constrainedPressureVelocityPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "primitives/functions/Polynomial/Polynomial.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class constrainedPressureVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class constrainedPressureVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

         autoPtr<Function1<scalar>> pressure_;

        scalar rho0_;
        scalar gamma_;
        scalar p0_;


        vector axis_;
        vector origin_;

        //main motion disk
        scalar sourceRadius_;

        //transition zone
        scalar bufferWidth_;

        //- transition polynomial
        Polynomial<4> cubicF_;


        //- calculate gradient of DataEntry
        scalar lineGrad
        (
            const Function1<scalar>&,
            scalar x,
            scalar dx
        );

public:

    //- Runtime type information
    TypeName("constrainedPressureVelocity");


    // Constructors

        //- Construct from patch and internal field
        constrainedPressureVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        constrainedPressureVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        constrainedPressureVelocityPointPatchVectorField
        (
            const constrainedPressureVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new constrainedPressureVelocityPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        constrainedPressureVelocityPointPatchVectorField
        (
            const constrainedPressureVelocityPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new constrainedPressureVelocityPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

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
