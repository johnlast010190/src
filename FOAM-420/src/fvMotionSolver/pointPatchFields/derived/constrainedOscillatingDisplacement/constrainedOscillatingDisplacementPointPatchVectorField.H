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
    Foam::constrainedOscillatingDisplacementPointPatchVectorField

Description
    Foam::constrainedOscillatingDisplacementPointPatchVectorField

    Oscillating surface constrained to a spherical subsection bounded
    by an annular buffer zone with cubic smoothing interface

SourceFiles
    constrainedOscillatingDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef constrainedOscillatingDisplacementPointPatchVectorField_H
#define constrainedOscillatingDisplacementPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class constrainedOscillatingDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class constrainedOscillatingDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        scalar amplitude_;
        scalar frequency_;
        vector axis_;
        vector origin_;

        //main motion disk
        scalar sourceRadius_;

        //transition zone
        scalar bufferWidth_;

        //- transition polynomial
        Polynomial<4> cubicF_;


public:

    //- Runtime type information
    TypeName("constrainedOscillatingDisplacement");


    // Constructors

        //- Construct from patch and internal field
        constrainedOscillatingDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        constrainedOscillatingDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        constrainedOscillatingDisplacementPointPatchVectorField
        (
            const constrainedOscillatingDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new constrainedOscillatingDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        constrainedOscillatingDisplacementPointPatchVectorField
        (
            const constrainedOscillatingDisplacementPointPatchVectorField&,
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
                new constrainedOscillatingDisplacementPointPatchVectorField
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
