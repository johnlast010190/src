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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::solidBodyMotionDisplacementPointPatchVectorField

Description
    Enables the specification of a fixed value boundary condition using the
    solid body motion functions.

SourceFiles
    solidBodyMotionDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionDisplacementPointPatchVectorField_H
#define solidBodyMotionDisplacementPointPatchVectorField_H

#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class solidBodyMotionDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionDisplacementPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
    // Private data

        //- Is the motion solid body motion?
        bool isSolidBody_;

        //- The motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Pointer to moving reference frame
        coordinateFrame* coorFramePtr_;

        mutable autoPtr<pointField> localPoints0Ptr_;


public:

    //- Runtime type information
    TypeName("solidBodyMotionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new solidBodyMotionDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&,
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
                new solidBodyMotionDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return old points
            const pointField& localPoints0() const;


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
