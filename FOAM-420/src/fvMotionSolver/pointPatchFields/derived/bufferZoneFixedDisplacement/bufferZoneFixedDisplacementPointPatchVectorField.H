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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::bufferZoneFixedDisplacementPointPatchVectorField

Description
    Defines a rigid body rotation within a buffer zone a specified
    distance form the wall

SourceFiles
    bufferZoneFixedDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bufferZoneFixedDisplacementPointPatchVectorField_H
#define bufferZoneFixedDisplacementPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bufferZoneFixedDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bufferZoneFixedDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        vector radialVelocity_;
        vector linearVelocity_;
        vector CofG_;
        scalar bufferSize_;

        pointField p0_;
        labelList bufferCells_;
        labelList bufferPoints_;

public:

    //- Runtime type information
    TypeName("bufferZoneFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const bufferZoneFixedDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new bufferZoneFixedDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        bufferZoneFixedDisplacementPointPatchVectorField
        (
            const bufferZoneFixedDisplacementPointPatchVectorField&,
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
                new bufferZoneFixedDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- set cells to constrain
            virtual void manipulateMatrix
            (
                const Field<vector>&,
                labelList&,
                Field<vector>&
            );

            //- set cells to constrain
            virtual void setField
            (
                const Field<vector>&,
                Field<vector>&
            );

            //- Return cells inside buffer region
           void calculateBufferCellsAndPoints();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
