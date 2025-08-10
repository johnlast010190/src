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
    Foam::bufferZoneFixedVelocityPointPatchVectorField

Description
    Defines a rigid body rotation within a buffer zone a specified
    distance form the wall

SourceFiles
    bufferZoneFixedVelocityPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bufferZoneFixedVelocityPointPatchVectorField_H
#define bufferZoneFixedVelocityPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bufferZoneFixedVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bufferZoneFixedVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        autoPtr<Function1<vector>> radialVelocity_;
        autoPtr<Function1<vector>> linearVelocity_;
        vector CofG_;
        scalar bufferSize_;

        labelList bufferCells_;
        labelList bufferPoints_;

public:

    //- Runtime type information
    TypeName("bufferZoneFixedVelocity");


    // Constructors

        //- Construct from patch and internal field
        bufferZoneFixedVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bufferZoneFixedVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        bufferZoneFixedVelocityPointPatchVectorField
        (
            const bufferZoneFixedVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new bufferZoneFixedVelocityPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        bufferZoneFixedVelocityPointPatchVectorField
        (
            const bufferZoneFixedVelocityPointPatchVectorField&,
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
                new bufferZoneFixedVelocityPointPatchVectorField
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
                labelList&,
                Field<vector>&
            );

            //- set cells to constrain
            virtual void setField
            (
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
