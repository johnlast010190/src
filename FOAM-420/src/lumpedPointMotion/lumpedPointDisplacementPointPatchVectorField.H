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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::lumpedPointDisplacementPointPatchVectorField

Description

    Interpolates pre-specified motion with motion specified as
    pointVectorFields.

    This is the point-patch responsible for managing the force
    integration on a 'lumped-point' basis, waiting for the external
    application, reading back the response from the external program
    and updating the locations of the associated patch points
    accordingly.

    The internal patch type name is 'lumpedPointDisplacement'.

    \heading Patch usage

    Example:
    \verbatim
    walls
    {
        type                lumpedPointDisplacement;
        value               uniform (0 0 0);
        fieldName           wantedDisplacement;
        interpolationScheme linear;
    }
    \endverbatim

    This will scan the case for \a wantedDisplacement pointVectorFields
    and interpolate those in time (using \c linear interpolation) to
    obtain the current displacement.
    The advantage of specifying displacement in this way is that it
    automatically works through decomposePar.

SourceFiles
    lumpedPointDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointDisplacementPointPatchVectorField_H
#define lumpedPointDisplacementPointPatchVectorField_H

#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "lumpedPointMovement.H"
#include "lumpedPointState.H"
#include "lumpedPointIOMovement.H"
#include "primitives/ints/lists/labelList.H"
#include "memory/tmp/tmp.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

/*---------------------------------------------------------------------------*\
        Class lumpedPointDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

    //- Convenience typedefs
    typedef lumpedPointDisplacementPointPatchVectorField patchType;
    typedef DimensionedField<vector, pointMesh> fieldType;

    // Private Member Functions

protected:

    //- The starting locations (obtained from the motionSolver).
    const pointField& points0() const;

    //- The auto-vivifying singleton for movement.
    const lumpedPointMovement& movement() const;


public:

    //- Runtime type information
    TypeName("lumpedPointDisplacement");


    // Constructors

        //- Construct from patch and internal field
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const dictionary& dict
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& pf,
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const pointPatchFieldMapper& mapper
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new lumpedPointDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& pf,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new lumpedPointDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    //- Destructor
    virtual ~lumpedPointDisplacementPointPatchVectorField();


    // Member functions

        //- The ids for all patches of this type
        static labelList patchIds(const pointVectorField& pvf);

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
